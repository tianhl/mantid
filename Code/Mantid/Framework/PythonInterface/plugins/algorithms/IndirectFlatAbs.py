from mantid.simpleapi import *
from mantid.api import PythonAlgorithm, AlgorithmFactory, MatrixWorkspaceProperty, WorkspaceGroupProperty, PropertyMode
from mantid.kernel import StringMandatoryValidator, Direction, logger
import math
import numpy as np


def fact(x_section, thickness, sec1, sec2):
    S = x_section * thickness * (sec1 - sec2)
    F = 1.0
    if S == 0.0:
        F = thickness
    else:
        S = (1 - math.exp(-S)) / S
        F = thickness * S
    return F


def calc_thickness_at_x_section(x_section, thickness, sec):
    sec1, sec2 = sec

    thickness_s1 = x_section * thickness * sec1
    thickness_s2 = x_section * thickness * sec2

    return thickness_s1, thickness_s2


def calc_flat_abs_can(ass, can_x_section, can_thickness_1, can_thickness_2, sampleSec1, sampleSec2, sec):
    assc = np.ones(ass.size)
    acsc = np.ones(ass.size)
    acc = np.ones(ass.size)

    sec1, sec2 = sec

    #vector version of fact
    vect_fact = np.vectorize(fact)
    f1 = vect_fact(can_x_section, can_thickness_1, sec1, sec2)
    f2 = vect_fact(can_x_section, can_thickness_2, sec1, sec2)

    can_thickness_1_s1, can_thickness_1_s2 = calc_thickness_at_x_section(can_x_section, can_thickness_1, sec)
    can_thickness_2_s1, can_thickness_2_s2 = calc_thickness_at_x_section(can_x_section, can_thickness_2, sec)

    if sec2 < 0.0:
        val = np.exp(-(can_thickness_1_s1 - can_thickness_1_s2))
        assc = ass * val

        acc1 = f1
        acc2 = f2 * val

        acsc1 = acc1
        acsc2 = acc2 * np.exp(-(sampleSec1 - sampleSec2))
    else:
        val = np.exp(-(can_thickness_1_s1 + can_thickness_2_s2))
        assc = ass * val

        acc1 = f1 * np.exp(-(can_thickness_1_s2 + can_thickness_2_s2))
        acc2 = f2 * val

        acsc1 = acc1 * np.exp(-sampleSec2)
        acsc2 = acc2 * np.exp(-sampleSec1)

    can_thichness = can_thickness_1 + can_thickness_2

    if(can_thichness > 0.0):
        acc = (acc1 + acc2) / can_thichness
        acsc = (acsc1 + acsc2) / can_thichness

    return assc, acsc, acc


class IndirectFlatAbs(PythonAlgorithm):

    def category(self):
        return "Workflow\\MIDAS;PythonAlgorithms"


    def PyInit(self):
        self.declareProperty(MatrixWorkspaceProperty('SampleWorkspace', '', direction=Direction.Input),
                             doc='Sample workspace.')

        self.declareProperty(name='SampleChemicalFormula', defaultValue='', validator=StringMandatoryValidator(),
                             doc='Sample chemical formula')
        self.declareProperty(name='SampleNumberDensity', defaultValue=0.1, doc='Sample number density')
        self.declareProperty(name='SampleThickness', defaultValue=0.0, doc='Sample thickness')
        self.declareProperty(name='SampleAngle', defaultValue=0.1, doc='Sample angle')

        self.declareProperty(MatrixWorkspaceProperty('CanWorkspace', '', optional=PropertyMode.Optional,
                                                     direction=Direction.Input),
                             doc='Container workspace.')

        self.declareProperty(name='CanChemicalFormula', defaultValue='', doc='Can chemical formula')
        self.declareProperty(name='CanNumberDensity', defaultValue=0.1, doc='Can number density')
        self.declareProperty(name='CanFrontThickness', defaultValue=0.0, doc='Can front thickness')
        self.declareProperty(name='CanBackThickness', defaultValue=0.0, doc='Can back thickness')
        self.declareProperty(name='CanScaleFactor', defaultValue=1.0, doc='Scale factor to multiply can data')

        self.declareProperty(name='Plot', defaultValue=False, doc='Plot options')

        self.declareProperty(WorkspaceGroupProperty('OutputWorkspace', '', direction=Direction.Output),
                             doc='Calculated correction functions.')


    def PyExec(self):
        from IndirectCommon import getEfixed, addSampleLogs

        self._setup()

        self._wave_range()

        sample_wave_ws = '__sam_wave'
        if self._diffraction:
            ConvertUnits(InputWorkspace=self._sample_ws, OutputWorkspace=sample_wave_ws, Target='Wavelength')
        else:
            ConvertUnits(InputWorkspace=self._sample_ws, OutputWorkspace=sample_wave_ws, Target='Wavelength',
                         EMode='Indirect', EFixed=self._efixed)

        SetSampleMaterial(sample_wave_ws, ChemicalFormula=self._sample_chemical_formula, SampleNumberDensity=self._sample_number_density)
        sample = mtd[sample_wave_ws].sample()
        sample_material = sample.getMaterial()
        can_material = None

        number_detectors = len(self._detector_angles)

        if self._can_ws is not None:
            can_wave_ws = '__can_wave'
            if self._diffraction:
                ConvertUnits(InputWorkspace=self._can_ws, OutputWorkspace=can_wave_ws, Target='Wavelength')
            else:
                ConvertUnits(InputWorkspace=self._can_ws, OutputWorkspace=can_wave_ws, Target='Wavelength',
                             EMode='Indirect', EFixed=self._efixed)

            SetSampleMaterial(InputWorkspace=can_wave_ws, ChemicalFormula=self._can_chemical_formula, SampleNumberDensity=self._can_number_density)
            can_sample = mtd[can_wave_ws].sample()
            can_material = can_sample.getMaterial()

        data_a1 = []
        data_a2 = []
        data_a3 = []
        data_a4 = []

        for det in range(number_detectors):
            det_angle = self._detector_angles[det]
            (a1, a2, a3, a4) = self._flat_abs(sample_material, can_material, det_angle)
            logger.information('Processed detector %d at angle %f' % (det, det_angle))

            data_a1 = np.append(data_a1, a1)
            data_a2 = np.append(data_a2, a2)
            data_a3 = np.append(data_a3, a3)
            data_a4 = np.append(data_a4, a4)

        DeleteWorkspace(sample_wave_ws)
        if self._can_ws is not None:
            DeleteWorkspace(can_wave_ws)

        sample_logs = {'sample_shape': 'flatplate',
                       'sample_filename': self._sample_ws,
                       'sample_thickness': self._sample_thickness,
                       'sample_angle': self._sample_angle}
        data_x = self._wavelengths * number_detectors

        if self._diffraction:
            v_axis_unit = 'dSpacing'
            v_axis_values = [1.0]
        else:
            v_axis_unit = 'MomentumTransfer'
            v_axis_values = self._q

        CreateWorkspace(OutputWorkspace=self._output_name('ass'),
                        DataX=data_x, DataY=data_a1,
                        NSpec=number_detectors, UnitX='Wavelength',
                        VerticalAxisUnit=v_axis_unit, VerticalAxisValues=v_axis_values)
        addSampleLogs(self._output_name('ass'), sample_logs)

        CreateWorkspace(OutputWorkspace=self._output_name('assc'),
                        DataX=data_x, DataY=data_a2,
                        NSpec=number_detectors, UnitX='Wavelength',
                        VerticalAxisUnit=v_axis_unit, VerticalAxisValues=v_axis_values)
        addSampleLogs(self._output_name('assc'), sample_logs)

        CreateWorkspace(OutputWorkspace=self._output_name('ascs'),
                        DataX=data_x, DataY=data_a3,
                        NSpec=number_detectors, UnitX='Wavelength',
                        VerticalAxisUnit=v_axis_unit, VerticalAxisValues=v_axis_values)
        addSampleLogs(self._output_name('ascs'), sample_logs)

        CreateWorkspace(OutputWorkspace=self._output_name('acc'),
                        DataX=data_x, DataY=data_a4,
                        NSpec=number_detectors, UnitX='Wavelength',
                        VerticalAxisUnit=v_axis_unit, VerticalAxisValues=v_axis_values)
        addSampleLogs(self._output_name('acc'), sample_logs)

        if self._can_ws is not None:
            workspaces = [self._output_name(ws_type) for ws_type in ['ass', 'assc', 'ascs', 'acc']]
            AddSampleLog(Workspace=self._output_name('ass'), LogName='can_filename', LogType='String', LogText=str(self._can_ws))
            AddSampleLog(Workspace=self._output_name('assc'), LogName='can_filename', LogType='String', LogText=str(self._can_ws))
            AddSampleLog(Workspace=self._output_name('ascs'), LogName='can_filename', LogType='String', LogText=str(self._can_ws))
            AddSampleLog(Workspace=self._output_name('acc'), LogName='can_filename', LogType='String', LogText=str(self._can_ws))
        else:
            workspaces = [self._output_name('ass')]

        GroupWorkspaces(InputWorkspaces=workspaces, OutputWorkspace=self._output_ws)
        self.setProperty('OutputWorkspace', self._output_ws)

        if self._plot:
            from IndirectImport import import_mantidplot
            mantid_plot = import_mantidplot()
            mantid_plot.plotSpectrum(workspaces, 0)
            time_bin = mantid_plot.plotTimeBin(workspaces, 0)
            time_bin.activeLayer().setAxisTitle(mantid_plot.Layer.Bottom, 'Angle')


    def _setup(self):
        """
        Get algorithm properties.
        """

        self._sample_ws = self.getPropertyValue('SampleWorkspace')
        self._sample_chemical_formula = self.getPropertyValue('SampleChemicalFormula')
        self._sample_number_density = self.getProperty('SampleNumberDensity').value
        self._sample_thickness = self.getProperty('SampleThickness').value
        self._sample_angle = self.getProperty('SampleAngle').value

        self._can_ws = self.getPropertyValue('CanWorkspace')
        if self._can_ws == '':
            self._can_ws = None

        self._can_chemical_formula = self.getPropertyValue('CanChemicalFormula')
        self._can_number_density = self.getProperty('CanNumberDensity').value
        self._can_front_thickness = self.getProperty('CanFrontThickness').value
        self._can_back_thickness = self.getProperty('CanBackThickness').value
        self._can_scale = self.getProperty('CanScaleFactor').value

        self._plot = self.getProperty('Plot').value
        self._output_ws = self.getPropertyValue('OutputWorkspace')

        self._output_name = lambda ws_type: self._output_ws + '_' + ws_type


    def _wave_range(self):
        """
        Create a list of 10 equi-spaced wavelengths spanning the input data.
        """

        from IndirectCommon import checkUnitIs, GetWSangles, getEfixed, GetThetaQ

        sample_spec_ws = '__WaveRange'
        ExtractSingleSpectrum(InputWorkspace=self._sample_ws, OutputWorkspace=sample_spec_ws, WorkspaceIndex=0)
        self._diffraction = checkUnitIs(self._sample_ws, 'dSpacing')

        if self._diffraction:
            self._detector_angles = GetWSangles(self._sample_ws)
            self._efixed = 0.0
            ConvertUnits(InputWorkspace=sample_spec_ws, OutputWorkspace=sample_spec_ws, Target='Wavelength',
                         EMode='Elastic')
        else:
            self._detector_angles, self._q = GetThetaQ(self._sample_ws)
            self._efixed = getEfixed(self._sample_ws)
            ConvertUnits(InputWorkspace=sample_spec_ws, OutputWorkspace=sample_spec_ws, Target='Wavelength',
                         EMode='Indirect', EFixed=self._efixed)

        data_x_in = mtd[sample_spec_ws].readX(0)
        x_min = mtd[sample_spec_ws].readX(0)[0]
        x_max = mtd[sample_spec_ws].readX(0)[len(data_x_in) - 1]
        ebin = 0.5
        nw1 = int(x_min / ebin)
        nw2 = int(x_max / ebin) + 1
        w1 = nw1 * ebin
        w2 = nw2 * ebin
        waves = []
        number_wavelengths = 10
        ebin = (w2 - w1) / (number_wavelengths - 1)
        for l in range(0, number_wavelengths):
            waves.append(w1 + l * ebin)
        DeleteWorkspace(sample_spec_ws)
        self._wavelengths = waves

        if self._diffraction:
            self._wavelas = waves[int(number_wavelengths / 2)]
        else:
            self._wavelas = math.sqrt(81.787 / self._efixed)  # Elastic wavelength

        logger.information('Elastic lambda: ' + str(self._wavelas))
        logger.information('Lambda: %d. From %d to %d' % (len(self._wavelengths), self._wavelengths[0], self._wavelengths[-1]))
        logger.information('Detector angles: %d. From %d to %d' %
                           (len(self._detector_angles), self._detector_angles[0], self._detector_angles[-1]))


    def _flat_abs(self, sample_material, can_material, detector_angle):
        """
        Calculate flat plate absorption factors

        For more information See:
          - MODES User Guide: http://www.isis.stfc.ac.uk/instruments/iris/data-analysis/modes-v3-user-guide-6962.pdf
          - C J Carlile, Rutherford Laboratory report, RL-74-103 (1974)

        @param sample_material Material sample is composed of
        @param can_material Material container is composed of
        @param detector_angle Angle of detector
        """

        PICONV = math.pi / 180.0

        # Can angle and detector angle
        can_angle = self._sample_angle * PICONV
        theta = detector_angle * PICONV

        # tsec is the angle the scattered beam makes with the normal to the sample surface.
        tsec = detector_angle - self._sample_angle

        nlam = len(self._wavelengths)

        ass = np.ones(nlam)
        assc = np.ones(nlam)
        acsc = np.ones(nlam)
        acc = np.ones(nlam)

        # Case where tsec is close to 90 degrees. CALCULATION IS UNRELIABLE
        if abs(abs(tsec) - 90.0) < 1.0:
            # Default to 1 for everything
            return ass, assc, acsc, acc

        else:
            tsec = tsec * PICONV

            sec1 = 1.0 / math.cos(can_angle)
            sec2 = 1.0 / math.cos(tsec)

            # List of wavelengths
            waves = np.array(self._wavelengths)

            # Sample cross section
            sample_x_section = (sample_material.totalScatterXSection() + sample_material.absorbXSection() * waves / 1.8) * self._sample_number_density

            # Vector version of fact
            vect_fact = np.vectorize(fact)
            fs = vect_fact(sample_x_section, self._sample_thickness, sec1, sec2)

            sampleSec1, sampleSec2 = calc_thickness_at_x_section(sample_x_section, self._sample_thickness, [sec1, sec2])

            if sec2 < 0.0:
                ass = fs / self._sample_thickness
            else:
                ass = np.exp(-sampleSec2) * fs / self._sample_thickness

            if self._can_ws is not None:
                # Calculate can cross section
                can_x_section = (can_material.totalScatterXSection() + can_material.absorbXSection() * waves / 1.8) * self._can_number_density
                assc, acsc, acc = calc_flat_abs_can(ass, can_x_section,
                                                    self._can_front_thickness, self._can_back_thickness,
                                                    sampleSec1, sampleSec2,
                                                    [sec1, sec2])

        return ass, assc, acsc, acc


# Register algorithm with Mantid
AlgorithmFactory.subscribe(IndirectFlatAbs)
