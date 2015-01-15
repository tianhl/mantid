from mantid.simpleapi import *
from mantid.api import PythonAlgorithm, AlgorithmFactory, MatrixWorkspaceProperty, WorkspaceGroupProperty, PropertyMode
from mantid.kernel import StringMandatoryValidator, Direction, logger
import math
import numpy as np


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
        from IndirectAbsCor import FlatAbs

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
        sam_mat = sample.getMaterial()
        # Total scattering x-section
        sigs = [sam_mat.totalScatterXSection()]
        # Absorption x-section
        siga = [sam_mat.absorbXSection()]
        size = [self._sample_thickness]
        density = [self._sample_number_density]
        ncan = 0
        ndet = len(self._detector_angles)

        if self._can_ws is not None:
            can_wave_ws = '__can_wave'
            if self._diffraction:
                ConvertUnits(InputWorkspace=self._can_ws, OutputWorkspace=can_wave_ws, Target='Wavelength')
            else:
                ConvertUnits(InputWorkspace=self._can_ws, OutputWorkspace=can_wave_ws, Target='Wavelength',
                             EMode='Indirect', EFixed=self._efixed)

            SetSampleMaterial(InputWorkspace=can_wave_ws, ChemicalFormula=self._can_chemical_formula, SampleNumberDensity=self._can_density)
            can_sample = mtd[can_wave_ws].sample()
            can_mat = can_sample.getMaterial()

            # Total scattering x-section for can
            sigs.append(can_mat.totalScatterXSection())
            sigs.append(can_mat.totalScatterXSection())
            # Absorption x-section for can
            siga.append(can_mat.absorbXSection())
            siga.append(can_mat.absorbXSection())
            size.append(self._can_front_thickness)
            size.append(self._can_back_thickness)
            density.append(self._can_density)
            density.append(self._can_density)
            ncan = 2

        data_a1 = []
        data_a2 = []
        data_a3 = []
        data_a4 = []

        for det in range(ndet):
            det_angle = self._detector_angles[det]
            angles = [self._sample_angle, det_angle]
            (a1, a2, a3, a4) = FlatAbs(ncan, size, density, sigs, siga, angles, self._wavelengths)
            logger.information('Processed detector %d at angle %f' % (det, det_angle))

            data_a1 = np.append(data_a1, a1)
            data_a2 = np.append(data_a2, a2)
            data_a3 = np.append(data_a3, a3)
            data_a4 = np.append(data_a4, a4)

        sample_logs = {'sample_shape': 'flatplate',
                       'sample_filename': self._sample_ws,
                       'sample_thickness': self._sample_thickness,
                       'sample_angle': self._sample_angle}
        data_x = self._wavelengths * ndet

        if self._diffraction:
            v_axis_unit = 'dSpacing'
            v_axis_values = [1.0]
        else:
            v_axis_unit = 'MomentumTransfer'
            v_axis_values = self._q

        CreateWorkspace(OutputWorkspace=self._output_name('ass'), DataX=data_x, DataY=data_a1,
                        NSpec=ndet, UnitX='Wavelength',
                        VerticalAxisUnit=v_axis_unit, VerticalAxisValues=v_axis_values)
        addSampleLogs(self._output_name('ass'), sample_logs)

        CreateWorkspace(OutputWorkspace=self._output_name('assc'), DataX=data_x, DataY=data_a2,
                        NSpec=ndet, UnitX='Wavelength',
                        VerticalAxisUnit=v_axis_unit, VerticalAxisValues=v_axis_values)
        addSampleLogs(self._output_name('assc'), sample_logs)

        CreateWorkspace(OutputWorkspace=self._output_name('ascs'), DataX=data_x, DataY=data_a3,
                        NSpec=ndet, UnitX='Wavelength',
                        VerticalAxisUnit=v_axis_unit, VerticalAxisValues=v_axis_values)
        addSampleLogs(self._output_name('ascs'), sample_logs)

        CreateWorkspace(OutputWorkspace=self._output_name('acc'), DataX=data_x, DataY=data_a4,
                        NSpec=ndet, UnitX='Wavelength',
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
        self._can_density = self.getProperty('CanNumberDensity').value
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


# Register algorithm with Mantid
AlgorithmFactory.subscribe(IndirectFlatAbs)
