# Algorithm to start Decon
from mantid.simpleapi import *
from mantid.api import PythonAlgorithm, AlgorithmFactory, MatrixWorkspaceProperty, FileProperty, FileAction, PropertyMode
from mantid.kernel import StringListValidator, StringMandatoryValidator, Direction, logger
from mantid import config
import math, os.path, numpy as np

class IndirectCylAbs(PythonAlgorithm):
 
    def category(self):
        return "Workflow\\MIDAS;PythonAlgorithms"

    def PyInit(self):
        self.declareProperty(name='Sample Input', defaultValue='Workspace', validator=StringListValidator(['Workspace','File']),
            doc='Sample input type')
        self.declareProperty(MatrixWorkspaceProperty('Sample Workspace', '', optional=PropertyMode.Optional, 
            direction=Direction.Input), doc="Name for the input Sample workspace.")
        self.declareProperty(FileProperty('Sample File', '', action=FileAction.OptionalLoad, extensions=["_red.nxs"]),
                             doc='File path for Sample file')
        self.declareProperty(name='Sample chemical formula', defaultValue='', doc = 'Sample chemical formula')
        self.declareProperty(name='Sample number density', defaultValue='', doc = 'Sample number density')
        self.declareProperty(name='Sample radius1', defaultValue='', doc = 'Sample radius1')
        self.declareProperty(name='Sample radius2', defaultValue='', doc = 'Sample radius2')

        self.declareProperty(name='Use Can', defaultValue=False, doc = 'Use Can')
        self.declareProperty(name='Can Input', defaultValue='Workspace', validator=StringListValidator(['Workspace','File']),
            doc='Can input type')
        self.declareProperty(MatrixWorkspaceProperty('Can Workspace', '', optional=PropertyMode.Optional, 
            direction=Direction.Input), doc="Name for the input Can workspace.")
        self.declareProperty(FileProperty('Can File', '', action=FileAction.OptionalLoad, extensions=["_red.nxs"]),
                             doc='File path for Can file')
        self.declareProperty(name='Can chemical formula', defaultValue='', doc = 'Can chemical formula')
        self.declareProperty(name='Can number density', defaultValue='', doc = 'Can number density')
        self.declareProperty(name='Can radius3', defaultValue='', doc = 'Can radius3')
        self.declareProperty(name='Can scale factor', defaultValue='1.0', doc = 'Scale factor to multiply can data')

        self.declareProperty(name='Step size', defaultValue='', doc = 'Step size')
        self.declareProperty(name='Beam width', defaultValue=2.0, doc = 'Beam width')
        self.declareProperty(name='Beam height', defaultValue=3.0, doc = 'Beam height')

        self.declareProperty(name='Verbose', defaultValue=False, doc = 'Switch Verbose Off/On')
        self.declareProperty(name='Plot', defaultValue=False, doc = 'Plot options')
        self.declareProperty(name='Save', defaultValue=False, doc = 'Switch Save result to nxs file Off/On')
 
    def PyExec(self):

        from IndirectCommon import StartTime, EndTime, getEfixed, addSampleLogs, getDefaultWorkingDirectory
        from IndirectAbsCor import FlatAbs
        from IndirectImport import import_mantidplot, is_supported_f2py_platform, import_f2py, unsupported_message
        if is_supported_f2py_platform():
            cylabs = import_f2py("cylabs")
        else:
            unsupported_message()
        mp = import_mantidplot()

        StartTime('Cylinder Absorption')
        workdir = config['defaultsave.directory']
        self._setup()
        self._waveRange()
        swaveWS = '__sam_wave'
        if self._diffraction:
            ConvertUnits(InputWorkspace=self._sam, OutputWorkspace=swaveWS, Target='Wavelength')
        else:
            ConvertUnits(InputWorkspace=self._sam, OutputWorkspace=swaveWS, Target='Wavelength',
            EMode='Indirect', EFixed=self._efixed)

        name = self._sam[:-4] + '_cyl'
#        wrk = getDefaultWorkingDirectory() + name
        wrk = 'CylAbs'
        wrk = wrk.ljust(120,' ')
        len_wrk = 0
        assWS = name + '_ass'
        SetSampleMaterial(swaveWS, ChemicalFormula=self._sam_chem, SampleNumberDensity=self._sam_density)
        sample = mtd[swaveWS].sample()
        sam_mat = sample.getMaterial()
        # total scattering x-section
        sigs = [sam_mat.totalScatterXSection()]
        # absorption x-section
        siga = [sam_mat.absorbXSection()]
        radii = [self._sam_rad1, self._sam_rad2]
        density = [self._sam_density]
        ndet = len(self._det)

        if self._usecan:
            cwaveWS = '__can_wave'
            if self._diffraction:
                ConvertUnits(InputWorkspace=self._can, OutputWorkspace=cwaveWS, Target='Wavelength')
            else:
                ConvertUnits(InputWorkspace=self._can, OutputWorkspace=cwaveWS, Target='Wavelength',
                    EMode='Indirect', EFixed=self._efixed)
            SetSampleMaterial(InputWorkspace=cwaveWS, ChemicalFormula=self._can_chem, SampleNumberDensity=self._can_density)
            can_sample = mtd[cwaveWS].sample()
            can_mat = can_sample.getMaterial()

        # total scattering x-section for can
            sigs.append(can_mat.totalScatterXSection())
            sigs.append(can_mat.totalScatterXSection())
        # absorption x-section for can
            siga.append(can_mat.absorbXSection())
            siga.append(can_mat.absorbXSection())
            radii.append(self._can_rad3)
            radii.append(self._can_rad3)
            density.append(self._can_density)
            density.append(self._can_density)
            ncan = 2
        else:
            ncan = 0
            sigs.append(0.0)
            sigs.append(0.0)
            siga.append(0.0)
            siga.append(0.0)
            radii.append(0.0)
            radii.append(0.0)
            density.append(0.0)
            density.append(0.0)

        dataA1 = []
        dataA2 = []
        dataA3 = []
        dataA4 = []

    #initially set errors to zero
        eZero = np.zeros(len(self._waves))

        for n in range(ndet):
            kill, A1, A2, A3, A4 = cylabs.cylabs(self._step_size, self._beam, ncan, radii,
                density, sigs, siga, self._det[n], self._wavelas, self._waves, n, wrk, len_wrk)
#      real step, beam(9), radii(4)
#      real density(3),sigs(3),siga(3)
#      real wavelas, angle, waves(10) 
#      integer ncan, n, len_wrk
#      character*120 wrk

            if self._verbose:
                logger.notice('Detector ' + str(n) + ' at angle : ' + str(self._det[n]) + ' * successful')

            dataA1 = np.append(dataA1, A1)
            dataA2 = np.append(dataA2, A2)
            dataA3 = np.append(dataA3, A3)
            dataA4 = np.append(dataA4, A4)

        sample_logs = {'sample_shape': 'cylinder', 'sample_filename': self._sam,
                        'sample_rad1': self._sam_rad1, 'sample_rad2': self._sam_rad2}
        dataX = self._waves * ndet

        if self._diffraction:
            v_axis_unit = 'dSpacing'
            v_axis_values = [1.0]
        else:
            v_axis_unit = 'MomentumTransfer'
            v_axis_values = self._q

    # Create the output workspaces
        assWS = name + '_ass'
        asscWS = name + '_assc'
        acscWS = name + '_acsc'
        accWS = name + '_acc'
        fname = name + '_abs'

        CreateWorkspace(OutputWorkspace=assWS, DataX=dataX, DataY=dataA1,
                    NSpec=ndet, UnitX='Wavelength',
                    VerticalAxisUnit=v_axis_unit, VerticalAxisValues=v_axis_values)
        addSampleLogs(assWS, sample_logs)

        CreateWorkspace(OutputWorkspace=asscWS, DataX=dataX, DataY=dataA2,
                    NSpec=ndet, UnitX='Wavelength',
                    VerticalAxisUnit=v_axis_unit, VerticalAxisValues=v_axis_values)
        addSampleLogs(asscWS, sample_logs)

        CreateWorkspace(OutputWorkspace=acscWS, DataX=dataX, DataY=dataA3,
                    NSpec=ndet, UnitX='Wavelength',
                    VerticalAxisUnit=v_axis_unit, VerticalAxisValues=v_axis_values)
        addSampleLogs(acscWS, sample_logs)

        CreateWorkspace(OutputWorkspace=accWS, DataX=dataX, DataY=dataA4,
                    NSpec=ndet, UnitX='Wavelength',
                    VerticalAxisUnit=v_axis_unit, VerticalAxisValues=v_axis_values)
        addSampleLogs(accWS, sample_logs)

        if self._usecan:
            workspaces = [assWS, asscWS, acscWS, accWS]
            AddSampleLog(Workspace=assWS, LogName='can_filename', LogType='String', LogText=str(self._can))
    	    AddSampleLog(Workspace=asscWS, LogName='can_filename', LogType='String', LogText=str(self._can))
    	    AddSampleLog(Workspace=acscWS, LogName='can_filename', LogType='String', LogText=str(self._can))
    	    AddSampleLog(Workspace=accWS, LogName='can_filename', LogType='String', LogText=str(self._can))
        else:
            workspaces = [assWS]

        group = assWS + ',' + asscWS + ',' + acscWS + ',' + accWS
        GroupWorkspaces(InputWorkspaces=group, OutputWorkspace=fname)

        if self._plot:
            graph1 = mp.plotSpectrum(workspaces, 0)
            graph2 = mp.plotTimeBin(workspaces, 0)
            graph2.activeLayer().setAxisTitle(mp.Layer.Bottom, 'Angle')

        if self._save:
            path = os.path.join(workdir,assWS + '.nxs')
            SaveNexusProcessed(InputWorkspace=assWS, Filename=path)
            if self._verbose:
                logger.notice('Output file created : '+path)

        EndTime('Cylinder Absorption')

    def _setup(self):
        self._verbose = self.getProperty('Verbose').value
        sInput = self.getPropertyValue('Sample Input')
        if sInput == 'Workspace':
            s_ws = self.getPropertyValue('Sample Workspace')
        else:
            s_ws = ''
        if sInput == 'File':
            s_file = self.getPropertyValue('Sample File')
        else:
            s_file = ''
        self._input = sInput
        self._path = s_file
        self._ws = s_ws
        self._getData()
        self._sam = self._name
        self._sam_chem = self.getPropertyValue('Sample chemical formula')
        self._sam_density = float(self.getPropertyValue('Sample number density'))
        self._sam_rad1 = float(self.getPropertyValue('Sample radius1'))
        self._sam_rad2 = float(self.getPropertyValue('Sample radius2'))

        self._usecan = self.getProperty('Use can').value
        if self._usecan:
            cInput = self.getPropertyValue('Can Input')
            if cInput == 'Workspace':
                c_ws = self.getPropertyValue('Can Workspace')
            else:
                c_ws = ''
            if cInput == 'File':
                c_file = self.getPropertyValue('Can File')
            else:
                c_file = ''
            self._input = cInput
            self._path = c_file
            self._ws = c_ws
            self._getData()
            self._can = self._name
            self._can_chem = self.getPropertyValue('Can chemical formula')
            self._can_density = float(self.getPropertyValue('Can number density'))
            self._can_rad3 = float(self.getPropertyValue('Can radius3'))
            self._can_scale = self.getPropertyValue('Can scale factor')

        self._step_size = float(self.getPropertyValue('Step size'))
        beam_width = float(self.getPropertyValue('Beam width'))
        beam_height = float(self.getPropertyValue('Beam height'))
        self._beam = [beam_height, 0.5 * beam_width, -0.5 * beam_width, (beam_width / 2), -(beam_width / 2), 0.0, beam_height, 0.0, beam_height]

        self._plot = self.getProperty('Plot').value
        self._save = self.getProperty('Save').value
		
    def _getData(self):   #get data
        if self._input == 'Workspace':
            inWS = self._ws
            self._name = inWS
            if self._verbose:
                logger.notice('Input from Workspace : '+inWS)
        elif self._input == 'File':
            self._getFileName()
            inWS = self._name
            LoadNexus(Filename=self._path, OutputWorkspace=inWS)
        else:
            raise ValueError('Input type not defined')

    def _getFileName(self):
        import os.path
        path = self._path
        if(os.path.isfile(path)): 
            base = os.path.basename(path)
            self._name = os.path.splitext(base)[0]
            ext = os.path.splitext(base)[1]
            if self._verbose:
                logger.notice('Input file : '+path)
        else:
            raise ValueError('Could not find file: ' + path)

    def _waveRange(self):
        from IndirectCommon import checkUnitIs, GetWSangles, getEfixed, GetThetaQ
# create a list of 10 equi-spaced wavelengths spanning the input data
        oWS = '__WaveRange'
        ExtractSingleSpectrum(InputWorkspace=self._sam, OutputWorkspace=oWS, WorkspaceIndex=0)
        self._diffraction = checkUnitIs(self._sam, 'dSpacing')

        if self._diffraction:
            self._det = GetWSangles(self._sam)
            self._efixed = 0.0
            ConvertUnits(InputWorkspace=oWS, OutputWorkspace=oWS, Target='Wavelength',
                     EMode='Elastic')
        else:
            self._det, self._q = GetThetaQ(self._sam)
            self._efixed = getEfixed(self._sam)
            ConvertUnits(InputWorkspace=oWS, OutputWorkspace=oWS, Target='Wavelength',
                     EMode='Indirect', EFixed=self._efixed)
        Xin = mtd[oWS].readX(0)
        xmin = mtd[oWS].readX(0)[0]
        xmax = mtd[oWS].readX(0)[len(Xin) - 1]
        ebin = 0.5
        nw1 = int(xmin/ebin)
        nw2 = int(xmax/ebin)+1
        w1 = nw1*ebin
        w2 = nw2*ebin
        waves = []
        nw = 10
        ebin = (w2-w1)/(nw-1)
        for l in range(0,nw):
            waves.append(w1+l*ebin)
        DeleteWorkspace(oWS)
        self._waves = waves

        if self._diffraction:
            self._wavelas = waves[int(nw / 2)]
        else:
            self._wavelas = math.sqrt(81.787/self._efixed) # elastic wavelength

        if self._verbose:
            logger.notice('Elastic lambda : ' + str(self._wavelas))
            nw = len(self._waves)
            message = 'Lambda : ' + str(nw) + ' values from ' + str(self._waves[0]) + ' to ' + str(self._waves[nw - 1])
            logger.notice(message)
            ndet = len(self._det)
            message = 'Detector angles : ' + str(ndet) + ' from ' + str(self._det[0]) + ' to ' + str(self._det[ndet - 1])
            logger.notice(message)

# Register algorithm with Mantid
AlgorithmFactory.subscribe(IndirectCylAbs)
#
