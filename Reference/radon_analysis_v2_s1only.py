#!/usr/bin/env python
# -*- coding: UTF-8 -*-
from radon_analysis_v2_s1s2 import *
#########################################################################################################################
class S1onlyRadonAnalysisMain(Analysis):
    """ Manage main radon analysis steps. """
    def __init__(self, config_dict, run_list, load_func, **kwargs):
        super().__init__(**kwargs)
        #self.logging.setLevel(logging.DEBUG)
        radon_analysis_class_initialization(self, config_dict, run_list, load_func, **kwargs)
        if (self.run_list is None) and (self.config['analysis_main']['load_from_run_list'] == True):
            self.run_list = read_run_list_file(self.config['analysis_main']['run_list_file'])
        self.step_functions = {'correction_derivation':       self._correction_derivation,
                               'estimate_fiducial_mass':      self._estimate_fiducial_mass,
                               'peak_region_rate_estimation': self._peak_region_rate_estimation,
                               'surface_rate_estimation':     self._surface_rate_estimation,
                               'custom_process':              self._custom_process}
        self.data_dirs  = self.config['analysis_main']['data_dir_names']  
        self.result_dir = self.config['analysis_main']['result_dir_name']
        self.plot_dir   = self.config['general']['plot_dir']
        pathlib.Path(self.result_dir).mkdir(parents=True, exist_ok=True)
        pathlib.Path(self.plot_dir)  .mkdir(parents=True, exist_ok=True)

        self.python_module = None
        if os.path.isfile(self.config['general']['python_config_file']) == True:
            python_config_file = self.config['general']['python_config_file'].replace('.py', '')
            python_config_path = os.path.join(*python_config_file.split('/')[:-1])
            if python_config_file[0] == '/': python_config_path = '/' + python_config_path
            sys.path.append(python_config_path)
            self.python_module = importlib.import_module(python_config_file.split('/')[-1])
            ALPHA_LINE_E                    = self.python_module.ALPHA_LINE_E
            ALPHA_LINE_C                    = self.python_module.ALPHA_LINE_C
            SR_TIMESTAMP_RANGE_TEXTS        = self.python_module.SR_TIMESTAMP_RANGE_TEXTS
            RN222CALI_TIMESTAMP_RANGE_TEXTS = self.python_module.RN222CALI_TIMESTAMP_RANGE_TEXTS
        else: print('\n\nNote: general:python_config_file not set\n\n')
    def run(self, args):
        """ Run requested analysis step. """
        for analysis_step_name in self.step_functions:
            if self.config['analysis_main']['process_step'][analysis_step_name] == True: 
                if self.verbosity >= 1: print('Proceeding with step: '+analysis_step_name)
                self.step_functions[analysis_step_name](args)
    def _load_data(self, args, load_pickle_data=None, data_dirs=None, iter_range=None, rateDict_path=None):
        """ Manage data loading. """
        if load_pickle_data is None: load_pickle_data = self.config['analysis_main']['load_pickle_data'] 
        if data_dirs        is None: data_dirs        = self.data_dirs
        updateNew, preSubN = False, 0
        if iter_range is None:
            iter_range = [0, None]
            if isinstance(self.config['analysis_main']['iter_range'], list) == True:
                iter_range = self.config['analysis_main']['iter_range']
                if iter_range[0] != None: iter_range[0] = int(iter_range[0])
                if iter_range[1] != None: iter_range[1] = int(iter_range[1])
            elif 'update' in self.config['analysis_main']['iter_range'].lower(): 
                updateNew = True
                if '-' in self.config['analysis_main']['iter_range']:
                    preSubN = int(self.config['analysis_main']['iter_range'].split('-')[-1])
        run_exclude = run_numbers_to_run_names(self.config['analysis_main']['run_exclude'])        
        run_include = None 
        if self.run_list is not None: run_include = run_numbers_to_run_names(self.run_list)
        date_range      = self.config['analysis_main']['date_range']
        timestamp_range = [YYMMDDtoTimestamp(date_range[0]), YYMMDDtoTimestamp(date_range[1])]
        
        runInfoDF_selection, procDF_selection = None, None
        if hasattr(self.python_module, 'runInfoDF_selection'):
            runInfoDF_selection = self.python_module.runInfoDF_selection
        else: print('\n\nNote: runInfoDF_selection() does not exist in general:python_config_file\n\n')
        if hasattr(self.python_module, 'procDF_selection'):
            procDF_selection = self.python_module.procDF_selection
        else: print('\n\nNote: procDF_selection() does not exist in general:python_config_file\n\n')
 
        existingRunChunks = None
        if rateDict_path is not None:
            if os.path.isfile(rateDict_path) == True:
                rateDictPickle = None
                with open(rateDict_path, 'rb') as handle:
                    rateDictPickle = pickle.load(handle)
                existingRunChunks = np.array(rateDictPickle['runChunks'], dtype=list)

        if self.verbosity >= 1: 
            print('Loading data:')
            print('  load_pickle_data =', self.config['analysis_main']['load_pickle_data'])
            print('  updateNew        =', updateNew)
            print('  run_exclude      =', run_exclude)
            print('  run_include      =', run_include)
            print('  date_range       =', date_range)
            print('  timestamp_range  =', timestamp_range)
        runInfoData, procData = None, None
        if self.config['analysis_main']['load_pickle_data'] == True: 
            _, runInfoData, _, procData =loadPickleRunEventData(data_dirs,iter_range=iter_range,run_exclude=run_exclude,\
                                                                run_include=run_include,timestamp_range=timestamp_range,\
                                                                runInfoDF_selection=runInfoDF_selection,\
                                                                run_chunks=existingRunChunks,\
                                                                procDF_selection=procDF_selection,\
                                                                verbosityLevel=(self.verbosity-1)*10)
        else:
            load_context_with_run_info_s1only(self)
            if updateNew == True:
                if (existingRunChunks is not None):
                    iter_range = [-np.sum(np.array(self.run_info['number']) > int(existingRunChunks[-1][-1]))-1, None]
                    iter_range[0] -= preSubN
                    if self.verbosity >= 1: print('\n\n_load_data(): updateNew => iter_range =', iter_range, '\n\n') 
                else:print('\n\nNote: updateNew does nothing because there are no pickle file to extract runs from\n\n')
            if iter_range[0] is None:
                iter_range[0] = -np.sum( np.array(self.run_info['number']) > int(existingRunChunks[-1][-1]) ) - 1
            runInfoDF, runInfoData, procDF, procData = load_run_event_data(\
                self, iter_range=iter_range, run_exclude=run_exclude, run_include=run_include,\
                timestamp_range=timestamp_range, runInfoDF_selection=runInfoDF_selection,\
                run_chunks=existingRunChunks, verbosityLevel=(self.verbosity-1)*10)
        if self.verbosity >= 1: print('  generating user defined data via gen_procData...')
        self.python_module.gen_procData(procData)
        return runInfoData, procData
    def _correction_derivation(self, args):
        """ Perform S1_area corrections for the alphas. """
        self.config['analysis_main']['correction_derivation']['generalBinN']   = self.config['general']['generalBinN']
        runInfoData, procData = self._load_data(args)
        mask_dict = self.python_module.gen_masks(procData, verbosity=self.verbosity)
        correction_derivation(procData, mask_dict, output_dir=self.result_dir, plot_dir=self.plot_dir,\
                              save_plot=self.save_plot, plot_format=self.plot_format,\
                              output_name_tag=self.config['analysis_main']['output_name_tag'],\
                              config_dict=self.config['analysis_main']['correction_derivation'],\
                              plot_extra=self.python_module.plot_extra, verbosity=self.verbosity)
    def _estimate_fiducial_mass(self, args):
        """ Estimate the s1only fiducial volume by comparing to the s1s2 data. """
        self.config['analysis_main']['estimate_fiducial_mass']['generalBinN']   = self.config['general']['generalBinN']
        runInfoData, procData = self._load_data(args)
        _,           s1s2Data = self._load_data(args, load_pickle_data=True,\
            data_dirs =self.config['analysis_main']['estimate_fiducial_mass']['s1s2_data_dir_names'],\
            iter_range=self.config['analysis_main']['estimate_fiducial_mass']['s1s2_data_iter_range'])
        if self.verbosity >= 1: print('  merging s1only and s1s2 data...')
        runInfoData['run_number'] = np.array([str(run_id).zfill(6) for run_id in runInfoData['number']]) 
        s1only_run_ids = np.unique(procData['run_number'])
        s1s2_run_ids   = np.unique(s1s2Data['run_number'])
        if 's1_time' not in s1s2Data:
            raise AssertionError('_estimate_fiducial_mass(): is s1s2Data valid? The s1s2Data keys are:\n'+\
                                 str([key for key in s1s2Data]))

        overlap_run_ids = np.array([run_id for run_id in s1only_run_ids if run_id in s1s2_run_ids])
        for run_id in s1only_run_ids:
            if run_id not in overlap_run_ids:
                warnings.warn('_estimate_fiducial_mass(): run '+str(run_id)+' is not found in s1s2 dataset', Warning)
        for run_id in s1s2_run_ids:
            if run_id not in overlap_run_ids:
                warnings.warn('_estimate_fiducial_mass(): run '+str(run_id)+' is not found in s1only dataset', Warning)
        sel_run_ids = copy.copy(overlap_run_ids)

        if self.config['analysis_main']['estimate_fiducial_mass']['remove_s1only_runs'] == True:
            if self.verbosity >= 1: print('  loading slow control data to remove s1only runs:')
            if self.context is None: load_context_with_run_info_s1only(self)
            sel_run_ids = select_runs_anode_voltage(overlap_run_ids, runInfoData, self.context, selectLow=True,\
                                                    verbosity=self.verbosity)
            for run_id in overlap_run_ids:
                if run_id not in sel_run_ids: print('    s1only run '+str(run_id)+' is removed') 
        if self.verbosity >= 1: print('  the remaining overlapped runs are:', list(sel_run_ids))

        runMask = np.isin(runInfoData['run_number'], sel_run_ids)
        for key in runInfoData: runInfoData[key] = runInfoData[key][runMask]
        runMask = np.isin(procData['run_number'], sel_run_ids)
        for key in procData: procData[key] = procData[key][runMask]
        runMask = np.isin(s1s2Data['run_number'], sel_run_ids)
        for key in s1s2Data: s1s2Data[key] = s1s2Data[key][runMask]
        
        ### s1s2_ContainsEventsNotIn_s1only WARNING
        eventMask = np.isin(procData['time'], s1s2Data['s1_time'])
        if len(procData['time'][eventMask]) != len(s1s2Data['s1_time']):
            if self.verbosity >= 1: 
                print('WARNING: exist runs that have events in s1s2 that are not in s1only')
                print('processing to find such runs...')
            s1s2_ContainsEventsNotIn_s1only = []
            for run_id in tqdm(sel_run_ids, disable=(self.verbosity<=0)):
                runMasks1only = (procData['run_number'] == run_id)
                runMasks1s2   = (s1s2Data['run_number'] == run_id)
                eventMask = np.isin(procData['time'][runMasks1only], s1s2Data['s1_time'][runMasks1s2])
                if len(procData['time'][runMasks1only][eventMask]) != len(s1s2Data['s1_time'][runMasks1s2]):
                    s1s2_ContainsEventsNotIn_s1only.append(run_id)
                    warnings.warn('_estimate_fiducial_mass(): run '+str(run_id)+' has events in s1s2 '+\
                                  'that are not in s1only', Warning)
            sel_run_ids = [run_id for run_id in sel_run_ids if run_id not in s1s2_ContainsEventsNotIn_s1only]
            runMask = np.isin(runInfoData['run_number'], sel_run_ids)
            for key in runInfoData: runInfoData[key] = runInfoData[key][runMask]
            runMask = np.isin(procData['run_number'], sel_run_ids)
            for key in procData: procData[key] = procData[key][runMask]
            runMask = np.isin(s1s2Data['run_number'], sel_run_ids)
            for key in s1s2Data: s1s2Data[key] = s1s2Data[key][runMask]
        #######################################################################
        s1corr_paths = self.config['analysis_main']['estimate_fiducial_mass']['s1corr_paths']
        if self.verbosity >= 1: print('  loading correction_derivation() files:', s1corr_paths)
        corrParPickle = None
        for s1corrPath in s1corr_paths:
            with open(s1corrPath, 'rb') as handle:
                corrParPickle = pickle.load(handle)
        if 'r2CorrPars' in corrParPickle:
            procData['r2_corr']  = np.array(procData['r2_naive'], dtype='float64')
            procData['r2_corr'] *= pow(60.0, 2)/polyFunc(procData['area_fraction_top'], corrParPickle['r2CorrPars'])
        elif 'r2InterpInputs' in corrParPickle:
            interpFunc = scipy.interpolate.interp1d(corrParPickle['zSpan'], corrParPickle['r2InterpInputs'])
            procData['r2_corr']  = np.array(procData['r2_s1_cnn'], dtype='float64')
            procData['r2_corr'] *= corrParPickle['r2WallRef']/interpFunc(procData['z_s1_cnn'])
        procData['s1_corr2']  = np.array(procData['area'], dtype='float64')
        procData['s1_corr2'] *= corrParPickle['areaCorr2Pars'][0]
        procData['s1_corr2'] /= linearFunc(procData['r2_corr'], *corrParPickle['areaCorr2Pars'])
        procData['s1_corrF']  = np.array(procData['s1_corr2'], dtype='float64')
        procData['s1_corrF'] *= corrParPickle['areaCorr3Pars'][-1]
        procData['s1_corrF'] /= exp2Func(procData['area_fraction_top'], *corrParPickle['areaCorr3Pars'])

        altS1ExistMask = (s1s2Data['alt_s1_time'] > 0)
        eventMask    = np.isin(procData['time'], s1s2Data['s1_time'])
        eventAltMask = np.isin(procData['time'], s1s2Data['alt_s1_time'][altS1ExistMask])

        procData_S1onlyOnly = {}
        for key in procData:
            procData_S1onlyOnly[key] = copy.deepcopy(procData[key][~eventMask]) #~(eventMask | eventAltMask)])
            procData[key]            = copy.deepcopy(procData[key][eventMask])
        for key in s1s2Data:
            if key not in procData: procData[key] = s1s2Data[key]
        if self.verbosity >= 1: 
            print('    final keys for procData:', [key for key in procData])
            print('      eventMask   :', np.sum(eventMask))
            print('      eventAltMask:', np.sum(eventAltMask))
        if (len(procData['time']) != len(procData['s1_time'])) or (np.sum(procData['time'] - procData['s1_time'])!= 0):
            raise AssertionError('_estimate_fiducial_mass(): s1only/s1s2 data merge dimension mismatch')
 
        mask_dict_S1onlyOnly = self.python_module.gen_masks(procData_S1onlyOnly, verbosity=self.verbosity)
        mask_dict            = self.python_module.gen_masks(procData,            verbosity=self.verbosity)
        s1s2_mask_dict = None
        s1s2_python_config_file = self.config['analysis_main']['estimate_fiducial_mass']['s1s2_python_config_file']
        if os.path.isfile(s1s2_python_config_file) == False:
            raise AssertionError('_estimate_fiducial_mass(): s1s2_python_config_file not found:\n  '+\
                                 s1s2_python_config_file)
        s1s2_python_config_file = s1s2_python_config_file.replace('.py', '')
        python_config_path = os.path.join(*s1s2_python_config_file.split('/')[:-1])
        if s1s2_python_config_file[0] == '/': python_config_path = '/' + python_config_path
        sys.path.append(python_config_path)
        s1s2_python_module = importlib.import_module(s1s2_python_config_file.split('/')[-1])
        s1s2_mask_dict = s1s2_python_module.gen_masks(procData, verbosity=self.verbosity)
        estimate_fiducial_mass(runInfoData, procData, mask_dict, procData_S1onlyOnly, mask_dict_S1onlyOnly,s1corr_paths,\
                               s1s2_mask_dict=s1s2_mask_dict, output_dir=self.result_dir, plot_dir=self.plot_dir,\
                               save_plot=self.save_plot, plot_format=self.plot_format,\
                               output_name_tag=self.config['analysis_main']['output_name_tag'],\
                               config_dict=self.config['analysis_main']['estimate_fiducial_mass'],\
                               plot_extra=self.python_module.plot_extra, verbosity=self.verbosity)
    def _peak_region_rate_estimation(self, args):
        """ Estimate the alpha rate through the counting within peak regions. """
        self.config['analysis_main']['peak_region_rate_estimation']['deadtime_type'] = \
            self.config['general']['deadtime_type']
        self.config['analysis_main']['peak_region_rate_estimation']['generalBinN'] =self.config['general']['generalBinN']
        rateDict_path = os.path.join(self.result_dir,\
                                     self.config['analysis_main']['output_name_tag']+\
                                     's1only_peak_region_rate_estimation.pkl')
        runInfoData, procData = self._load_data(args, rateDict_path=rateDict_path)
        runInfoData['run_number'] = np.array([str(run_id).zfill(6) for run_id in runInfoData['number']])
      
        if self.config['analysis_main']['peak_region_rate_estimation']['remove_s1s2_runs'] == True:
            if self.verbosity >= 1: print('  loading slow control data to remove s1s2 runs:')
            if self.context is None: load_context_with_run_info_s1only(self)
            run_ids = np.unique(procData['run_number'])
            sel_run_ids = select_runs_anode_voltage(run_ids, runInfoData, self.context, selectLow=False,\
                                                    verbosity=self.verbosity)
            for run_id in run_ids:
                if run_id not in sel_run_ids: print('    s1s2 run '+str(run_id)+' is removed')
            runMask = np.isin(runInfoData['run_number'], sel_run_ids)
            for key in runInfoData: runInfoData[key] = runInfoData[key][runMask]
            runMask = np.isin(procData['run_number'], sel_run_ids)
            for key in procData: procData[key] = procData[key][runMask]
        
        s1corr_paths = self.config['analysis_main']['peak_region_rate_estimation']['s1corr_paths']
        if self.verbosity >= 1: print('  loading correction_derivation() files:', s1corr_paths)
        corrParPickle = None
        for s1corrPath in s1corr_paths:
            with open(s1corrPath, 'rb') as handle:
                corrParPickle = pickle.load(handle)
        if 'r2CorrPars' in corrParPickle:
            procData['r2_corr']  = np.array(procData['r2_naive'], dtype='float64')
            procData['r2_corr'] *= pow(60.0, 2)/polyFunc(procData['area_fraction_top'], corrParPickle['r2CorrPars'])
        elif 'r2InterpInputs' in corrParPickle:
            interpFunc = scipy.interpolate.interp1d(corrParPickle['zSpan'], corrParPickle['r2InterpInputs'])
            procData['r2_corr']  = np.array(procData['r2_s1_cnn'], dtype='float64')
            procData['r2_corr'] *= corrParPickle['r2WallRef']/interpFunc(procData['z_s1_cnn'])
        procData['s1_corr2']  = np.array(procData['area'], dtype='float64')
        procData['s1_corr2'] *= corrParPickle['areaCorr2Pars'][0]
        procData['s1_corr2'] /= linearFunc(procData['r2_corr'], *corrParPickle['areaCorr2Pars'])
        procData['s1_corrF']  = np.array(procData['s1_corr2'], dtype='float64')
        procData['s1_corrF'] *= corrParPickle['areaCorr3Pars'][-1]
        procData['s1_corrF'] /= exp2Func(procData['area_fraction_top'], *corrParPickle['areaCorr3Pars'])

        fidMass_path = self.config['analysis_main']['peak_region_rate_estimation']['fidMass_path'] 
        mask_dict = self.python_module.gen_masks(procData, verbosity=self.verbosity)
        peak_region_rate_estimation(runInfoData, procData, mask_dict, s1corr_paths, fidMass_path,\
                                    output_dir=self.result_dir, plot_dir=self.plot_dir, save_plot=self.save_plot,\
                                    plot_format=self.plot_format,\
                                    output_name_tag=self.config['analysis_main']['output_name_tag'],\
                                    config_dict=self.config['analysis_main']['peak_region_rate_estimation'],\
                                    plot_extra=self.python_module.plot_extra, verbosity=self.verbosity)
    def _surface_rate_estimation(self, args):
        """ Estimate the surface rate through the counting within specified energy regions and FVs. """
        self.config['analysis_main']['surface_rate_estimation']['deadtime_type'] =self.config['general']['deadtime_type']
        self.config['analysis_main']['surface_rate_estimation']['generalBinN']   =self.config['general']['generalBinN']
        rateDict_path = os.path.join(self.result_dir,\
                                     self.config['analysis_main']['output_name_tag']+\
                                     's1only_surface_rate_estimation.pkl')
        runInfoData, procData = self._load_data(args, rateDict_path=rateDict_path)
        runInfoData['run_number'] = np.array([str(run_id).zfill(6) for run_id in runInfoData['number']])
 
        s1corr_paths = self.config['analysis_main']['surface_rate_estimation']['s1corr_paths']
        if self.verbosity >= 1: print('  loading correction_derivation() files:', s1corr_paths)
        corrParPickle = None
        for s1corrPath in s1corr_paths:
            with open(s1corrPath, 'rb') as handle:
                corrParPickle = pickle.load(handle)
        if 'r2CorrPars' in corrParPickle:
            procData['r2_corr']  = np.array(procData['r2_naive'], dtype='float64')
            procData['r2_corr'] *= pow(60.0, 2)/polyFunc(procData['area_fraction_top'], corrParPickle['r2CorrPars'])
        elif 'r2InterpInputs' in corrParPickle: 
            interpFunc = scipy.interpolate.interp1d(corrParPickle['zSpan'], corrParPickle['r2InterpInputs'])
            procData['r2_corr']  = np.array(procData['r2_s1_cnn'], dtype='float64')
            procData['r2_corr'] *= corrParPickle['r2WallRef']/interpFunc(procData['z_s1_cnn'])
        procData['s1_corr2']  = np.array(procData['area'], dtype='float64')
        procData['s1_corr2'] *= corrParPickle['areaCorr2Pars'][0]
        procData['s1_corr2'] /= linearFunc(procData['r2_corr'], *corrParPickle['areaCorr2Pars'])
        procData['s1_corrF']  = np.array(procData['s1_corr2'], dtype='float64')
        procData['s1_corrF'] *= corrParPickle['areaCorr3Pars'][-1]
        procData['s1_corrF'] /= exp2Func(procData['area_fraction_top'], *corrParPickle['areaCorr3Pars'])

        fidMass_path = self.config['analysis_main']['surface_rate_estimation']['fidMass_path']
        mask_dict = self.python_module.gen_masks(procData, verbosity=self.verbosity)
        surface_rate_estimation(runInfoData, procData, mask_dict, s1corr_paths, fidMass_path,output_dir=self.result_dir,\
                                plot_dir=self.plot_dir, save_plot=self.save_plot, plot_format=self.plot_format,\
                                output_name_tag=self.config['analysis_main']['output_name_tag'],\
                                config_dict=self.config['analysis_main']['surface_rate_estimation'],\
                                plot_extra=self.python_module.plot_extra, verbosity=self.verbosity)
    def _custom_process(self, args):
        self.config['analysis_main']['custom_process']['generalBinN']   = self.config['general']['generalBinN']
        self.config['analysis_main']['custom_process']['deadtime_type'] = self.config['general']['deadtime_type']
        runInfoData, procData, mask_dict = None, None, None
        if self.config['analysis_main']['custom_process']['load_data_with_mask'] == True: 
            runInfoData, procData = self._load_data(args)
            mask_dict = self.python_module.gen_masks(procData, verbosity=self.verbosity)
        self.python_module.custom_process(config_dict=self.config,\
                                          runInfoData=runInfoData, procData=procData, mask_dict=mask_dict)
#########################################################################################################################
def correction_derivation(procData, mask_dict, output_dir=None, plot_dir=None, save_plot=None, plot_format=None,\
                          output_name_tag=None, config_dict=None, plot_extra=None, verbosity=None):
    if isinstance(procData, pd.DataFrame) == True:
        procData = procData.to_dict('series')
        for key in procData: procData[key] = np.array(procData[key])
    if output_dir      is None: output_dir  = './'
    if plot_dir        is None: plot_dir    = './'
    if save_plot       is None: save_plot   = True
    if plot_format     is None: plot_format = 'png'
    if output_name_tag is None: output_name_tag = ''
    if verbosity       is None: verbosity   = 1
    dateStr = date.today().strftime('%y%m%d')
    #####################################################################################################################
    zName, r2name, r2Range = None, None, None
    if 'r2_naive'  in procData:
        zName     = 'area_fraction_top' 
        r2name    = 'r2_naive'
        thetaName = 'theta_naive'
        zRange  = [0.0, 0.6]
        r2Range = [0.0, 3.0E2]
    if 'r2_s1_cnn' in procData:
        zName     = 'z_s1_cnn'
        r2name    = 'r2_s1_cnn'
        thetaName = 'theta_s1_cnn'
        zRange  = [-150.0, 0.0]
        r2Range = [0.0, 6.0E3] 
    r2_corrRange = [0.0, 6.0E3]
    if r2name is None: 
        raise AssertionError('correction_derivation(): either r2_naive or r2_s1_cnn are required for s1only analysis')
    
    #####################################################################################################################
    if verbosity >= 1: print('  base data:')
    plot_name = 's1only_correction_derivation_0base' 
    mask0 = np.copy(mask_dict['mask_base'])
    mask1 = np.copy(mask_dict['mask_final_AFT'])
    if save_plot == True:
        plt, ax = plotStandard_s1only_highA(procData, mask0, mask1, verbosity=verbosity,binN=config_dict['generalBinN'],\
                                            r2Key=r2name, r2Range=r2Range, zKey=zName, zRange=zRange)
        if plot_extra is not None: plot_extra(plot_name, ax, mask_dict) 
        plot_path = os.path.join(plot_dir, dateStr+output_name_tag+plot_name+'.'+plot_format)
        plt.savefig(plot_path, format=plot_format)
        if verbosity >= 1: print('    saving:', plot_path)
        plt.close('all')

    #####
    # Correction: AFT vs area
    #####################################################################################################################
    plot_name = 's1only_correction_derivation_1base_fit'
    selAFTupperBound = [30.0E4, -0.2, -0.1, 2.5E4] if config_dict is None else config_dict['selAFTupperBound']
    selAFTlowerBound = [30.0E4, -0.1, -0.1, 5.5E4] if config_dict is None else config_dict['selAFTlowerBound']
    if verbosity >= 1: print('  first correction:')
    ###Corr1 Fit Cuts########################################################################
    ### mask_area_AFT_sel
    mask_gen = np.copy(mask_dict['mask_final_AFT'])
    mask_gen &= (exp2Func(procData['area_fraction_top'], *selAFTupperBound) <= procData['area'])
    mask_gen &= (procData['area'] < exp2Func(procData['area_fraction_top'], *selAFTlowerBound))
    mask_area_AFT_sel = np.copy(mask_gen)
    if verbosity >= 1: print('    mask_area_AFT_sel,', np.sum(mask_area_AFT_sel), ', done')
    ###Corr1 Fit##############################################################################
    maskAFT = np.copy(mask_area_AFT_sel)
    dataX = np.array(procData['area'][maskAFT],              dtype='float64')        
    dataY = np.array(procData['area_fraction_top'][maskAFT], dtype='float64')
    costFunc = iminuit.cost.LeastSquares(dataY, dataX, np.ones(len(dataX)), exp2Func)
    objMinuit = iminuit.Minuit(costFunc, *((np.array(selAFTupperBound) + np.array(selAFTlowerBound))/2.0))
    #for parIdx in range(len(selAFTupperBound)):
    #    objMinuit.limits[parIdx] = [min(selAFTupperBound[parIdx], selAFTlowerBound[parIdx]),\
    #                                max(selAFTupperBound[parIdx], selAFTlowerBound[parIdx])]
    optResult_areaCorr1 = objMinuit.migrad()
    areaCorr1Pars = [val for val in optResult_areaCorr1.values]
    areaCorr1Errs = [val for val in optResult_areaCorr1.errors]
    if verbosity >= 1: print('    area_fraction_top fit:\n', optResult_areaCorr1)

    procData['s1_corr1']  = np.array(procData['area'], dtype='float64')
    procData['s1_corr1'] *= areaCorr1Pars[-1]
    procData['s1_corr1'] /= exp2Func(procData['area_fraction_top'], *areaCorr1Pars)
    mask0 = np.copy(mask_dict['mask_final_AFT'])
    mask1 = np.copy(maskAFT)
    if save_plot == True:
        plt, ax = plotStandard_s1only_highA(procData, mask0, mask1, verbosity=verbosity,binN=config_dict['generalBinN'],\
                                            r2Key=r2name, r2Range=r2Range, zKey=zName, zRange=zRange)
        funcY = np.linspace(*ax[3].get_ylim(), 1000)[:-1]
        ax[3].plot(exp2Func(funcY, *selAFTupperBound), funcY, color='orange', linestyle='dashed', linewidth=1)
        ax[3].plot(exp2Func(funcY, *selAFTlowerBound), funcY, color='orange', linestyle='dashed', linewidth=1)
        ax[3].plot(exp2Func(funcY, *areaCorr1Pars),    funcY, color='orange', linestyle='solid',  linewidth=2)
        if plot_extra is not None: plot_extra(plot_name, ax, mask_dict) 
        plot_path = os.path.join(plot_dir, dateStr+output_name_tag+plot_name+'.'+plot_format)
        plt.savefig(plot_path, format=plot_format)
        if verbosity >= 1: print('    saving:', plot_path)
        plt.close('all')
    ###########################################################
    plot_name = 's1only_correction_derivation_2corr1'
    mask0 = np.copy(mask_dict['mask_base'])
    mask1 = np.copy(mask_dict['mask_final_AFT'])
    if save_plot == True:
        plt, ax = plotStandard_s1only_highA(procData, mask0, mask1, verbosity=verbosity,binN=config_dict['generalBinN'],\
                                            ampKey='s1_corr1', r2Key=r2name, r2Range=r2Range, zKey=zName, zRange=zRange)
        if plot_extra is not None: plot_extra(plot_name, ax, mask_dict) 
        plot_path = os.path.join(plot_dir, dateStr+output_name_tag+plot_name+'.'+plot_format)
        plt.savefig(plot_path, format=plot_format)
        if verbosity >= 1: print('    saving:', plot_path)
        plt.close('all')

    #####
    # Correction: r2_naive
    #####################################################################################################################
    if verbosity >= 1: print('  second correction:') 
    applyR2Correction = True if config_dict is None else config_dict['applyR2Correction']
    selAreaRange = [3.0E4, 4.0E4]                           if config_dict is None else config_dict['r2SelAreaRange']
    selAFTRange  = [0.0, 0.6]                               if config_dict is None else config_dict['r2SelAFTRange']
    r2CorrInit   = [0.35, 30.0,0.0,1.0E3,0.0,1.0E5,0.0,0.0] if config_dict is None else config_dict['r2CorrInit']
    ### mask_area_r2_sel to select Po210
    mask_gen = np.copy(mask_dict['mask_base'])
    mask_gen &= (selAreaRange[0] <= procData['s1_corr1'])          & (procData['s1_corr1']          < selAreaRange[1])
    mask_gen &= (selAFTRange[0]  <= procData['area_fraction_top']) & (procData['area_fraction_top'] < selAFTRange[1])
    if r2name == 'r2_naive':
        mask_gen &= (r2Range[0] <= procData['r2_naive']) & (procData['r2_naive'] < r2Range[1])
    elif r2name == 'r2_s1_cnn':
        mask_gen &= (r2Range[0] <= procData['r2_s1_cnn']) & (procData['r2_s1_cnn'] < r2Range[1])
    mask_area_r2_sel = np.copy(mask_gen)
    if verbosity >= 1: print('    mask_area_r2_sel,', np.sum(mask_area_r2_sel), ', done')
    if applyR2Correction == True:
        ###r2Corr Fit Cuts###########################################################################
        if r2name == 'r2_naive':
            plot_name = 's1only_correction_derivation_3r2corr_fit'
            #Corr2 Fit###############################################################################



### NOTE: WARNING: this parts ignored due to how bad the fits get
            '''
            dataX = np.array(procData['r2_naive'][mask_area_r2_sel],          dtype='float64')        
            dataY = np.array(procData['area_fraction_top'][mask_area_r2_sel], dtype='float64')
            costFunc = iminuit.cost.LeastSquares(dataY, dataX, np.ones(len(dataX)), polyFunc)
            objMinuit = iminuit.Minuit(costFunc, r2CorrInit)
            objMinuit.limits[0] = [mask_dict['S1AFT_cutRange'][0], mask_dict['S1AFT_cutRange'][1]]
            objMinuit.limits[1] = [0.0, np.max(procData['r2_naive'][mask_area_r2_sel])]
            objMinuit.fixed[2] = True
            objMinuit.fixed[4] = True
            objMinuit.fixed[6] = True
            optResult_r2Corr = objMinuit.migrad()
            r2CorrPars = [val for val in optResult_r2Corr.values]
            r2CorrErrs = [val for val in optResult_r2Corr.errors]
            if verbosity >= 1: print('    r2_naive fit:\n', optResult_r2Corr)
            r2WallRef = np.mean(polyFunc(procData['area_fraction_top'], r2CorrPars))
            if 'tpc_r_cm' in mask_dict: r2WallRef = pow(mask_dict['tpc_r_cm'], 2)
            procData['r2_corr']  = np.array(procData['r2_naive'], dtype='float64')
            procData['r2_corr'] *= r2WallRef/polyFunc(procData['area_fraction_top'], r2CorrPars)
            '''
###############################################################################################




            dataX = np.array(procData['area_fraction_top'][mask_area_r2_sel], dtype='float64')
            dataY = np.array(procData['r2_naive'][mask_area_r2_sel],          dtype='float64')
            zSpan, r2InterpInputs = getProfileQuantile(dataX, dataY, quantileRatio=0.5, xRange=zRange, xBinN=40)
            validMask = (0.0 < r2InterpInputs)
            zSpan          = zSpan[validMask]
            r2InterpInputs = r2InterpInputs[validMask]
            zSpan = np.append(-1.0E12, zSpan)
            zSpan = np.append(zSpan, 1.0E12)
            r2InterpInputs = np.append(r2InterpInputs[0], r2InterpInputs)
            r2InterpInputs = np.append(r2InterpInputs, r2InterpInputs[-1])
            interpFunc = scipy.interpolate.interp1d(zSpan, r2InterpInputs)
            r2WallRef = np.mean(r2InterpInputs[1:-1])
            if 'tpc_r_cm' in mask_dict: r2WallRef = pow(mask_dict['tpc_r_cm'], 2)
            procData['r2_corr']  = np.array(procData['r2_naive'], dtype='float64')
            procData['r2_corr'] *= r2WallRef/interpFunc(procData['area_fraction_top']) 

            mask0 = np.copy(mask_dict['mask_final_AFT'])
            mask1 = np.copy(mask_area_r2_sel)
            if save_plot == True:
                plt, ax = plotStandard_s1only_highA(procData, mask0, mask1, verbosity=verbosity,\
                                                    binN=config_dict['generalBinN'], ampKey='s1_corr1',\
                                                    r2Key=r2name, r2Range=r2Range, zKey=zName, zRange=zRange)
                funcY = np.linspace(*ax[5].get_ylim(), 1000)[:-1]
                #ax[5].plot(polyFunc(funcY, r2CorrInit), funcY, color='orange', linestyle='--', linewidth=2)
                #ax[5].plot(polyFunc(funcY, r2CorrPars), funcY, color='orange', linestyle='solid', linewidth=2)
                ax[5].scatter(r2InterpInputs, zSpan, color='orange')
                if plot_extra is not None: plot_extra(plot_name, ax, mask_dict) 
                plot_path = os.path.join(plot_dir, dateStr+output_name_tag+plot_name+'.'+plot_format)
                plt.savefig(plot_path, format=plot_format)
                if verbosity >= 1: print('    saving:', plot_path)
                plt.close('all')
        elif r2name == 'r2_s1_cnn':
            plot_name = 's1only_correction_derivation_3r2corr_interp'
            #Corr2 interpolation ########################################################### 
            dataX = np.array(procData['z_s1_cnn'][mask_area_r2_sel],  dtype='float64')
            dataY = np.array(procData['r2_s1_cnn'][mask_area_r2_sel], dtype='float64')
            zSpan, r2InterpInputs = getProfileQuantile(dataX, dataY, quantileRatio=0.5, xRange=zRange, xBinN=40)
            validMask = (0.0 < r2InterpInputs)
            zSpan          = zSpan[validMask]
            r2InterpInputs = r2InterpInputs[validMask]
            zSpan = np.append(-1.0E12, zSpan)
            zSpan = np.append(zSpan, 1.0E12)
            r2InterpInputs = np.append(r2InterpInputs[0], r2InterpInputs)
            r2InterpInputs = np.append(r2InterpInputs, r2InterpInputs[-1])
            interpFunc = scipy.interpolate.interp1d(zSpan, r2InterpInputs)
            r2WallRef = np.mean(r2InterpInputs[1:-1])   
            if 'tpc_r_cm' in mask_dict: r2WallRef = pow(mask_dict['tpc_r_cm'], 2) 
            procData['r2_corr']  = np.array(procData['r2_s1_cnn'], dtype='float64')
            procData['r2_corr'] *= r2WallRef/interpFunc(procData['z_s1_cnn'])
            
            mask0 = np.copy(mask_dict['mask_base'])
            mask1 = np.copy(mask_area_r2_sel)
            if save_plot == True:
                plt, ax = plotStandard_s1only_highA(procData, mask0, mask1, verbosity=verbosity,\
                                                    binN=config_dict['generalBinN'], ampKey='s1_corr1',\
                                                    r2Key=r2name, r2Range=r2Range, zKey=zName, zRange=zRange)
                ax[5].scatter(r2InterpInputs, zSpan, color='orange')
                if plot_extra is not None: plot_extra(plot_name, ax, mask_dict)
                plot_path = os.path.join(plot_dir, dateStr+output_name_tag+plot_name+'.'+plot_format)
                plt.savefig(plot_path, format=plot_format)
                if verbosity >= 1: print('    saving:', plot_path)
                plt.close('all')
        ###########################################################
        plot_name = 's1only_correction_derivation_4r2corr'
        mask0 = np.copy(mask_dict['mask_base'])
        mask1 = np.copy(mask_dict['mask_final_AFT'])
        mask1 &= (2.5E4 < procData['s1_corr1'])                                                 #NOTE: for display
        if save_plot == True:
            plt, ax = plotStandard_s1only_highA(procData, mask0, mask1, verbosity=verbosity,\
                                                binN=config_dict['generalBinN'], ampKey='s1_corr1',\
                                                r2Key='r2_corr', r2Range=r2_corrRange, zKey=zName, zRange=zRange)
            ######
            axIdx = 1
            hist = ax[axIdx].hist(procData['r2_corr'][mask1], bins=config_dict['generalBinN'], range=r2_corrRange,\
                                  color='steelblue', linewidth=1)
            ax[axIdx].set_xlabel('r2_corr', fontsize=18)
            ax[axIdx].set_ylabel('counts', fontsize=18)

            if plot_extra is not None: plot_extra(plot_name, ax, mask_dict)
            plot_path = os.path.join(plot_dir, dateStr+output_name_tag+plot_name+'.'+plot_format)
            plt.savefig(plot_path, format=plot_format)
            if verbosity >= 1: print('    saving:', plot_path)
            plt.close('all')

    else: 
        procData['r2_corr']  = np.array(procData[r2name], dtype='float64')
        if verbosity >= 1: print('  skipped...')
    #####
    # Correction: r2 vs area
    #####################################################################################################################
    plot_name = 's1only_correction_derivation_5corr1_fit'
    if verbosity >= 1: print('  third correction:')
    selAreaRange    = [4.5E4,  6.0E4] if config_dict is None else config_dict['selAreaRange']
    selR2range      = [0.0,    2.8E3] if config_dict is None else config_dict['selR2range']
    selR2upperBound = [4.0E4, -1.4]   if config_dict is None else config_dict['selR2upperBound']
    selR2lowerBound = [5.8E4, -1.4]   if config_dict is None else config_dict['selR2lowerBound']
    ###Corr2 Fit Cuts########################################################################
    ### mask_area_r2_sel
    mask_gen = np.copy(mask_dict['mask_final_AFT'])
    mask_gen &= (selAreaRange[0] <= procData['s1_corr1']) & (procData['s1_corr1'] < selAreaRange[1])
    mask_gen &= (selR2range[0]   <= procData['r2_corr'])  & (procData['r2_corr']  < selR2range[1])
    mask_gen &= (linearFunc(procData['r2_corr'], *selR2upperBound) <= procData['s1_corr1'])
    mask_gen &= (procData['s1_corr1'] < linearFunc(procData['r2_corr'], *selR2lowerBound))
    mask_area_r2_sel = np.copy(mask_gen)
    if verbosity >= 1: print('    mask_area_r2_sel,', np.sum(mask_area_r2_sel), ', done')
    #Corr2 Fit###############################################################################
    maskR2 = np.copy(mask_area_r2_sel)
    dataX = np.array(procData['s1_corr1'][maskR2], dtype='float64')        
    #dataX = np.array(procData['area'][maskR2],    dtype='float64')  
    dataY = np.array(procData['r2_corr'][maskR2], dtype='float64')
    costFunc = iminuit.cost.LeastSquares(dataY, dataX, np.ones(len(dataX)), linearFunc)
    objMinuit = iminuit.Minuit(costFunc, 4.8E4, -1.4)
    optResult_areaCorr2 = objMinuit.migrad()
    areaCorr2Pars = [val for val in optResult_areaCorr2.values]
    areaCorr2Errs = [val for val in optResult_areaCorr2.errors]
    if verbosity >= 1: print('    z fit:\n', optResult_areaCorr2)

    #NOTE: using 's1_corr1' correction to correct 'area'
    #procData['s1_corr2']  = np.array(procData['s1_corr1'], dtype='float64')
    procData['s1_corr2']  = np.array(procData['area'], dtype='float64')
    procData['s1_corr2'] *= areaCorr2Pars[0]/linearFunc(procData['r2_corr'], *areaCorr2Pars)
    mask0 = np.copy(mask_dict['mask_final'])
    mask1 = np.copy(maskR2)
    if save_plot == True:
        plt, ax = plotStandard_s1only_highA(procData, mask0, mask1, verbosity=verbosity,binN=config_dict['generalBinN'],\
                                            ampKey='s1_corr1', r2Key='r2_corr', r2Range=r2_corrRange, zKey=zName,\
                                            zRange=zRange)
        funcY = np.linspace(*ax[4].get_ylim(), 1000)[:-1]
        ax[4].plot(linearFunc(funcY, *selR2upperBound), funcY, color='orange', linestyle='dashed', linewidth=1)
        ax[4].plot(linearFunc(funcY, *selR2lowerBound), funcY, color='orange', linestyle='dashed', linewidth=1)
        ax[4].plot(linearFunc(funcY, *areaCorr2Pars),   funcY, color='orange', linestyle='solid',  linewidth=1)
        if plot_extra is not None: plot_extra(plot_name, ax, mask_dict) 
        plot_path = os.path.join(plot_dir, dateStr+output_name_tag+plot_name+'.'+plot_format)
        plt.savefig(plot_path, format=plot_format)
        if verbosity >= 1: print('    saving:', plot_path)
        plt.close('all')
    ###########################################################
    plot_name = 's1only_correction_derivation_6corr2'
    mask0 = np.copy(mask_dict['mask_base'])
    mask1 = np.copy(mask_dict['mask_final_AFT'])
    if save_plot == True:
        plt, ax = plotStandard_s1only_highA(procData, mask0, mask1, verbosity=verbosity,binN=config_dict['generalBinN'],\
                                            ampKey='s1_corr2', r2Key='r2_corr', r2Range=r2_corrRange, zKey=zName,\
                                            zRange=zRange)
        if plot_extra is not None: plot_extra(plot_name, ax, mask_dict) 
        plot_path = os.path.join(plot_dir, dateStr+output_name_tag+plot_name+'.'+plot_format)
        plt.savefig(plot_path, format=plot_format)
        if verbosity >= 1: print('    saving:', plot_path) 
        plt.close('all')
  
    #####
    # Correction: AFT vs area redone
    #####################################################################################################################
    plot_name = 's1only_correction_derivation_7corr2_fit'
    selAFTupperBound = [30.0E4, -0.2, -0.1, 4.0E4] if config_dict is None else config_dict['selAFTupperBound_redone']
    selAFTlowerBound = [30.0E4, -0.1, -0.1, 5.5E4] if config_dict is None else config_dict['selAFTlowerBound_redone']
    if verbosity >= 1: print('  forth correction:')
    ###Corr3 Fit Cuts########################################################################
    ### mask_area_AFT_sel
    mask_gen = np.copy(mask_dict['mask_final_AFT'])
    mask_gen &= (exp2Func(procData['area_fraction_top'], *selAFTupperBound) <= procData['s1_corr2'])
    mask_gen &= (procData['s1_corr2'] < exp2Func(procData['area_fraction_top'], *selAFTlowerBound))
    mask_gen &= (selR2range[0]   <= procData['r2_corr'])  & (procData['r2_corr']  < selR2range[1])
    mask_area_AFT_sel = np.copy(mask_gen)
    if verbosity >= 1: print('    mask_area_AFT_sel,', np.sum(mask_area_AFT_sel), ', done')
    ###Corr3 Fit##############################################################################
    maskAFT = np.copy(mask_area_AFT_sel)
    dataX = np.array(procData['s1_corr2'][maskAFT],          dtype='float64')        
    dataY = np.array(procData['area_fraction_top'][maskAFT], dtype='float64')
    costFunc = iminuit.cost.LeastSquares(dataY, dataX, np.ones(len(dataX)), exp2Func)
    objMinuit = iminuit.Minuit(costFunc, *((np.array(selAFTupperBound) + np.array(selAFTlowerBound))/2.0))
    #for parIdx in range(len(selAFTupperBound)):
    #    objMinuit.limits[parIdx] = [min(selAFTupperBound[parIdx], selAFTlowerBound[parIdx]),\
    #                                max(selAFTupperBound[parIdx], selAFTlowerBound[parIdx])]
    optResult_areaCorr3 = objMinuit.migrad()
    areaCorr3Pars = [val for val in optResult_areaCorr3.values]
    areaCorr3Errs = [val for val in optResult_areaCorr3.errors]
    if verbosity >= 1: print('    area_fraction_top fit:\n', optResult_areaCorr3)

    procData['s1_corr3']  = np.array(procData['s1_corr2'], dtype='float64')
    procData['s1_corr3'] *= areaCorr3Pars[-1]
    procData['s1_corr3'] /= exp2Func(procData['area_fraction_top'], *areaCorr3Pars)
    mask0 = np.copy(mask_dict['mask_final_AFT'])
    mask1 = np.copy(maskAFT)
    if save_plot == True:
        plt, ax = plotStandard_s1only_highA(procData, mask0, mask1, verbosity=verbosity,binN=config_dict['generalBinN'],\
                                            ampKey='s1_corr2', r2Key='r2_corr', r2Range=r2_corrRange, zKey=zName,\
                                            zRange=zRange)
        funcY = np.linspace(*ax[3].get_ylim(), 1000)[:-1]
        ax[3].plot(exp2Func(funcY, *selAFTupperBound), funcY, color='orange', linestyle='dashed', linewidth=1)
        ax[3].plot(exp2Func(funcY, *selAFTlowerBound), funcY, color='orange', linestyle='dashed', linewidth=1)
        ax[3].plot(exp2Func(funcY, *areaCorr3Pars),    funcY, color='orange', linestyle='solid',  linewidth=2)
        if plot_extra is not None: plot_extra(plot_name, ax, mask_dict) 
        plot_path = os.path.join(plot_dir, dateStr+output_name_tag+plot_name+'.'+plot_format)
        plt.savefig(plot_path, format=plot_format)
        if verbosity >= 1: print('    saving:', plot_path)
        plt.close('all')
    ###########################################################
    plot_name = 's1only_correction_derivation_8corr3'
    mask0 = np.copy(mask_dict['mask_base'])
    mask1 = np.copy(mask_dict['mask_final_AFT'])
    if save_plot == True:
        plt, ax = plotStandard_s1only_highA(procData, mask0, mask1, verbosity=verbosity,binN=config_dict['generalBinN'],\
                                            ampKey='s1_corr3', r2Key='r2_corr', r2Range=r2_corrRange, zKey=zName,\
                                            zRange=zRange)
        if plot_extra is not None: plot_extra(plot_name, ax, mask_dict) 
        plot_path = os.path.join(plot_dir, dateStr+output_name_tag+plot_name+'.'+plot_format)
        plt.savefig(plot_path, format=plot_format)
        if verbosity >= 1: print('    saving:', plot_path)
        plt.close('all')
    #####################################################################################################################
    procData['s1_corrF'] = procData['s1_corr3'] 

    #####
    # Multi-Peak Fit
    #####################################################################################################################
    plot_name = 's1only_correction_derivation_9peak_fit'
    if verbosity >= 1: print('  multi-peak fit:') 
    s1_ampBinN     = 800            if config_dict is None else config_dict['s1_ampBinN']
    s1_ampRange    = [1.0E4, 1.0E5] if config_dict is None else config_dict['s1_ampRange']
    multipeak_type = 'Rn222'        if config_dict is None else config_dict['multipeak_type']
    fitRangeSet = {}
    if multipeak_type == 'Rn222':
        fitRangeSet['Rn222'] = [4.2E4,  4.7E4]  if config_dict is None else config_dict['fitRangeSet_Rn222']['Rn222']
        fitRangeSet['Po218'] = [4.75E4, 5.15E4] if config_dict is None else config_dict['fitRangeSet_Rn222']['Po218']
        fitRangeSet['Po214'] = [6.1E4,  6.8E4]  if config_dict is None else config_dict['fitRangeSet_Rn222']['Po214']

        fitRangeSet['Rn220'] = [5.2E4, 5.3E4]   if config_dict is None else config_dict['fitRangeSet_Rn222']['Rn220']
        fitRangeSet['Po216'] = [5.5E4, 5.7E4]   if config_dict is None else config_dict['fitRangeSet_Rn222']['Po216']
    elif multipeak_type == 'Rn220':
        fitRangeSet['Rn222'] = [4.2E4, 4.55E4]  if config_dict is None else config_dict['fitRangeSet_Rn220']['Rn222']
        fitRangeSet['Po214'] = [6.1E4, 6.4E4]   if config_dict is None else config_dict['fitRangeSet_Rn220']['Po214']

        fitRangeSet['Bi212'] = [4.7E4, 5.0E4]   if config_dict is None else config_dict['fitRangeSet_Rn220']['Bi212']
        fitRangeSet['Rn220'] = [5.0E4, 5.2E4]   if config_dict is None else config_dict['fitRangeSet_Rn220']['Rn220']
        fitRangeSet['Po216'] = [5.3E4, 5.6E4]   if config_dict is None else config_dict['fitRangeSet_Rn220']['Po216']
        fitRangeSet['Po212'] = [6.8E4, 7.4E4]   if config_dict is None else config_dict['fitRangeSet_Rn220']['Po212']
    else: raise AssertionError('correction_derivation(): the available multipeak_type is [Rn222, Rn220]')
    ###Multi-peak Fit Cuts########################################################################
    mask0 = mask_dict['mask_final_AFT']
    mask1 = mask0 & (selR2range[0] <= procData['r2_corr'])& (procData['r2_corr'] < selR2range[1])
    ###Multi-peak Fit#############################################################################
    fitPars, fitErrs, fitCov = None, None, None
    
    data_s1_area = np.array(procData['s1_corrF'][mask1], dtype='float64')
    s1_areaHist, s1_areaX = np.histogram(data_s1_area, bins=s1_ampBinN, range=s1_ampRange)
    s1_areaXcenter = 0.5*(s1_areaX[1:] + s1_areaX[:-1])
    s1_areaX = s1_areaX[:-1]
    s1_areaErr = np.array([np.sqrt(count) if count != 0 else 1.0 for count in s1_areaHist])
    ###########################################################
    peakFitFunc = multiSigLinkedGausPDF_scipy
    fitDataSet, initPars = {}, []
    fitRangeMask = np.full(len(s1_areaX), False)
    for peakName in fitRangeSet:
        fitDataSet[peakName] = data_s1_area[(fitRangeSet[peakName][0] <= data_s1_area)&
                                            (data_s1_area < fitRangeSet[peakName][1])] 
        fitRangeMask |= ((fitRangeSet[peakName][0] < s1_areaX) & (s1_areaX < fitRangeSet[peakName][1]))
        initPars.extend([len(fitDataSet[peakName])*(s1_areaX[1]-s1_areaX[0]), np.mean(fitDataSet[peakName])])
    initPars.append(0.8*np.std(fitDataSet[list(fitRangeSet.keys())[-1]]))

    costFunc = iminuit.cost.LeastSquares(s1_areaX, s1_areaHist, s1_areaErr, peakFitFunc)
    costFunc.mask = fitRangeMask
    objMinuit = iminuit.Minuit(costFunc, initPars)
    for peakIdx, peakName in enumerate(fitRangeSet):
        objMinuit.limits[2*peakIdx]   = [0, None]
        objMinuit.limits[2*peakIdx+1] = fitRangeSet[peakName]
    objMinuit.limits[-1] = [0, fitRangeSet[peakName][1]-fitRangeSet[peakName][0]]
    optResult = objMinuit.migrad()
    try:
        fitPars = [val for val in optResult.values]
        fitErrs = [err for err in optResult.errors]
        fitCov  = [row for row in optResult.covariance]
    except:
        pass
    if verbosity >= 1:
        print('    multi-peak fit:\n      initPars =', initPars, '\n', optResult)
        for peakName in fitRangeSet: print('    '+peakName+':', fitRangeSet[peakName])
    ###########################################################
    peak_amps, peak_width_percent, peak_normPar, peak_normErr, cali_factors = {}, {}, {}, {}, []
    for peakIdx, peakName in enumerate(fitRangeSet):
        peak_amps[peakName]          = fitPars[2*peakIdx+1]
        peak_width_percent[peakName] = fitPars[-1]*np.sqrt(fitPars[2*peakIdx+1]/fitPars[1])
        peak_width_percent[peakName] = 100.0*peak_width_percent[peakName]/fitPars[2*peakIdx+1]
        peak_normSig  = ufloat(fitPars[-1], fitErrs[-1])
        peak_normSig *= ufloat(fitPars[2*peakIdx+1], fitErrs[2*peakIdx+1])/ufloat(fitPars[1], fitErrs[1])
        peak_normPar[peakName] = [fitPars[2*peakIdx+1], peak_normSig.n]
        peak_normErr[peakName] = [fitErrs[2*peakIdx+1], peak_normSig.s]
        cali_factors.append(ufloat(fitPars[2*peakIdx+1]/ALPHA_LINE_E[peakName],\
                                   fitErrs[2*peakIdx+1]/ALPHA_LINE_E[peakName]))
    cali_factor = weightedAverage(cali_factors).n
    if verbosity >= 1: print('    cali_factor =', cali_factor, '/MeV')
    ###########################################################
    fitX     = np.copy(s1_areaX)
    fitY     = peakFitFunc(fitX, fitPars)
    fitYinit = peakFitFunc(fitX, initPars)
    fitResidual = np.array([np.sqrt(iminuit.cost.poisson_chi2(np.array([1.0*histYval]), np.array([1.0*fitYval])))
                            for histYval, fitYval in zip(s1_areaHist, peakFitFunc(s1_areaX, fitPars))])
    ###########################################################
    if save_plot == True:
        gridSpec = [1, 1]   
        figSize, marginRatio = getSizeMargin(gridSpec, subplotSize=[15.0, 5.0], plotMargin=[0.08, 0.14, 0.04, 0.06])
        fig = plt.figure(figsize=figSize); fig.subplots_adjust(*marginRatio)
        gs = gridspec.GridSpec(*gridSpec)
        matplotlib.rc('xtick', labelsize=18)
        matplotlib.rc('ytick', labelsize=18)
        ax = []
        for axIdx in range(gs.nrows*gs.ncols):
            ax.append(fig.add_subplot(gs[axIdx]));
            ax[-1].ticklabel_format(style='sci', scilimits=(-2, 2))
            ax[-1].xaxis.set_tick_params(width=1.5, length=4)
            ax[-1].xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
            ax[-1].yaxis.set_tick_params(width=1.5, length=4)
            ax[-1].yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
        #####
        axIdx = 0
        hist = ax[axIdx].hist(procData['s1_corrF'][mask0], bins=s1_ampBinN, range=s1_ampRange,
                              color='gainsboro', linewidth=1)
        hist = ax[axIdx].hist(procData['s1_corrF'][mask1], bins=s1_ampBinN, range=s1_ampRange,
                              color='steelblue', linewidth=1)
        ax[axIdx].get_xaxis().set_visible(False)
        ax[axIdx].set_ylabel('counts', fontsize=18)
        ax[axIdx].set_xlim(s1_ampRange)
        ax[axIdx].set_ylim(bottom=0.5)
        ax[axIdx].set_yscale('log')

        divider = make_axes_locatable(ax[axIdx])
        axRes = divider.append_axes('bottom', size='20%', pad='2%')
        ax[axIdx].figure.add_axes(axRes)
        axRes.set_xlim(s1_ampRange)
        axRes.set_xlabel('s1_corrF', fontsize=18)
        axRes.set_ylabel('fit\nresidual', fontsize=18)
        axRes.ticklabel_format(style='sci', scilimits=(-2, 2))
        axRes.xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
        axRes.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
        if fitCov is not None:
            linePlot, linePlotText= {}, {}
            for peakName in fitRangeSet:
                fitMask = ((fitRangeSet[peakName][0] < fitX) & (fitX < fitRangeSet[peakName][1]))
                # ax[axIdx].plot(fitX[fitMask], fitYinit[fitMask], color='orange', linewidth=1.0, alpha=1.0, zorder=2)
                ax[axIdx].plot(fitX[fitMask], fitY[fitMask],     color='red',    linewidth=1.5, alpha=0.8, zorder=3)
                linePlot[peakName] = ax[axIdx].axvline(x=peak_amps[peakName], ymin=0.0, ymax=ax[axIdx].get_ylim()[1],
                                                       color=ALPHA_LINE_C[peakName], linestyle='dashed', linewidth=1.5)
                linePlotText[peakName] = peakName + ', w=' + scientificStr(peak_width_percent[peakName]) + '%'
            for peakName in ALPHA_LINE_E:
                if peakName not in fitRangeSet:
                    peak_amp = cali_factor*ALPHA_LINE_E[peakName]
                    linePlot[peakName] = ax[axIdx].axvline(x=peak_amp, ymin=0.0, ymax=ax[axIdx].get_ylim()[1],
                                                           color=ALPHA_LINE_C[peakName],linestyle='dashed',linewidth=1.5)
                    linePlotText[peakName] = peakName + ''
            legObj = ax[axIdx].legend([linePlot[peakName]     for peakName in ALPHA_LINE_E],\
                                      [linePlotText[peakName] for peakName in ALPHA_LINE_E],\
                                      loc='upper left', fontsize=16)
            # axRes.errorbar(fitX[costFunc.mask], fitResidual[costFunc.mask], yerr=np.sqrt(s1_areaHist[costFunc.mask]),
            #                 color='black', fmt='s', **ERROR_PLOT_PARS)
            for peakName in fitRangeSet:
                fitMask = ((fitRangeSet[peakName][0] < s1_areaX) & (s1_areaX < fitRangeSet[peakName][1]))
                axRes.plot(s1_areaX[fitMask], fitResidual[fitMask], color='black')
        if plot_extra is not None: plot_extra(plot_name, ax, mask_dict) 
        plot_path = os.path.join(plot_dir, dateStr+output_name_tag+plot_name+'.'+plot_format)
        plt.savefig(plot_path, format=plot_format)
        if verbosity >= 1: print('    saving:', plot_path)
        plt.close('all')
    #####################################################################################################################
    corrParDict = {}
    corrParDict['areaCorr1Pars'] = areaCorr1Pars
    corrParDict['areaCorr1Errs'] = areaCorr1Errs
    if applyR2Correction == True:
        if r2name == 'r2_naive':
            corrParDict['r2CorrPars'] = r2CorrPars
            corrParDict['r2CorrErrs'] = r2CorrErrs
        elif r2name == 'r2_s1_cnn':
            corrParDict['zSpan']          = zSpan
            corrParDict['r2InterpInputs'] = r2InterpInputs
            corrParDict['r2WallRef']      = r2WallRef
    corrParDict['areaCorr2Pars'] = areaCorr2Pars
    corrParDict['areaCorr2Errs'] = areaCorr2Errs
    corrParDict['areaCorr3Pars'] = areaCorr3Pars
    corrParDict['areaCorr3Errs'] = areaCorr3Errs
    corrParDict['fitRangeSet'], corrParDict['peak_normPar'], corrParDict['peak_normErr'] = {}, {}, {}
    for peakName in fitRangeSet:
        corrParDict['fitRangeSet'][peakName]  = fitRangeSet[peakName]
        corrParDict['peak_normPar'][peakName] = peak_normPar[peakName]
        corrParDict['peak_normErr'][peakName] = peak_normErr[peakName]
    corrParDict['cali_factor'] = cali_factor
    corrParDict['parameters'] = {}
    for key in mask_dict:
        if 'mask' not in key: corrParDict['parameters'][key] = mask_dict[key]
    output_path = os.path.join(output_dir, output_name_tag+'s1only_correction_derivation.pkl')
    with open(output_path, 'wb') as handle:
        pickle.dump(corrParDict, handle, protocol=pickle.HIGHEST_PROTOCOL)
    if verbosity >= 1: print('    saving:', output_path)
#########################################################################################################################
ONE_SIMGA_LEVEL = 0.6826894921370859 
def estimate_fiducial_mass(runInfoData, procData, mask_dict, procData_S1onlyOnly, mask_dict_S1onlyOnly, s1corr_paths,\
                           s1s2_mask_dict=None, output_dir=None, plot_dir=None, save_plot=None, plot_format=None,\
                           output_name_tag=None, config_dict=None, plot_extra=None, verbosity=None):
    if isinstance(runInfoData, pd.DataFrame) == True:
        runInfoData = runInfoData.to_dict('series')
        for key in runInfoData: runInfoData[key] = np.array(runInfoData[key])
    if isinstance(procData, pd.DataFrame) == True:
        procData = procData.to_dict('series')
        for key in procData: procData[key] = np.array(procData[key])
    if output_dir      is None: output_dir      = './'
    if plot_dir        is None: plot_dir        = './'
    if save_plot       is None: save_plot       = True
    if plot_format     is None: plot_format     = 'png'
    if output_name_tag is None: output_name_tag = ''
    if verbosity       is None: verbosity       = 1
    dateStr = date.today().strftime('%y%m%d')
    #####################################################################################################################
    zName, r2name, r2Range = None, None, None
    if 'r2_naive'  in procData:
        zName   = 'area_fraction_top'
        r2name  = 'r2_naive'
        zRange  = [0.0, 1.0]
        r2Range = [0.0, 3.0E2]
    if 'r2_s1_cnn' in procData:
        zName  = 'z_s1_cnn'
        r2name = 'r2_s1_cnn'
        zRange = [-150.0, 0.0]
        r2Range= [0.0, 6.0E3]
    r2_corrRange = [0.0, 6.0E3]
    if r2name is None:
        raise AssertionError('correction_derivation(): either r2_naive or r2_s1_cnn are required for s1only analysis')
    #####################################################################################################################
    if verbosity >= 1: print('  base data:')
    plot_name = 's1only_estimate_fiducial_mass_0base' 
    mask0 = copy.deepcopy(mask_dict['mask_base'])
    mask1 = copy.deepcopy(mask_dict['mask_final_FV_mass_est'])
    mask0_ = copy.deepcopy(mask_dict_S1onlyOnly['mask_base'])
    mask1_ = copy.deepcopy(mask_dict_S1onlyOnly['mask_final_FV_mass_est'])
    mask2_ = copy.deepcopy(mask_dict_S1onlyOnly['mask_final_FV'])
    S1woS2_altS1_ratio = np.sum(mask1_)/np.sum(mask1)                   #NOTE: S1woS2_altS1_ratio included as systematics
    if save_plot == True:
        plt, ax = plotStandard_s1only_highA(procData, mask0, mask1, with_s1s2=True, verbosity=verbosity,\
                                            binN=config_dict['generalBinN'], r2Key=r2name, r2Range=r2Range,\
                                            zKey=zName, zRange=zRange)
        axIdx = 0
        ax[axIdx].hist(procData_S1onlyOnly['area'][mask0_], bins=config_dict['generalBinN'], 
                       range=ax[0].get_xlim(), color='green', linewidth=2, histtype='step', alpha=0.2)
        ax[axIdx].hist(procData_S1onlyOnly['area'][mask1_], bins=config_dict['generalBinN'],\
                       range=ax[0].get_xlim(), color='green', linewidth=2, histtype='step')
        legObj = ax[axIdx].legend([patches.Rectangle((0,0), 1, 1, color='green')], ['S1woS2/altS1'],\
                                  loc='center right', fontsize=16)
        ax[axIdx].text(0.97, 0.38, str(np.sum(mask0_)), fontsize=18, color='green', transform=ax[0].transAxes, \
                       horizontalalignment='right', verticalalignment='top', zorder=16, alpha=0.2)
        ax[axIdx].text(0.97, 0.44, str(np.sum(mask1_)), fontsize=18, color='green', transform=ax[0].transAxes, \
                       horizontalalignment='right', verticalalignment='top', zorder=16)

        if plot_extra is not None: plot_extra(plot_name, ax, mask_dict) 
        plot_path = os.path.join(plot_dir, dateStr+output_name_tag+plot_name+'.'+plot_format)
        plt.savefig(plot_path, format=plot_format)
        if verbosity >= 1: print('    saving:', plot_path)
        plt.close('all')

    plot_name = 's1only_estimate_fiducial_mass_0base_s1onlyOnly_mass_est_'
    if save_plot == True:
        plt, ax = plotStandard_s1only_highA(procData_S1onlyOnly, mask0_, mask1_, verbosity=verbosity,\
                                            binN=config_dict['generalBinN'], r2Key=r2name, r2Range=r2Range,\
                                            zKey=zName, zRange=zRange)
        if plot_extra is not None: plot_extra(plot_name, ax, mask_dict)
        plot_path = os.path.join(plot_dir, dateStr+output_name_tag+plot_name+'.'+plot_format)
        plt.savefig(plot_path, format=plot_format)
        if verbosity >= 1: print('    saving:', plot_path)
        plt.close('all')

    plot_name = 's1only_estimate_fiducial_mass_0base_s1onlyOnly'
    if save_plot == True:
        plt, ax = plotStandard_s1only_highA(procData_S1onlyOnly, mask0_, mask2_, verbosity=verbosity,\
                                            binN=config_dict['generalBinN'], r2Key=r2name, r2Range=r2Range,\
                                            zKey=zName, zRange=zRange)
        if plot_extra is not None: plot_extra(plot_name, ax, mask_dict)
        plot_path = os.path.join(plot_dir, dateStr+output_name_tag+plot_name+'.'+plot_format)
        plt.savefig(plot_path, format=plot_format)
        if verbosity >= 1: print('    saving:', plot_path)
        plt.close('all')

    ###S1 Corr###############################################################################
    if verbosity >= 1: print('  s1 corrected data:')
    plot_name = 's1only_estimate_fiducial_mass_1s1corr'
    mask0 = np.copy(mask_dict['mask_base'])
    mask1 = np.copy(mask_dict['mask_final_FV_mass_est'])
    mask0_ = np.copy(mask_dict_S1onlyOnly['mask_base'])
    mask1_ = np.copy(mask_dict_S1onlyOnly['mask_final_FV_mass_est'])
    if save_plot == True:
        plt, ax = plotStandard_s1only_highA(procData, mask0, mask1, with_s1s2=True, verbosity=verbosity,\
                                            binN=config_dict['generalBinN'], ampKey='s1_corrF',\
                                            r2Key='r2_corr', r2Range=r2_corrRange, zKey=zName, zRange=zRange)
        axIdx = 0
        ax[axIdx].hist(procData_S1onlyOnly['s1_corrF'][mask0_], bins=config_dict['generalBinN'],\
                       range=ax[0].get_xlim(), color='green', linewidth=2, histtype='step', alpha=0.2)
        ax[axIdx].hist(procData_S1onlyOnly['s1_corrF'][mask1_], bins=config_dict['generalBinN'],\
                       range=ax[0].get_xlim(), color='green', linewidth=2, histtype='step')
        legObj = ax[axIdx].legend([patches.Rectangle((0,0), 1, 1, color='green')], ['S1woS2/altS1'],\
                                  loc='center right', fontsize=16)
        ax[axIdx].text(0.97, 0.38, str(np.sum(mask0_)), fontsize=18, color='green', transform=ax[0].transAxes, \
                       horizontalalignment='right', verticalalignment='top', zorder=16, alpha=0.2)
        ax[axIdx].text(0.97, 0.44, str(np.sum(mask1_)), fontsize=18, color='green', transform=ax[0].transAxes, \
                       horizontalalignment='right', verticalalignment='top', zorder=16)
 
        if plot_extra is not None: plot_extra(plot_name, ax, mask_dict)
        plot_path = os.path.join(plot_dir, dateStr+output_name_tag+plot_name+'.'+plot_format)
        plt.savefig(plot_path, format=plot_format)
        if verbosity >= 1: print('    saving:', plot_path)
        plt.close('all') 
 
    plot_name = 's1only_estimate_fiducial_mass_1s1corr_s1onlyOnly_mass_est'
    if save_plot == True:
        plt, ax = plotStandard_s1only_highA(procData_S1onlyOnly, mask0_, mask1_, verbosity=verbosity,\
                                            binN=config_dict['generalBinN'], ampKey='s1_corrF',\
                                            r2Key='r2_corr', r2Range=r2_corrRange, zKey=zName, zRange=zRange)
        if plot_extra is not None: plot_extra(plot_name, ax, mask_dict)
        plot_path = os.path.join(plot_dir, dateStr+output_name_tag+plot_name+'.'+plot_format)
        plt.savefig(plot_path, format=plot_format)
        if verbosity >= 1: print('    saving:', plot_path)
        plt.close('all')
    plot_name = 's1only_estimate_fiducial_mass_1s1corr_s1onlyOnly'
    if save_plot == True:
        plt, ax = plotStandard_s1only_highA(procData_S1onlyOnly, mask0_, mask2_, verbosity=verbosity,\
                                            binN=config_dict['generalBinN'], ampKey='s1_corrF',\
                                            r2Key='r2_corr', r2Range=r2_corrRange, zKey=zName, zRange=zRange)
        if plot_extra is not None: plot_extra(plot_name, ax, mask_dict)
        plot_path = os.path.join(plot_dir, dateStr+output_name_tag+plot_name+'.'+plot_format)
        plt.savefig(plot_path, format=plot_format)
        if verbosity >= 1: print('    saving:', plot_path)
        plt.close('all')
    #####
    # S1S2 Reference Volume
    #####################################################################################################################
    if verbosity >= 1: print('  s1s2 reference volume:')
    plot_name = 's1only_estimate_fiducial_mass_2refVol'
    mask0 = np.copy(mask_dict['mask_base'])
    mask1 = np.copy(mask_dict['mask_ref_FV']) & s1s2_mask_dict['mask_final']

    ################################## reference volume rz uniformity error
    refVol_r2z_hist2D, r2_edges, z_edges = np.histogram2d(np.power(procData['r'], 2)[mask1], procData['z'][mask1],\
                                                          bins=int(config_dict['generalBinN']/10),\
                                                          range=[mask_dict['r2_cutRange'], mask_dict['z_cutRange']])
    ### removing the edge bins
    refVol_r2z_hist2D_rmEdges = []
    for row in refVol_r2z_hist2D[1:-1]:
        refVol_r2z_hist2D_rmEdges.append(np.array(row[:-1]))
    refVol_r2z_hist2D = np.array(refVol_r2z_hist2D_rmEdges)
    r2_edges = r2_edges[1:-1]
    z_edges  = z_edges[1:-1]
    ###
    refVol_r2z_hist2D_flat = refVol_r2z_hist2D.flatten()
    refVol_r2z_mean = np.mean(refVol_r2z_hist2D_flat)
    refVol_r2z_std  = np.std(refVol_r2z_hist2D_flat)
    refVol_r2z_err = 0
    if refVol_r2z_std > np.sqrt(refVol_r2z_mean): refVol_r2z_err = np.sqrt(pow(refVol_r2z_std, 2) - refVol_r2z_mean)
    refVol_r2z_err_sqrtN = refVol_r2z_err/np.sqrt(len(refVol_r2z_hist2D_flat))                  #sum of gaussians
    refVol_r2z_multErr = ufloat(1.0, refVol_r2z_err_sqrtN/refVol_r2z_mean)
    if verbosity >= 1: 
        print('  reference volume rz uniformity error:')
        print('    refVol_r2z_mean       =', refVol_r2z_mean)
        print('    sqrt(refVol_r2z_mean) =', np.sqrt(refVol_r2z_mean))
        print('    refVol_r2z_std        =', refVol_r2z_std)
        print('    refVol_r2z_err        =', refVol_r2z_err)
        print('    number of bins        =', len(refVol_r2z_hist2D_flat))
        print('    refVol_r2z_err_sqrtN  =', refVol_r2z_err_sqrtN)
        print('    refVol_r2z_multErr    =', refVol_r2z_multErr)
    ################################## reference volume xy uniformity error
    r_range = [-np.sqrt(mask_dict['r2_cutRange'][1]), np.sqrt(mask_dict['r2_cutRange'][1])]
    refVol_xy_hist2D, x_edges, y_edges = np.histogram2d(procData['x'][mask1], procData['y'][mask1],\
                                                        bins=int(config_dict['generalBinN']/10),range=[r_range, r_range])
    ### removing the edge bins
    x_mesh, y_mesh = np.meshgrid(x_edges[:-1]+(x_edges[1]-x_edges[0])/2.0,\
                                 y_edges[:-1]+(y_edges[1]-y_edges[0])/2.0)
    for xIdx, val_x in enumerate(x_mesh[0]):
        for yIdx, y_row in enumerate(y_mesh):
            val_y = copy.copy(y_row[0])
            val_r = np.sqrt(pow(val_x, 2) + pow(val_y, 2))
            if val_r > (r_range[1] - np.sqrt(pow(x_edges[1]-x_edges[0], 2)+pow(y_edges[1]-y_edges[0], 2))):
                refVol_xy_hist2D[xIdx][yIdx] = 0
            if refVol_xy_hist2D[xIdx][yIdx] == 0: refVol_xy_hist2D[xIdx][yIdx] = np.nan
    ###
    refVol_xy_hist2D_flat = refVol_xy_hist2D.flatten()
    refVol_xy_mask = ~np.isnan(refVol_xy_hist2D_flat)
    refVol_xy_mean = np.mean(refVol_xy_hist2D_flat[refVol_xy_mask])
    refVol_xy_std  = np.std(refVol_xy_hist2D_flat[refVol_xy_mask])
    refVol_xy_err = 0
    if refVol_xy_std > np.sqrt(refVol_xy_mean): refVol_xy_err = np.sqrt(pow(refVol_xy_std, 2) - refVol_xy_mean)
    refVol_xy_err_sqrtN = refVol_xy_err/np.sqrt(len(refVol_xy_hist2D_flat[refVol_xy_mask]))          #sum of gaussians
    refVol_xy_multErr   = ufloat(1.0, refVol_xy_err_sqrtN/refVol_xy_mean)
    if verbosity >= 1:
        print('  reference volume xy uniformity error:')
        print('    refVol_xy_mean       =', refVol_xy_mean)
        print('    sqrt(refVol_xy_mean) =', np.sqrt(refVol_xy_mean))
        print('    refVol_xy_std        =', refVol_xy_std)
        print('    refVol_xy_err        =', refVol_xy_err)
        print('    number of bins       =', len(refVol_xy_hist2D_flat[refVol_xy_mask]))
        print('    refVol_xy_err_sqrtN  =', refVol_xy_err_sqrtN)
        print('    refVol_xy_multErr    =', refVol_xy_multErr)

    if save_plot == True:
        plt, ax = plotStandard_s1only_highA(procData, mask0, mask1, with_s1s2=True, verbosity=verbosity,\
                                            binN=config_dict['generalBinN'], ampKey='s1_corrF',\
                                            r2Key='r2_corr', r2Range=r2_corrRange, zKey=zName, zRange=zRange)
        if plot_extra is not None: plot_extra(plot_name, ax, mask_dict)
        plot_path = os.path.join(plot_dir, dateStr+output_name_tag+plot_name+'.'+plot_format)
        plt.savefig(plot_path, format=plot_format)
        if verbosity >= 1: print('    saving:', plot_path)
        plt.close('all')        
        ##################################
        gridSpec = [2, 2]
        figSize, marginRatio = getSizeMargin(gridSpec, subplotSize=[6.0, 5.0], plotMargin=[0.15, 0.14, 0.12, 0.06])
        fig = plt.figure(figsize=figSize); fig.subplots_adjust(*marginRatio)
        gs = gridspec.GridSpec(*gridSpec)
        matplotlib.rc('xtick', labelsize=18)
        matplotlib.rc('ytick', labelsize=18)
        ax = list()
        for axIdx in range(gs.nrows*gs.ncols):
            ax.append(fig.add_subplot(gs[axIdx]));
            ax[-1].ticklabel_format(style='sci', scilimits=(-2, 2))
            ax[-1].xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
            ax[-1].yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
        ###
        axIdx = 2
        rangeXYZ = [[0.0, 5.0E3], [-150.0, 0.0], [None, None]]
        hist = ax[axIdx].imshow(refVol_r2z_hist2D.T, extent=[r2_edges[0], r2_edges[-1], z_edges[0], z_edges[-1]],\
                                cmap=CMAP, aspect='auto')
        hist.set_clim(*rangeXYZ[2])
        divider = make_axes_locatable(ax[axIdx])
        cax = divider.append_axes('right', size='5%', pad=0.1)
        if len(procData['cs1'][mask1]) > 10:
            cbar = fig.colorbar(hist, ax=ax[axIdx], cax=cax)
            cbar.mappable.set_clim(*rangeXYZ[2])
        ax[axIdx].set_xlabel('r2', fontsize=18)
        ax[axIdx].set_xlim(*rangeXYZ[0])
        ax[axIdx].set_ylabel('z', fontsize=18)
        ax[axIdx].set_ylim(*rangeXYZ[1])
        ax[axIdx].text(0.98, 0.98, 'binN='+str(len(refVol_r2z_hist2D_flat)), fontsize=18, color='black',\
                       transform=ax[axIdx].transAxes, horizontalalignment='right', verticalalignment='top',zorder=16)
        ###
        axIdx = 3
        hist = ax[axIdx].hist(refVol_r2z_hist2D_flat, bins=int(config_dict['generalBinN']/10))
        ax[axIdx].set_xlabel('r2z bin count', fontsize=18)
        ax[axIdx].set_ylabel('count', fontsize=18)
        ###
        axIdx = 0
        rangeXYZ = [[-80, 80], [-80, 80], [None, None]]
        hist = ax[axIdx].imshow(refVol_xy_hist2D.T, extent=[x_edges[0], x_edges[-1], y_edges[0], y_edges[-1]], cmap=CMAP)
        hist.set_clim(*rangeXYZ[2])
        divider = make_axes_locatable(ax[axIdx])
        cax = divider.append_axes('right', size='5%', pad=0.1)
        if len(procData['cs1'][mask1]) > 10:
            cbar = fig.colorbar(hist, ax=ax[axIdx], cax=cax)
            cbar.mappable.set_clim(*rangeXYZ[2])
        ax[axIdx].set_xlabel('x', fontsize=18)
        ax[axIdx].set_xlim(*rangeXYZ[0])
        ax[axIdx].set_ylabel('y', fontsize=18)
        ax[axIdx].set_ylim(*rangeXYZ[1])
        ax[axIdx].set_aspect(1, adjustable="box")
        ax[axIdx].text(0.98, 0.98, 'binN='+str(len(refVol_xy_hist2D_flat[refVol_xy_mask])), fontsize=18, color='black',\
                       transform=ax[axIdx].transAxes, horizontalalignment='right', verticalalignment='top', zorder=16)
        ###
        axIdx = 1
        hist = ax[axIdx].hist(refVol_xy_hist2D_flat, bins=int(config_dict['generalBinN']/10))
        ax[axIdx].set_xlabel('xy bin count', fontsize=18)
        ax[axIdx].set_ylabel('count', fontsize=18)
        ###
        plot_path = os.path.join(plot_dir, dateStr+output_name_tag+plot_name+'_.'+plot_format)
        plt.savefig(plot_path, format=plot_format)
        if verbosity >= 1: print('    saving:', plot_path)
        plt.close('all')
    #####
    # S1 Only Fiducial Volume Error Bar using Reference Volume
    #####################################################################################################################
    if verbosity >= 1: print('  s1only fiducial volume:')
    plot_name = 's1only_estimate_fiducial_mass_3fidVol'
    mask0 = np.copy(mask_dict['mask_ref_FV']) & s1s2_mask_dict['mask_final']
    mask1 = np.copy(mask_dict['mask_final_FV_mass_est'] & mask_dict['mask_ref_FV']) & s1s2_mask_dict['mask_final']

    statObj = stats.binomtest(np.sum(mask1), np.sum(mask0))
    insideRefFVRatio     = statObj.proportion_estimate
    insideRefFVRatio_CI  = statObj.proportion_ci(confidence_level=ONE_SIMGA_LEVEL)
    insideRefFVRatio_err = max(insideRefFVRatio-insideRefFVRatio_CI.low, insideRefFVRatio_CI.high-insideRefFVRatio)
    insideRefFVRatio     = ufloat(insideRefFVRatio, insideRefFVRatio_err)
    XeFidVol_noUniErr  = insideRefFVRatio*mask_dict['XeFidVol_ref']
    XeFidMass_noUniErr = insideRefFVRatio*mask_dict['XeFidMass_ref']
    XeFidVol           = insideRefFVRatio*mask_dict['XeFidVol_ref']*refVol_r2z_multErr*refVol_xy_multErr
    XeFidMass          = insideRefFVRatio*mask_dict['XeFidMass_ref']*refVol_r2z_multErr*refVol_xy_multErr
    if verbosity >= 1: 
        print('    insideRefFVRatio     =', insideRefFVRatio)
        print('    insideRefFVRatio_CI  =', insideRefFVRatio_CI)
        print('    insideRefFVRatio_err =', insideRefFVRatio_err)
        print('    XeFidVol_ref  =', mask_dict['XeFidVol_ref'],  'cm^3')
        print('    XeFidMass_ref =', mask_dict['XeFidMass_ref'], 'kg')
        print('    XeFidVol_noUniErr  =', XeFidVol_noUniErr,  'cm^3')
        print('    XeFidMass_noUniErr =', XeFidMass_noUniErr, 'kg')
        print('    XeFidVol           =', XeFidVol,  'cm^3')
        print('    XeFidMass          =', XeFidMass, 'kg')
    if save_plot == True:
        plt, ax = plotStandard_s1only_highA(procData, mask0, mask1, with_s1s2=True, verbosity=verbosity,\
                                            binN=config_dict['generalBinN'], ampKey='s1_corrF',\
                                            r2Key='r2_corr', r2Range=r2_corrRange, zKey=zName, zRange=zRange)
        if plot_extra is not None: plot_extra(plot_name, ax, mask_dict)
        plot_path = os.path.join(plot_dir, dateStr+output_name_tag+plot_name+'.'+plot_format)
        plt.savefig(plot_path, format=plot_format)
        if verbosity >= 1: print('    saving:', plot_path)    
    #####
    # S1 Only Cut Associated Error Bar from Comparing to S1S2 Reference Volume and Cuts
    #####################################################################################################################
    validS1S2Ratio = None                               #NOTE: important correction
    if s1s2_mask_dict is not None:
        if verbosity >= 1: print('  evaluate s1only FV event outside the cuts from s1s2:')

        plot_name = 's1only_estimate_fiducial_mass_4outside'
        mask0 = np.copy(mask_dict['mask_final_FV_mass_est']) 
        mask1 = mask_dict['mask_final_FV_mass_est'] & mask_dict['mask_ref_FV_outside']
        mask2 = mask_dict['mask_final_FV_mass_est'] & (~s1s2_mask_dict['mask_final'])

        statObj = stats.binomtest(np.sum(mask1|mask2), np.sum(mask0))
        outsideS1S2Ratio      = statObj.proportion_estimate
        outsideS1S2Ratio_CI   = statObj.proportion_ci(confidence_level=ONE_SIMGA_LEVEL)
        outsideS1S2Ratio_err  = max(outsideS1S2Ratio-outsideS1S2Ratio_CI.low, outsideS1S2Ratio_CI.high-outsideS1S2Ratio)
        outsideS1S2Ratio_err += np.sum(mask2)/np.sum(mask0)     #NOTE: xenon:xenonnt:lin:sr1:240721updates#section240722
        outsideS1S2Ratio      = ufloat(outsideS1S2Ratio, outsideS1S2Ratio_err)
        validS1S2Ratio        = 1.0 - outsideS1S2Ratio
        validS1S2Ratio       *= ufloat(1.0, S1woS2_altS1_ratio) #NOTE: xenon:xenonnt:lin:sr1:240721updates#section240726
        if verbosity >= 1:
            print('    outsideS1S2Ratio     =', outsideS1S2Ratio)
            print('    outsideS1S2Ratio_CI  =', outsideS1S2Ratio_CI)
            print('    outsideS1S2Ratio_err =', outsideS1S2Ratio_err)
        if save_plot == True:
            plt, ax = plotStandard_s1only_highA(procData, mask0, mask1|mask2, with_s1s2=True, verbosity=verbosity,\
                                                binN=int(config_dict['generalBinN']/4), ampKey='s1_corrF',\
                                                r2Key='r2_corr', r2Range=r2_corrRange, zKey=zName, zRange=zRange)
            if plot_extra is not None: plot_extra(plot_name, ax, mask_dict)
            plot_path = os.path.join(plot_dir, dateStr+output_name_tag+plot_name+'S1S2FV_quality.'+plot_format)
            plt.savefig(plot_path, format=plot_format)
            if verbosity >= 1: print('    saving:', plot_path)
            plt.close('all')
            
            plt, ax = plotStandard_s1only_highA(procData, mask0, mask1, with_s1s2=True, verbosity=verbosity,\
                                                binN=int(config_dict['generalBinN']/4), ampKey='s1_corrF',\
                                                r2Key='r2_corr', r2Range=r2_corrRange, zKey=zName, zRange=zRange)
            if plot_extra is not None: plot_extra(plot_name, ax, mask_dict)
            plot_path = os.path.join(plot_dir, dateStr+output_name_tag+plot_name+'S1S2FV.'+plot_format)
            plt.savefig(plot_path, format=plot_format)
            if verbosity >= 1: print('    saving:', plot_path)
            plt.close('all')

            plt, ax = plotStandard_s1only_highA(procData, mask0, mask2 & mask_dict['mask_ref_FV'], with_s1s2=True,\
                                                verbosity=verbosity, binN=int(config_dict['generalBinN']/4), \
                                                ampKey='s1_corrF', r2Key='r2_corr', r2Range=r2_corrRange, zKey=zName,\
                                                zRange=zRange)
            if plot_extra is not None: plot_extra(plot_name, ax, mask_dict)
            plot_path = os.path.join(plot_dir, dateStr+output_name_tag+plot_name+'S1S2quality.'+plot_format)
            plt.savefig(plot_path, format=plot_format)
            if verbosity >= 1: print('    saving:', plot_path)
            plt.close('all')
    if verbosity >= 1: print('    validS1S2Ratio =', validS1S2Ratio)
    #####
    # S1 Only Z Span Error Bar using Reference Volume
    #####################################################################################################################
    if verbosity >= 1: print('  s1only z span:')
    plot_name = 's1only_estimate_fiducial_mass_5zspan'
    mask0 = np.copy(mask_dict['mask_ref_FV'] & s1s2_mask_dict['mask_final'])
    mask1 = np.copy(mask_dict['mask_final_FV_z_est'] & mask_dict['mask_ref_FV'] & s1s2_mask_dict['mask_final'])

    statObj = stats.binomtest(np.sum(mask1), np.sum(mask0))
    insideRefZRatio     = statObj.proportion_estimate
    insideRefZRatio_CI  = statObj.proportion_ci(confidence_level=ONE_SIMGA_LEVEL)
    insideRefZRatio_err = max(insideRefZRatio-insideRefZRatio_CI.low, insideRefZRatio_CI.high-insideRefZRatio)
    insideRefZRatio     = ufloat(insideRefZRatio, insideRefZRatio_err)
    XeFidZ = insideRefZRatio*mask_dict['XeFidZ_ref']
    
    #fit for fiducial surface on the wall
    dataX = np.array(np.power(procData['r'], 2)[mask1], dtype='float64')
    dataY = np.array(procData['z'][mask1], dtype='float64')
    r2Span, z_99pQuantile = getProfileQuantile(dataX, dataY, quantileRatio=0.99,xRange=mask_dict['r2_cutRange'],xBinN=20)
    _,      z_01pQuantile = getProfileQuantile(dataX, dataY, quantileRatio=0.01,xRange=mask_dict['r2_cutRange'],xBinN=20)
    s1s2_z_err = mask_dict['s1s2_z_err'] + 0.0*r2Span

    costFunc  = iminuit.cost.LeastSquares(r2Span, z_99pQuantile, s1s2_z_err, linearFunc)
    objMinuit = iminuit.Minuit(costFunc, *[np.mean(z_99pQuantile), 0.0])
    optResult = objMinuit.migrad()
    fitPars_99pQ = [ufloat(val, err) for (val, err) in zip(optResult.values, optResult.errors)]
    costFunc  = iminuit.cost.LeastSquares(r2Span, z_01pQuantile, s1s2_z_err, linearFunc)
    objMinuit = iminuit.Minuit(costFunc, *[np.mean(z_01pQuantile), 0.0])
    optResult = objMinuit.migrad() 
    fitPars_01pQ = [ufloat(val, err) for (val, err) in zip(optResult.values, optResult.errors)]
    r2_wall = pow(mask_dict['tpc_r_cm'], 2)
    XeFidZ_wall = linearFunc(r2_wall, *fitPars_99pQ) - linearFunc(r2_wall, *fitPars_01pQ)
    if verbosity >= 1: 
        print('    insideRefZRatio     =', insideRefZRatio)
        print('    insideRefZRatio_CI  =', insideRefZRatio_CI)
        print('    insideRefZRatio_err =', insideRefZRatio_err)
        print('    XeFidZ =', XeFidZ,'cm')
        print('    fitPars_99pQ =', fitPars_99pQ)
        print('    fitPars_01pQ =', fitPars_01pQ)
        print('    r2_wall      =', r2_wall)
        print('    XeFidZ_wall  =', XeFidZ_wall, 'cm')        
    if save_plot == True:
        plt, ax = plotStandard_s1only_highA(procData, mask0, mask1, with_s1s2=True, verbosity=verbosity,\
                                            binN=config_dict['generalBinN'], ampKey='s1_corrF',\
                                            r2Key='r2_corr', r2Range=r2_corrRange, zKey=zName, zRange=zRange)
        ax[6].axvline(x=r2_wall, ymin=0.0, ymax=1.0, color='orange', linestyle='dashed', linewidth=1.5)
        ax[7].scatter(r2Span, z_99pQuantile, color='orange')
        ax[7].scatter(r2Span, z_01pQuantile, color='orange')
        ax[7].axvline(x=r2_wall, ymin=0.0, ymax=1.0, color='orange', linestyle='dashed', linewidth=1.5)
        funcX = np.linspace(*ax[7].get_xlim(), 1000)[:-1]
        ax[7].plot(funcX, [item.n for item in linearFunc(funcX, *fitPars_99pQ)], color='orange', linewidth=2)
        ax[7].plot(funcX, [item.n for item in linearFunc(funcX, *fitPars_01pQ)], color='orange', linewidth=2)
        if plot_extra is not None: plot_extra(plot_name, ax, mask_dict)
        plot_path = os.path.join(plot_dir, dateStr+output_name_tag+plot_name+'.'+plot_format)
        plt.savefig(plot_path, format=plot_format)
        if verbosity >= 1: print('    saving:', plot_path)    
    #####
    # S1 Only Z Span Cut Associated Error Bar from Comparing to S1S2 Reference Volume and Cuts
    #####################################################################################################################
    validS1S2ZRatio = None
    if s1s2_mask_dict is not None:
        if verbosity >= 1: print('  evaluate s1only FV event outside the cuts from s1s2:')
        plot_name = 's1only_estimate_fiducial_mass_6outsideS1S2_zspan'
        mask0 = np.copy(mask_dict['mask_final_FV_z_est'] & mask_dict['mask_ref_FV'])
        mask1 = np.copy(mask_dict['mask_final_FV_z_est'] & mask_dict['mask_ref_FV'] & (~s1s2_mask_dict['mask_final']))

        statObj = stats.binomtest(np.sum(mask1), np.sum(mask0))
        outsideS1S2ZRatio     = statObj.proportion_estimate
        outsideS1S2ZRatio_CI  = statObj.proportion_ci(confidence_level=ONE_SIMGA_LEVEL)
        outsideS1S2ZRatio_err = max(outsideS1S2ZRatio-outsideS1S2ZRatio_CI.low,\
                                    outsideS1S2ZRatio_CI.high-outsideS1S2ZRatio)
        outsideS1S2ZRatio     = ufloat(outsideS1S2ZRatio, outsideS1S2ZRatio_err)
        validS1S2ZRatio       = 1.0 - outsideS1S2ZRatio
        validS1S2ZRatio      *= ufloat(1.0, S1woS2_altS1_ratio)
        if verbosity >= 1:
            print('    outsideS1S2ZRatio     =', outsideS1S2ZRatio)
            print('    outsideS1S2ZRatio_CI  =', outsideS1S2ZRatio_CI)
            print('    outsideS1S2ZRatio_err =', outsideS1S2ZRatio_err)
        if save_plot == True:
            plt, ax = plotStandard_s1only_highA(procData, mask0, mask1, with_s1s2=True, verbosity=verbosity,\
                                                binN=int(config_dict['generalBinN']/4), ampKey='s1_corrF',\
                                                r2Key='r2_corr', r2Range=r2_corrRange, zKey=zName, zRange=zRange)
            if plot_extra is not None: plot_extra(plot_name, ax, mask_dict)
            plot_path = os.path.join(plot_dir, dateStr+output_name_tag+plot_name+'.'+plot_format)
            plt.savefig(plot_path, format=plot_format)
            if verbosity >= 1: print('    saving:', plot_path)
            plt.close('all')
    if verbosity >= 1: print('    validS1S2ZRatio =', validS1S2ZRatio)
    #####
    # S1 Only R2 Span Estimated Error Bar using Reference Volume
    #####################################################################################################################
    if verbosity >= 1: print('  s1only r2 span:')
    plot_name = 's1only_estimate_fiducial_mass_7r2Span'
    mask0 = np.copy(mask_dict['mask_ref_FV'] & s1s2_mask_dict['mask_final'])
    mask1 = np.copy(mask_dict['mask_final_FV_r2_est'] & mask_dict['mask_ref_FV'] & s1s2_mask_dict['mask_final'])

    statObj = stats.binomtest(np.sum(mask1), np.sum(mask0))
    insideRefR2Ratio     = statObj.proportion_estimate
    insideRefR2Ratio_CI  = statObj.proportion_ci(confidence_level=ONE_SIMGA_LEVEL)
    insideRefR2Ratio_err = max(insideRefR2Ratio-insideRefR2Ratio_CI.low, insideRefR2Ratio_CI.high-insideRefR2Ratio)
    insideRefR2Ratio     = ufloat(insideRefR2Ratio, insideRefR2Ratio_err)
    XeFidR2 = insideRefR2Ratio*mask_dict['XeFidR2_ref']
    if verbosity >= 1: 
        print('    insideRefR2Ratio     =', insideRefR2Ratio)
        print('    insideRefR2Ratio_CI  =', insideRefR2Ratio_CI)
        print('    insideRefR2Ratio_err =', insideRefR2Ratio_err)
        print('    XeFidR2 =', XeFidR2,  'cm^2')
    if save_plot == True:
         plt, ax = plotStandard_s1only_highA(procData, mask0, mask1, with_s1s2=True, verbosity=verbosity,\
                                             binN=config_dict['generalBinN'], ampKey='s1_corrF',\
                                             r2Key='r2_corr', r2Range=r2_corrRange, zKey=zName, zRange=zRange)
         if plot_extra is not None: plot_extra(plot_name, ax, mask_dict)
         plot_path = os.path.join(plot_dir, dateStr+output_name_tag+plot_name+'.'+plot_format)
         plt.savefig(plot_path, format=plot_format)
         if verbosity >= 1: print('    saving:', plot_path)    
    #####
    # S1 Only R2 Span Cut Associated Error Bar from Comparing to S1S2 Reference Volume and Cuts
    #####################################################################################################################
    validS1S2R2Ratio = None
    if s1s2_mask_dict is not None:
        if verbosity >= 1: print('  evaluate s1only FV event outside the cuts from s1s2:')
        plot_name = 's1only_estimate_fiducial_mass_8outsideS1S2_r2Span'
        mask0 = np.copy(mask_dict['mask_final_FV_r2_est'] & mask_dict['mask_ref_FV'])
        mask1 = np.copy(mask_dict['mask_final_FV_r2_est'] & mask_dict['mask_ref_FV'] & (~s1s2_mask_dict['mask_final']))

        statObj = stats.binomtest(np.sum(mask1), np.sum(mask0))
        outsideS1S2R2Ratio     = statObj.proportion_estimate
        outsideS1S2R2Ratio_CI  = statObj.proportion_ci(confidence_level=ONE_SIMGA_LEVEL)
        outsideS1S2R2Ratio_err = max(outsideS1S2R2Ratio-outsideS1S2R2Ratio_CI.low,\
                                     outsideS1S2R2Ratio_CI.high-outsideS1S2R2Ratio)
        outsideS1S2R2Ratio     = ufloat(outsideS1S2R2Ratio, outsideS1S2R2Ratio_err)
        validS1S2R2Ratio       = 1.0 - outsideS1S2R2Ratio
        validS1S2R2Ratio      *= ufloat(1.0, S1woS2_altS1_ratio)
        if verbosity >= 1:
            print('    outsideS1S2R2Ratio     =', outsideS1S2R2Ratio)
            print('    outsideS1S2R2Ratio_CI  =', outsideS1S2R2Ratio_CI)
            print('    outsideS1S2R2Ratio_err =', outsideS1S2R2Ratio_err)
        if save_plot == True:
            plt, ax = plotStandard_s1only_highA(procData, mask0, mask1, with_s1s2=True, verbosity=verbosity,\
                                                binN=int(config_dict['generalBinN']/4), ampKey='s1_corrF',\
                                                r2Key='r2_corr', r2Range=r2_corrRange, zKey=zName, zRange=zRange)
            if plot_extra is not None: plot_extra(plot_name, ax, mask_dict)
            plot_path = os.path.join(plot_dir, dateStr+output_name_tag+plot_name+'.'+plot_format)
            plt.savefig(plot_path, format=plot_format)
            if verbosity >= 1: print('    saving:', plot_path)
            plt.close('all')
    if verbosity >= 1: print('    validS1S2R2Ratio =', validS1S2R2Ratio)

    #####
    # S1 Only BiPi Acceptance from Comparing to S1S2 Alt_S1
    #####################################################################################################################
    '''
    time_diffBinN          = 200            
    time_diffRange_BiPo212 = [0.0E3, 5.0E3] 
    time_diffRange_BiPo214 = [0.0E3, 2.0E6]
    if config_dict is not None:
        if 'time_diffBinN'          in config_dict: time_diffBinN          = config_dict['time_diffBinN']
        if 'time_diffRange_BiPo212' in config_dict: time_diffRange_BiPo212 = config_dict['time_diffRange_BiPo212']
        if 'time_diffRange_BiPo214' in config_dict: time_diffRange_BiPo214 = config_dict['time_diffRange_BiPo214']
    if s1s2_mask_dict is not None:
        if verbosity >= 1: print('  evaluate s1only BiPo Po214 event efficiency from s1s2 alt_s1:')
        plot_name = 's1only_estimate_fiducial_mass_9BiPoPo214_altS1'
        mask0 = np.copy(mask_dict['mask_final_FV'])
        mask1 = np.copy(mask_dict['mask_final_FV']) & mask_dict['mask_s1only_Po214']
        if save_plot == True:
            plt, ax = plotStandard_s1only_highA(procData, mask0, mask1, with_s1s2=True, verbosity=verbosity,\
                                                binN=int(config_dict['generalBinN']/4), ampKey='s1_area',\
                                                r2Key='r2_corr', r2Range=r2_corrRange, zKey=zName, zRange=zRange)
            if 'alt_s1_area' in procData:
                axIdx = 0
                altAmpRange = [0.0, 5.0E4]
                hist0 = ax[axIdx].hist(procData['alt_s1_area'][mask0], bins=int(config_dict['generalBinN']/4),\
                                       range=altAmpRange, color='green', linewidth=2, histtype='step', alpha=0.2)
                hist1 = ax[axIdx].hist(procData['alt_s1_area'][mask1], bins=int(config_dict['generalBinN']/4),\
                                       range=altAmpRange, color='green', linewidth=2, histtype='step')
                ax[axIdx].set_xlabel('s1_area/alt_s1_area', fontsize=18) 
            if 's1_time_diff' in procData:
                axIdx = 1
                hist0 = ax[axIdx].hist(procData['s1_time_diff'][mask0], bins=time_diffBinN,\
                                       range=time_diffRange_BiPo214, color='lightgray', linewidth=2, histtype='step')
                hist1 = ax[axIdx].hist(procData['s1_time_diff'][mask1], bins=time_diffBinN,\
                                       range=time_diffRange_BiPo214, color='blue', linewidth=2, histtype='step')
                ax[axIdx].set_xlabel('s1_time - alt_s1_time', fontsize=18)
                ax[axIdx].set_ylabel('counts', fontsize=18)
                ax[axIdx].set_xlim(time_diffRange_BiPo214)
                ax[axIdx].set_ylim(bottom=0.5)
                ax[axIdx].set_yscale('log')

            if plot_extra is not None: plot_extra(plot_name, ax, mask_dict)
            plot_path = os.path.join(plot_dir, dateStr+output_name_tag+plot_name+'.'+plot_format)
            plt.savefig(plot_path, format=plot_format)
            if verbosity >= 1: print('    saving:', plot_path)
            plt.close('all')
    '''
    #####################################################################################################################
    FVDict = {}
    for key in mask_dict:
        if '_cutRange' in key: FVDict[key] = mask_dict[key]
    FVDict['count_final_FV_mass_est'] = np.sum(mask_dict['mask_final_FV_mass_est'])
    FVDict['count_ref_FV']            = np.sum(mask_dict['mask_ref_FV'])
    FVDict['count_in_ref_FV']         = np.sum(mask_dict['mask_final_FV_mass_est'] &   mask_dict['mask_ref_FV'])
    FVDict['count_out_ref_FV']        = np.sum(mask_dict['mask_final_FV_mass_est'] & (~mask_dict['mask_ref_FV']))
    FVDict['refVol_r2z_hist2D']  = refVol_r2z_hist2D
    FVDict['refVol_r2z_mean']    = refVol_r2z_mean
    FVDict['refVol_r2z_std']     = refVol_r2z_std
    FVDict['refVol_r2z_err']     = refVol_r2z_err
    FVDict['refVol_r2z_multErr'] = refVol_r2z_multErr
    FVDict['refVol_xy_hist2D']   = refVol_xy_hist2D
    FVDict['refVol_xy_mean']     = refVol_xy_mean
    FVDict['refVol_xy_std']      = refVol_xy_std
    FVDict['refVol_xy_err']      = refVol_xy_err
    FVDict['refVol_xy_multErr']  = refVol_xy_multErr
    FVDict['XeDensity']        = mask_dict['XeDensity']
    FVDict['XeFidVol_ref']     = mask_dict['XeFidVol_ref']
    FVDict['XeFidMass_ref']    = mask_dict['XeFidMass_ref']
    FVDict['insideRefFVRatio'] = insideRefFVRatio
    FVDict['XeFidVol_noUniErr']  = XeFidVol_noUniErr
    FVDict['XeFidMass_noUniErr'] = XeFidMass_noUniErr
    FVDict['XeFidVol']           = XeFidVol
    FVDict['XeFidMass']          = XeFidMass
    FVDict['validS1S2Ratio']   = validS1S2Ratio

    FVDict['XeFidZ_ref']      = mask_dict['XeFidZ_ref']
    FVDict['insideRefZRatio'] = insideRefZRatio
    FVDict['XeFidZ']          = XeFidZ
    FVDict['XeFidZ_wall']     = XeFidZ_wall
    FVDict['validS1S2ZRatio'] = validS1S2ZRatio

    FVDict['XeFidR2_ref']      = mask_dict['XeFidR2_ref']
    FVDict['insideRefR2Ratio'] = insideRefR2Ratio
    FVDict['XeFidR2']          = XeFidR2
    FVDict['validS1S2R2Ratio'] = validS1S2R2Ratio
    output_path = os.path.join(output_dir, output_name_tag+'s1only_estimate_fiducial_mass.pkl')
    with open(output_path, 'wb') as handle:
        pickle.dump(FVDict, handle, protocol=pickle.HIGHEST_PROTOCOL)
    if verbosity >= 1: print('    saving:', output_path) 
#########################################################################################################################
def genCountRange(peakName, mu, sig, sigN):
    if peakName == 'Rn220': return [mu, mu+sigN*sig]
    if peakName == 'Po212': return [mu-sigN*sig, mu+5*sigN*sig]
    return [mu-sigN*sig, mu+sigN*sig]
def get_activity(peakName, countDataSet, countRangeSet, peakMuSigSet, dataTime, XeFidMass,
                 timeDiff_Po212_efficiency=None, Bi212betaBranchingRatio=None):
    activity = ufloat(countDataSet[peakName], max(1.0, np.sqrt(countDataSet[peakName])))
    if peakMuSigSet is not None:
        if peakMuSigSet[peakName] is None:
            if timeDiff_Po212_efficiency is None: timeDiff_Po212_efficiency = 1.0
            if Bi212betaBranchingRatio   is None: Bi212betaBranchingRatio   = 1.0
            activity /= timeDiff_Po212_efficiency
            activity /= Bi212betaBranchingRatio
        elif peakName in ['Rn220', 'Po218', 'Bi212']:
            closePeakName = None
            if peakName in ['Rn220']:
                if   'Bi212' in countDataSet: closePeakName = 'Bi212'
                elif 'Po218' in countDataSet: closePeakName = 'Po218'
            elif peakName in ['Po218', 'Bi212']:
                if 'Rn220' in countDataSet: closePeakName = 'Rn220'
            closePeakActivity = 0
            if closePeakName is not None:
                closePeakActivity  = ufloat(countDataSet[closePeakName], np.sqrt(countDataSet[closePeakName]))
                closePeakActivity /= (gausCDF_err(countRangeSet[closePeakName][1], 1.0, *peakMuSigSet[closePeakName]) 
                                     -gausCDF_err(countRangeSet[closePeakName][0], 1.0, *peakMuSigSet[closePeakName]))
                closePeakActivity *= (gausCDF_err(countRangeSet[peakName][1],      1.0, *peakMuSigSet[closePeakName]) 
                                     -gausCDF_err(countRangeSet[peakName][0],      1.0, *peakMuSigSet[closePeakName]))
                activity -= closePeakActivity
            activity /= (gausCDF_err(countRangeSet[peakName][1], 1.0, *peakMuSigSet[peakName]) 
                        -gausCDF_err(countRangeSet[peakName][0], 1.0, *peakMuSigSet[peakName]))
        else:
            activity /= (gausCDF_err(countRangeSet[peakName][1], 1.0, *peakMuSigSet[peakName])
                        -gausCDF_err(countRangeSet[peakName][0], 1.0, *peakMuSigSet[peakName]))
    activity *= 1.0E6/(XeFidMass*dataTime)                                                      #uBq/kg
    return activity
def peak_region_rate_estimation(runInfoData, procData, mask_dict, s1corr_paths, fidMass_path, output_dir=None,\
                                plot_dir=None, save_plot=None, plot_format=None, output_name_tag=None, config_dict=None,\
                                plot_extra=None, verbosity=None):
    if isinstance(runInfoData, pd.DataFrame) == True:
        runInfoData = runInfoData.to_dict('series')
        for key in runInfoData: runInfoData[key] = np.array(runInfoData[key])
    if isinstance(procData, pd.DataFrame) == True:
        procData = procData.to_dict('series')
        for key in procData: procData[key] = np.array(procData[key])
    if output_dir      is None: output_dir      = './'
    if plot_dir        is None: plot_dir        = './'
    if save_plot       is None: save_plot       = True
    if plot_format     is None: plot_format     = 'png'
    if output_name_tag is None: output_name_tag = ''
    if verbosity       is None: verbosity       = 1
    dateStr = date.today().strftime('%y%m%d')

    lightYieldRefPath = s1corr_paths[0]
    deadtime_type    = None
    generalBinN      = 400
    countRangeSigN   = 3
    chunkCheckTags   = False
    combTimeThres    = 6*3600.0
    s1_ampBinN       = 800
    s1_ampRange      = [1.0E4, 1.0E5]
    chi2DisplayRange = [0.0, 'inf']
    dateDisplayRange = ['-inf', 'inf']
    peakDisplaySubset = ['Rn222']
    if config_dict is not None:
        if 'light_yield_ref_path' in config_dict: lightYieldRefPath = config_dict['light_yield_ref_path']
        if 'deadtime_type'     in config_dict: deadtime_type     = config_dict['deadtime_type']
        if 'generalBinN'       in config_dict: generalBinN       = config_dict['generalBinN']
        if 'countRangeSigN'    in config_dict: countRangeSigN    = config_dict['countRangeSigN']
        if 'chunkCheckTags'    in config_dict: chunkCheckTags = config_dict['chunkCheckTags']
        if 'combTimeThres'     in config_dict: combTimeThres     = config_dict['combTimeThres']
        if 's1_ampBinN'        in config_dict: s1_ampBinN        = config_dict['s1_ampBinN']
        if 's1_ampRange'       in config_dict: s1_ampRange       = config_dict['s1_ampRange']
        if 'chi2DisplayRange'  in config_dict: chi2DisplayRange  = config_dict['chi2DisplayRange']
        if 'dateDisplayRange'  in config_dict: dateDisplayRange  = config_dict['dateDisplayRange']
        if 'peakDisplaySubset' in config_dict: peakDisplaySubset = config_dict['peakDisplaySubset']
    inputPars = {'s1corr_paths':         s1corr_paths,
                 'light_yield_ref_path': lightYieldRefPath,
                 'deadtime_type':     deadtime_type,
                 'generalBinN':       generalBinN,
                 'countRangeSigN':    countRangeSigN,
                 'chunkCheckTags':    chunkCheckTags,
                 'combTimeThres':     combTimeThres,
                 's1_ampBinN':        s1_ampBinN,
                 's1_ampRange':       s1_ampRange,
                 'chi2DisplayRange':  chi2DisplayRange,
                 'dateDisplayRange':  dateDisplayRange,
                 'peakDisplaySubset': peakDisplaySubset}
    for chi2DispIdx in range(len(chi2DisplayRange)):
        if chi2DisplayRange[chi2DispIdx] =='-inf': chi2DisplayRange[chi2DispIdx] =-np.inf
        if chi2DisplayRange[chi2DispIdx] == 'inf': chi2DisplayRange[chi2DispIdx] = np.inf
 
    if verbosity >= 1: print('  loading light yield reference file:', lightYieldRefPath)
    lightYieldRefPickle, cali_factor = None, None
    with open(lightYieldRefPath, 'rb') as handle:
        lightYieldRefPickle = pickle.load(handle)
    cali_factor = copy.deepcopy(lightYieldRefPickle['cali_factor'])
    if verbosity >= 1: print('    cali_factor =', cali_factor, '/MeV')

    if verbosity >= 1: print('  loading correction_derivation() files:', s1corr_paths) 
    corrPar, corrParPickle = {}, None
    fitRangeSetPickle, peak_normParPickle, peak_normErrPickle = {}, {}, {}
    for s1corrPath in s1corr_paths:
        with open(s1corrPath, 'rb') as handle:
            corrParPickle = pickle.load(handle)
        norm_factor = cali_factor/corrParPickle['cali_factor']
        for key in corrParPickle:
            if key not in corrPar:
                if isinstance(corrParPickle[key], dict) == False: 
                    corrPar[key] = np.array(corrParPickle[key])
                    if verbosity >= 1: print('   ', key, '=', corrPar[key])
        for peakName in corrParPickle['fitRangeSet']:
            if peakName not in fitRangeSetPickle:
                if (peakName not in ['Po218', 'Bi212']) or\
                   ((peakName == 'Po218') and ('Bi212' not in fitRangeSetPickle)) or\
                   ((peakName == 'Bi212') and ('Po218' not in fitRangeSetPickle)):
                    fitRangeSetPickle [peakName] = np.array(corrParPickle['fitRangeSet'][peakName]) *norm_factor
                    peak_normParPickle[peakName] = np.array(corrParPickle['peak_normPar'][peakName])*norm_factor
                    peak_normErrPickle[peakName] = np.array(corrParPickle['peak_normErr'][peakName])*norm_factor
                    if verbosity >= 1: 
                        print('    fitRangeSet['+peakName+']  =', fitRangeSetPickle[peakName])
                        print('    peak_normPar['+peakName+'] =', peak_normParPickle[peakName])
                        print('    peak_normErr['+peakName+'] =', peak_normErrPickle[peakName])

    if verbosity >= 1: print('  loading estimate_fiducial_mass() file:', fidMass_path)
    fidMassPar = None
    with open(fidMass_path, 'rb') as handle:
        fidMassPar = pickle.load(handle)
    if verbosity >= 1:
        for key in fidMassPar: print('   ', key, '=', fidMassPar[key])

    fitLightYieldCorrRanges   = None        # for scaling to find the best fit ranges, default=fitLightYieldCorrPeaks
    fitLightYieldCorrPeaks    = ['Rn222', 'Po218', 'Rn220', 'Po216', 'Po214']
    lightYieldCorrCountRanges = []          # for count in range not included in multipeak fit
    fitLightYieldCorrCountThres = 20
    fitInitRangeFromConfig      = False
    fitInitRangeSet             = {}
    if config_dict is not None: 
        if 'fitLightYieldCorrRanges' in config_dict: fitLightYieldCorrRanges = config_dict['fitLightYieldCorrRanges']
        fitLightYieldCorrPeaks = config_dict['fitLightYieldCorrPeaks']
        if 'lightYieldCorrCountRanges' in config_dict: lightYieldCorrCountRanges=config_dict['lightYieldCorrCountRanges']
        if 'fitLightYieldCorrCountThres' in config_dict:      
            fitLightYieldCorrCountThres = config_dict['fitLightYieldCorrCountThres']
        if 'fitInitRangeFromConfig' in config_dict: fitInitRangeFromConfig = config_dict['fitInitRangeFromConfig']
        if 'fitInitRangeSet'        in config_dict: fitInitRangeSet        = config_dict['fitInitRangeSet']
    if fitInitRangeFromConfig == True:
        if fitInitRangeSet == {}:
            raise AssertionError('peak_region_rate_estimation(): fitInitRangeFromConfig requires setting '+\
                                 'fitInitRangeSet in the configuration json file')
        for peakName in fitInitRangeSet: fitRangeSetPickle[peakName] = fitInitRangeSet[peakName]
        if verbosity >= 1:
            print('  fitInitRangeFromConfig is enabled, the new fitRangeSet is:\n')
            for peakName in fitRangeSetPickle: print('    fitRangeSet['+peakName+']  =', fitRangeSetPickle[peakName])
    if fitLightYieldCorrPeaks != None:
        for peakName in fitLightYieldCorrPeaks:
            if peakName not in fitRangeSetPickle:
                raise AssertionError('peak_region_rate_estimation(): '+peakName+' form fitLightYieldCorrPeaks not '
                                     'found in fitRangeSetPickle')
    #####################################################################################################################
    zName, r2name, r2Range = None, None, None
    if 'r2_naive'  in procData:
        zName   = 'area_fraction_top'
        r2name  = 'r2_naive'
        zRange  = [0.0, 1.0]
        r2Range = [0.0, 3.0E2]
    if 'r2_s1_cnn' in procData:
        zName  = 'z_s1_cnn'
        r2name = 'r2_s1_cnn'
        zRange  = [-150.0, 0.0]
        r2Range = [0.0, 6.0E3]
    r2_corrRange = [0.0, 6.0E3]
    if r2name is None:
        raise AssertionError('peak_region_rate_estimation(): either r2_naive or r2_s1_cnn are required for s1only '+\
                             'analysis')
    #####################################################################################################################
    mask0 = np.copy(mask_dict['mask_base'])
    mask1 = np.copy(mask_dict['mask_final_FV'])
    #####################################################################################################################
    if verbosity >= 1: print('  base data:')
    plot_name = 's1only_peak_region_rate_estimation_0base' 
    if save_plot == True:
        plt, ax = plotStandard_s1only_highA(procData, mask0, mask1, verbosity=verbosity,binN=generalBinN,\
                                            r2Key=r2name, r2Range=r2Range, zKey=zName, zRange=zRange)
        if plot_extra is not None: plot_extra(plot_name, ax, mask_dict) 
        plot_path = os.path.join(plot_dir, dateStr+output_name_tag+plot_name+'.'+plot_format)
        plt.savefig(plot_path, format=plot_format)
        if verbosity >= 1: print('    saving:', plot_path)
        plt.close('all')
    ###S1 Corr###############################################################################
    plot_name = 's1only_peak_region_rate_estimation_1s1corr'
    if save_plot == True:
        plt, ax = plotStandard_s1only_highA(procData, mask0, mask1, verbosity=verbosity, binN=generalBinN,\
                                            ampKey='s1_corrF', r2Key='r2_corr', r2Range=r2_corrRange, zKey=zName,\
                                            zRange=zRange)
        if plot_extra is not None: plot_extra(plot_name, ax, mask_dict)
        plot_path = os.path.join(plot_dir, dateStr+output_name_tag+plot_name+'.'+plot_format)
        plt.savefig(plot_path, format=plot_format)
        if verbosity >= 1: print('    saving:', plot_path)
        plt.close('all') 
    #####
    # Chunk Multi-Peak Rate with No Light Yeild Correction
    #####################################################################################################################
    plot_name = 's1only_peak_region_rate_estimation_2peak_fit0'
    if verbosity >= 1: print('  multi-peak fit, no LY corr:')
    if len(procData['run_number'][mask0]) == 0: 
        raise AssertionError('peak_region_rate_estimation(): mask0 cuts out all events')
    time_init = toBeginningOfDay(1.0E-9*procData['time'][mask0][0])
    time_fin  = 1.0E-9*procData['endtime'][mask0][-1]

    peak_amps, peak_width_percent = {}, {}
    for peakIdx, peakName in enumerate(peak_normParPickle):
        peak_amps[peakName]          = peak_normParPickle[peakName][0]
        peak_width_percent[peakName] = 100.0*peak_normParPickle[peakName][1]/peak_normParPickle[peakName][0]
    
    countRangeSet = {}
    for peakName in peak_normParPickle:
        fitMu  = peak_normParPickle[peakName][0]
        fitSig = peak_normParPickle[peakName][1]
        countRangeSet[peakName] = genCountRange(peakName, fitMu, fitSig, countRangeSigN)
    if verbosity >= 1:
        for peakName in peak_normParPickle: print('    countRangeSet['+peakName+'] =',  countRangeSet[peakName])
    ###########################################################
    runChunks, runModes = None, None
    if (config_dict is not None) and ('reference_runChunk_paths' in config_dict):
        reference_runChunk_paths = config_dict['reference_runChunk_paths']
        runChunks, runModes = [], []
        for reference_runChunk_path in reference_runChunk_paths:
            if os.path.isfile(reference_runChunk_path) == False:
                raise AssertionError('peak_region_rate_estimation(): pickle file:\n  '+reference_runChunk_path+\
                                     '\ndoes nto exist')
            runChunks_, runModes_ = genRunChunksWithRef(np.unique(procData['run_number']), reference_runChunk_path,\
                                                        verbosity=(verbosity-1))
            runChunks.extend(runChunks_)
            runModes .extend(runModes_)
    else:
        runChunks, runModes = genRunChunks(np.unique(procData['run_number']), runInfoData, combTimeThres,\
                                           checkTags=chunkCheckTags, verbosity=(verbosity-1))
    ###########################################################
    timestampRanges, liveTimes, dayArr = [], [], []
    eventDataS1ampArr, eventDataS1ampArr_multipeakFit, eventData_chunk = [], [], []
    eventCounts = {peakName:[] for peakName in peak_normParPickle}
    eventRates  = {peakName:[] for peakName in peak_normParPickle}
    if (deadtime_type is not None) and ((deadtime_type in runInfoData) == False): deadtime_type == None
    for runChunk in tqdm(runChunks, disable=(verbosity<=0)):
        timestampRangeChunk, liveTimeChunk = [], ufloat(0, 0)
        mask_run_sel = (procData['run_number'] == -1)
        for run_number in runChunk:      
            timestampRangeChunk.append([runInfoData['start_timestamp'][runInfoData['name']==run_number][0],
                                        runInfoData['end_timestamp'][runInfoData['name']==run_number][0]])
            liveTime  = pd.Timedelta(runInfoData['livetime'][runInfoData['name']==run_number][0]).total_seconds()
            if deadtime_type is not None: liveTime-=runInfoData[deadtime_type][runInfoData['name']==run_number][0]*1.0E-9
            else: warnings.warn('peak_region_rate_estimation(): deadtime plugin NOT implemented', Warning)
            liveTimeChunk += ufloat(liveTime, 0.5)          #NOTE: 0.5 gps uncertainty for livetime
            mask_run_sel |= (procData['run_number'] == run_number)
        timestampRanges.append(timestampRangeChunk)
        liveTimes      .append(liveTimeChunk)
        dayArr         .append((timestampRangeChunk[0][0] - time_init)/(3600.0*24.0))

        eventDataS1amp = np.array(procData['s1_corrF'][mask1 & mask_run_sel], dtype='float64')
        eventDataS1ampArr.append(eventDataS1amp)
        eventData_chunk.append({key:procData[key][mask1 & mask_run_sel] for key in procData})
        if 'mask_multipeakFit' in mask_dict:
            eventDataS1amp_multipeakFit = np.array(procData['s1_corrF'][mask_dict['mask_multipeakFit'] & mask_run_sel],\
                                                   dtype='float64')
            eventDataS1ampArr_multipeakFit.append(eventDataS1amp_multipeakFit)
        countDataSet = {}
        for peakName in peak_normParPickle:
            peakData = eventDataS1amp[(countRangeSet[peakName][0] <= eventDataS1amp) &
                                      (eventDataS1amp < countRangeSet[peakName][1])]
            countDataSet[peakName] = len(peakData)
            eventCounts[peakName].append(len(peakData))
        activityPars = [countDataSet, countRangeSet, peak_normParPickle, liveTimeChunk, fidMassPar['XeFidMass']]
        activity  = get_activity(peakName, *activityPars)
        activity *= fidMassPar['validS1S2Ratio']                            #NOTE: s1only rate correction from s1s2
        for peakName in peak_normParPickle: eventRates[peakName].append(activity)
    ###########################################################
    if save_plot == True:
        gridSpec = [2, 1]
        figSize, marginRatio = getSizeMargin(gridSpec, subplotSize=[18.0, 5.0], plotMargin=[0.08, 0.14, 0.04, 0.06])
        fig = plt.figure(figsize=figSize); fig.subplots_adjust(*marginRatio)
        gs = gridspec.GridSpec(*gridSpec)
        matplotlib.rc('xtick', labelsize=18)
        matplotlib.rc('ytick', labelsize=18)
        ax = []
        for axIdx in range(gs.nrows*gs.ncols):
            ax.append(fig.add_subplot(gs[axIdx]));
            ax[-1].ticklabel_format(style='sci', scilimits=(-2, 2))
            ax[-1].xaxis.set_tick_params(which='major', width=2.0, length=4.0)
            ax[-1].xaxis.set_tick_params(which='minor', width=1.0, length=2.0)
            ax[-1].xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
            ax[-1].yaxis.set_tick_params(which='major', width=2.0, length=4.0)
            ax[-1].yaxis.set_tick_params(which='minor', width=1.0, length=2.0)
            ax[-1].yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
        ######
        axIdx = 0
        hist = ax[axIdx].hist(procData['s1_corrF'][mask1], bins=s1_ampBinN, range=s1_ampRange, color='steelblue',\
                              linewidth=1)
        ax[axIdx].set_xlabel('s1_corrF', fontsize=18)
        ax[axIdx].set_ylabel('counts', fontsize=18)
        ax[axIdx].set_xlim(s1_ampRange)
        ax[axIdx].set_ylim(bottom=0.5)
        ax[axIdx].set_yscale('log')

        for peakName in peak_normParPickle:
            mask_isotope  = (countRangeSet[peakName][0] <= procData['s1_corrF']) 
            mask_isotope &= (procData['s1_corrF'] < countRangeSet[peakName][1])
            mask_isotope &= mask1
            hist = ax[axIdx].hist(procData['s1_corrF'][mask_isotope], bins=s1_ampBinN, range=s1_ampRange,
                                  color='gainsboro', linewidth=1, alpha=0.3)

        linePlot, linePlotText= {}, {}
        for peakName in peak_normParPickle:
            linePlot[peakName] = ax[axIdx].axvline(x=peak_amps[peakName], ymin=0.0, ymax=ax[axIdx].get_ylim()[1],
                                                   color=ALPHA_LINE_C[peakName], linestyle='dashed', linewidth=1.5)
            linePlotText[peakName] = peakName + ', w=' + scientificStr(peak_width_percent[peakName]) + '%'
        for peakName in ALPHA_LINE_E:
            if peakName not in peak_normParPickle:
                peak_amp = cali_factor*ALPHA_LINE_E[peakName]
                linePlot[peakName] = ax[axIdx].axvline(x=peak_amp, ymin=0.0, ymax=ax[axIdx].get_ylim()[1],
                                                       color=ALPHA_LINE_C[peakName], linestyle='dashed', linewidth=1.5)
                linePlotText[peakName] = peakName + ''
        legObj = ax[axIdx].legend([linePlot[peakName] for peakName in ALPHA_LINE_E],\
                                  [linePlotText[peakName] for peakName in ALPHA_LINE_E], loc='upper left', fontsize=16)
        ######
        axIdx = 1
        tracePlot, tracePlotText = {}, {}
        for peakName in peak_normParPickle:
            tracePlot[peakName] = ax[axIdx].errorbar(dayArr, [item.n for item in eventRates[peakName]], 
                                                             [item.s for item in eventRates[peakName]], 
                                                     color=ALPHA_LINE_C[peakName], fmt='s', **ERROR_PLOT_PARS)[0]
            tracePlotText[peakName]  = peakName
        setXaxisDateLabel(ax[axIdx], time_init, [time_init, time_fin])
        ax[axIdx].set_ylabel('activity concentration (uBq/kg_Xe)', fontsize=18)
        ax[axIdx].set_yscale('log')
        ax[axIdx].secondary_yaxis('right')
        if 'Po218' in tracePlotText: tracePlotText['Po218'] = tracePlotText['Po218'].replace('Po218', 'Po218/Bi212')
        if 'Bi212' in tracePlotText: tracePlotText['Bi212'] = tracePlotText['Bi212'].replace('Bi212', 'Po218/Bi212')
        legObj=ax[axIdx].legend([tracePlot[peakName]     for peakName in ALPHA_LINE_E if peakName in peak_normParPickle],
                                [tracePlotText[peakName] for peakName in ALPHA_LINE_E if peakName in peak_normParPickle],
                                loc='upper left', fontsize=10, markerscale=1.5)
        if plot_extra is not None: plot_extra(plot_name, ax, mask_dict) 
        plot_path = os.path.join(plot_dir, dateStr+output_name_tag+plot_name+'.'+plot_format)
        plt.savefig(plot_path, format=plot_format)
        if verbosity >= 1: print('    saving:', plot_path)
        plt.close('all')
    #####
    # Chunk Multi-peak Fit, Light Yeild Correction
    #####################################################################################################################
    plot_dir_name   = 's1only_peak_region_rate_estimation_3chunk_fit'
    plot_name_chunk = 's1only_peak_region_rate_estimation_3chunk_multipeak_fit'

    peakFitFunc = multiSigLinkedGausPDF_scipy
    output_path = os.path.join(output_dir, output_name_tag+'s1only_peak_region_rate_estimation.pkl')
    ###########################################################
    fitLightYieldCorrPeaks = [fitLightYieldCorrPeak for fitLightYieldCorrPeak in fitLightYieldCorrPeaks
                              if fitLightYieldCorrPeak in fitRangeSetPickle]
    if fitLightYieldCorrRanges is None: fitLightYieldCorrRanges = fitLightYieldCorrPeaks
    fitRanges_gainCorr, optFitPars, peakFitPars, relative_light_yields, corr_factors, chi2_ndofs = {}, [], [], [], [], []
    eventData = {}
    ###########################################################
    RnRatePickle, inputDictToExtend = None, None    
    inputDictToExtend = {'runChunks': runChunks, 'runModes': runModes, 'timestampRanges': timestampRanges,\
                         'liveTimes': liveTimes, 'eventCounts': eventCounts, 'eventRates': eventRates}
    RnRatePickle = loadPickleRateData(output_path, dayArr, time_init, inputDictToExtend=inputDictToExtend)
    if RnRatePickle is not None:
        if fitLightYieldCorrPeaks != RnRatePickle['fitLightYieldCorrPeaks']:
            raise AssertionError('peak_region_rate_estimation(): fitLightYieldCorrPeaks =', fitLightYieldCorrPeaks,\
                                 '\nwhile', output_path, 'provides fitLightYieldCorrPeaks =',\
                                 RnRatePickle['fitLightYieldCorrPeaks'])
        for peakName in fitLightYieldCorrPeaks:
            fitRanges_gainCorr[peakName] = copy.copy(RnRatePickle['fitRanges_gainCorr'][peakName])
        optFitPars            = copy.copy(RnRatePickle['optFitPars'])
        peakFitPars           = copy.copy(RnRatePickle['peakFitPars'])
        relative_light_yields = copy.copy(RnRatePickle['relative_light_yields'])
        corr_factors          = copy.copy(RnRatePickle['corr_factors'])
        chi2_ndofs            = copy.copy(RnRatePickle['chi2_ndofs'])
        if 'multipeak_chi2_ndof' in RnRatePickle['eventData']:
            eventData = copy.deepcopy(RnRatePickle['eventData'])
            timeMask = (eventData['run_number'] < runChunks[0][0])
            eventData = {key:eventData[key][timeMask] for key in eventData}
    ###########################################################
    if verbosity >= 1: 
        print('  fitLightYieldCorrPeaks:', fitLightYieldCorrPeaks)
        print('  chunk fits to process:', len(eventDataS1ampArr), '\n')
    if len(eventDataS1ampArr) == 0: raise AssertionError('peak_region_rate_estimation(): no run chunks to process')
    for chunkIdx, eventDataS1amp in enumerate(eventDataS1ampArr):
        eventDataS1ampC = np.array(eventDataS1amp, dtype='float64')
        eventDataS1ampFit = copy.deepcopy(eventDataS1ampC)
        if 'mask_multipeakFit' in mask_dict:
            eventDataS1ampFit = np.array(eventDataS1ampArr_multipeakFit[chunkIdx], dtype='float64')
        s1ampHist, s1ampX = np.histogram(eventDataS1ampFit, bins=s1_ampBinN, range=s1_ampRange)
        s1ampXcenter = 0.5*(s1ampX[1:] + s1ampX[:-1])
        s1ampX = s1ampX[:-1]
        s1ampErr = np.array([np.sqrt(count) if count != 0 else 1.0 for count in s1ampHist])   
        ###########################scaling for the best fit range
        scaleStep   = 0.02
        segIdxRange = [-25, 25]
        if config_dict is not None:
            if 'scaleStep'   in config_dict: scaleStep   = config_dict['scaleStep']
            if 'segIdxRange' in config_dict: segIdxRange = config_dict['segIdxRange']
        fitRangeChunk = {}
        optCount, optIdx, optIdx_sameN = -1, 0, 0
        if verbosity >= 3: print('  find peak fit range shift:')
        for segIdx in range(*segIdxRange, 1):
            countInRange = 0
            for peakName in fitLightYieldCorrRanges:
                fitRangeChunk[peakName] = np.array(fitRangeSetPickle[peakName])*(1.0 + segIdx*scaleStep)
                countInRange += len(eventDataS1ampFit[(fitRangeChunk[peakName][0] <= eventDataS1ampFit) &
                                                      (eventDataS1ampFit < fitRangeChunk[peakName][1])])
            if countInRange >= optCount:
                if countInRange == optCount: optIdx_sameN += 1
                else:                        optIdx_sameN  = 0
                optCount = np.copy(countInRange)
                optIdx   = np.copy(segIdx)
            if verbosity >= 3: print('   ', segIdx, segIdx*scaleStep, countInRange, optCount)
        optIdx -= int(optIdx_sameN/2)
        for peakName in fitLightYieldCorrPeaks:
            fitRangeChunk[peakName] = np.array(fitRangeSetPickle[peakName])*(1.0 + optIdx*scaleStep)
            if peakName not in fitRanges_gainCorr: fitRanges_gainCorr[peakName] = []
            fitRanges_gainCorr[peakName].append(fitRangeChunk[peakName])
        ##################################################################
        if verbosity >= 1: 
            print(chunkIdx, '/', len(eventDataS1ampArr), runChunks[chunkIdx])
            print('  day:     ', dayArr[chunkIdx])
            print('  eventN:  ', len(eventDataS1ampC))
            print('  livetime:', liveTimes[chunkIdx])
        countInRange = 0
        for peakName in fitLightYieldCorrPeaks:
            if verbosity >= 1: print('   ', peakName, fitRangeSetPickle[peakName], fitRangeChunk[peakName])
            countInRange += len(eventDataS1ampFit[(fitRangeChunk[peakName][0] <= eventDataS1ampFit) &
                                                  (eventDataS1ampFit < fitRangeChunk[peakName][1])])
        if verbosity >= 2: print('    optIdx =', optIdx, ', countInRange =', countInRange)
        ##################################################################
        fitPars, fitErrs, fitCov = None, None, None
        if countInRange >= fitLightYieldCorrCountThres:
            fitDataSet, initPars = {}, []
            fitRangeMask = np.full(len(s1ampX), False)
            for peakName in fitLightYieldCorrPeaks:
                fitDataSet[peakName] = eventDataS1ampFit[(fitRangeChunk[peakName][0] <= eventDataS1ampFit)&
                                                         (eventDataS1ampFit < fitRangeChunk[peakName][1])] 
                fitRangeMask |= ((fitRangeChunk[peakName][0] < s1ampX) & (s1ampX < fitRangeChunk[peakName][1]))
                initPars.append(len(fitDataSet[peakName])*(s1ampX[1]-s1ampX[0]))
                if len(fitDataSet[peakName]) != 0: initPars.append(np.mean(fitDataSet[peakName]))
                else:                              initPars.append(np.mean(fitRangeChunk[peakName]))
            if fitLightYieldCorrPeaks[0] in peak_normParPickle:
                initPars.append(peak_normParPickle[fitLightYieldCorrPeaks[0]][1])
            else:
                initPars.append(np.std(fitDataSet[fitLightYieldCorrPeaks[0]]))

            costFunc = iminuit.cost.LeastSquares(s1ampX, s1ampHist, s1ampErr, peakFitFunc)
            costFunc.mask = fitRangeMask
            objMinuit = iminuit.Minuit(costFunc, initPars)
            for peakIdx, peakName in enumerate(fitLightYieldCorrPeaks):
                objMinuit.limits[2*peakIdx]   = [0.1*initPars[2*peakIdx], 10.0*initPars[2*peakIdx]]
                objMinuit.limits[2*peakIdx+1] = fitRangeChunk[peakName]
                if initPars[2*peakIdx] == 0:
                    objMinuit.fixed[2*peakIdx]   = True
                    objMinuit.fixed[2*peakIdx+1] = True
            objMinuit.limits[-1] = [0.8*initPars[-1], 1.2*initPars[-1]]         #NOTE the fit limit for the peak width
            if verbosity >= 2: print('    initial values:', initPars, '\n    fit limits:', objMinuit.limits)
            optResult = objMinuit.migrad()
            try:
                fitPars = [val for val in optResult.values]
                fitErrs = [err for err in optResult.errors]
                fitCov  = [row for row in optResult.covariance]
            except:
                pass
            if verbosity >= 1: print('    multi-peak fit:\n', optResult)
            ############################
            if optResult.valid == True:
                peak_normPar_chunk, peak_locs = {}, []
                for peakIdx, peakName in enumerate(fitLightYieldCorrPeaks):
                    peak_normSig  = ufloat(fitPars[-1], fitErrs[-1])
                    peak_normSig *= ufloat(fitPars[2*peakIdx+1], fitErrs[2*peakIdx+1])/ufloat(fitPars[1],fitErrs[1])
                    peak_normPar_chunk[peakName] = [ufloat(fitPars[2*peakIdx+1], fitErrs[2*peakIdx+1]),peak_normSig]
                    peak_locs.append(ufloat(fitPars[2*peakIdx+1], fitErrs[2*peakIdx+1])/ALPHA_LINE_E[peakName])
                peak_loc_ave = weightedAverage(peak_locs)
                optFitPars           .append({'fitPars': fitPars, 'fitErrs': fitErrs, 'fitCov': fitCov})
                peakFitPars          .append(peak_normPar_chunk)
                relative_light_yields.append(peak_loc_ave/cali_factor)
                corr_factors         .append(peak_loc_ave)
                chi2_ndofs           .append(optResult.fval/optResult.ndof)
                if verbosity >= 1: 
                    print('    relative_light_yield =', relative_light_yields[-1])
                    print('    chi2_ndof =', chi2_ndofs[-1])
                fitX     = np.linspace(*s1_ampRange, 10*s1_ampBinN+1)[:-1]
                fitY     = peakFitFunc(fitX, fitPars)
                fitYinit = peakFitFunc(fitX, initPars)
                fitResidual = np.array([np.sqrt(iminuit.cost.poisson_chi2(np.array(\
                    [1.0*histYval]), np.array([1.0*fitYval])))\
                     for histYval, fitYval in zip(s1ampHist, peakFitFunc(s1ampX, fitPars))])
            else:
                optFitPars           .append({'fitPars': None, 'fitErrs': None, 'fitCov': None})
                peakFitPars          .append({peakName:None for peakName in fitLightYieldCorrPeaks})
                relative_light_yields.append(ufloat(-1, 0))
                corr_factors         .append(ufloat(-1, 0))
                chi2_ndofs           .append(-1)
        else:
            optFitPars           .append({'fitPars': None, 'fitErrs': None, 'fitCov': None})
            peakFitPars          .append({peakName:None for peakName in fitLightYieldCorrPeaks})
            relative_light_yields.append(ufloat(-1, 0))
            corr_factors         .append(ufloat(-1, 0))
            chi2_ndofs           .append(-1)
        ###########################################################
        eventData_chunk[chunkIdx]['relative_light_yield'] = 0*eventData_chunk[chunkIdx]['s1_corrF'] +\
                                                            relative_light_yields[-1].n
        eventData_chunk[chunkIdx]['multipeak_chi2_ndof']  = 0*eventData_chunk[chunkIdx]['s1_corrF'] + chi2_ndofs[-1]
        eventData_chunk[chunkIdx]['s1_corrF_LYcorr']      = eventData_chunk[chunkIdx]['s1_corrF']/\
                                                            relative_light_yields[-1].n
        if len(eventData) == 0: eventData = copy.deepcopy(eventData_chunk[chunkIdx])
        else:
            for key in eventData: 
                if key in eventData_chunk[chunkIdx]:
                    eventData[key] = np.append(eventData[key], eventData_chunk[chunkIdx][key]) 
                else:
                    eventData[key] =np.append(eventData[key],np.nan*np.empty(len(eventData_chunk[chunkIdx]['s1_corrF'])))
        ###########################################################
        if save_plot == True:
            gridSpec = [1, 1]
            figSize, marginRatio = getSizeMargin(gridSpec, subplotSize=[18.0, 5.0], plotMargin=[0.06, 0.16, 0.02, 0.12])
            fig = plt.figure(figsize=figSize); fig.subplots_adjust(*marginRatio)
            gs = gridspec.GridSpec(*gridSpec)
            matplotlib.rc('xtick', labelsize=18)
            matplotlib.rc('ytick', labelsize=18)
            ax = []
            for axIdx in range(gs.nrows*gs.ncols):
                ax.append(fig.add_subplot(gs[axIdx]));
                ax[-1].ticklabel_format(style='sci', scilimits=(-2, 2))
                ax[-1].xaxis.set_tick_params(width=1.5, length=4)
                ax[-1].xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
                ax[-1].yaxis.set_tick_params(width=1.5, length=4)
                ax[-1].yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
            ######
            axIdx = 0
            hist = ax[axIdx].hist(eventDataS1ampC,   bins=s1_ampBinN, range=s1_ampRange, color='steelblue', linewidth=1,\
                                  alpha=0.5)
            hist = ax[axIdx].hist(eventDataS1ampFit, bins=s1_ampBinN, range=s1_ampRange, color='steelblue', linewidth=1)
            runChunkTitle = [str(runChunk[runIdx:runIdx          + 12])
                             for runIdx in range(0, len(runChunk), 12)]
            dateInit = timestampToDate(timestampRanges[chunkIdx][0][0]).split('\n')
            dateInit = dateInit[0] + '_' + dateInit[1]
            dateFin  = timestampToDate(timestampRanges[chunkIdx][-1][1]).split('\n')
            dateFin  = dateFin[0]  + '_' + dateFin[1]
            titleStr = ''
            titleStr += dateInit + ' to ' + dateFin + '\n'
            # titleStr += 'Run: ' + '\n'.join(runChunkTitle) + '\n'
            titleStr += 'Run: ' + str(runChunks[chunkIdx][0]) + '-' + str(runChunks[chunkIdx][-1]) 
            titleStr += ' ' + runModes[chunkIdx][0] + ', '
            titleStr += 'Count: ' + str(len(eventDataS1ampC)) + ', '
            titleStr += 'Duration: ' + str(liveTimes[chunkIdx]) + ' s'
            ax[axIdx].set_title(titleStr, fontsize=16, y=1.03)
            ax[axIdx].get_xaxis().set_visible(False)
            ax[axIdx].set_xlabel('s1_corrF', fontsize=18)
            ax[axIdx].set_ylabel('counts', fontsize=18)
            ax[axIdx].set_xlim(s1_ampRange)
            ax[axIdx].set_ylim(bottom=0.5)
            ax[axIdx].set_yscale('log')
    
            divider = make_axes_locatable(ax[axIdx])
            axRes = divider.append_axes('bottom', size='20%', pad='2%')
            ax[axIdx].figure.add_axes(axRes)
            axRes.set_xlim(s1_ampRange)
            axRes.set_xlabel('s1_corrF', fontsize=18)
            axRes.set_ylabel('fit\nresidual', fontsize=18)
            axRes.ticklabel_format(style='sci', scilimits=(-2, 2))
            axRes.xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
            axRes.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())    
            if fitCov is not None:
                linePlot, linePlotText= {}, {}      
                for peakName in fitLightYieldCorrPeaks:
                    fitMask = ((fitRangeChunk[peakName][0] < fitX) & (fitX < fitRangeChunk[peakName][1]))
                    # ax[axIdx].plot(fitX[fitMask], fitYinit[fitMask], color='orange', linewidth=1.0, alpha=1.0,zorder=2)
                    ax[axIdx].plot(fitX[fitMask], fitY[fitMask], color='red', linewidth=1.5, alpha=0.8, zorder=3)
                    linePlot[peakName] = ax[axIdx].axvline(x=peak_normPar_chunk[peakName][0].n, 
                                                           color=ALPHA_LINE_C[peakName],linestyle='dashed',linewidth=1.5)
                    linePlotText[peakName] = peakName + ''
                for peakName in ALPHA_LINE_E:
                    if peakName not in fitLightYieldCorrPeaks:
                        linePlot[peakName] = ax[axIdx].axvline(x=corr_factors[-1].n*ALPHA_LINE_E[peakName],
                                                               color=ALPHA_LINE_C[peakName], linestyle='dashed', 
                                                               linewidth=1.5)
                        linePlotText[peakName] = peakName + ''        
                legObj = ax[axIdx].legend([linePlot[peakName] for peakName in ALPHA_LINE_E],\
                                          [linePlotText[peakName] for peakName in ALPHA_LINE_E],\
                                          title=('chi2/ndof=' + scientificStr(chi2_ndofs[-1])),\
                                          loc='upper left', fontsize=16)
                # axRes.errorbar(fitX[costFunc.mask],fitResidual[costFunc.mask],yerr=np.sqrt(s1ampHist[costFunc.mask]),
                #                color='black', fmt='s', **ERROR_PLOT_PARS)
                for peakName in fitLightYieldCorrPeaks:
                    fitMask = ((fitRangeChunk[peakName][0] <= s1ampX) & (s1ampX < fitRangeChunk[peakName][1]))
                    axRes.plot(s1ampX[fitMask], fitResidual[fitMask], color='black')
                axRes.axhline(y=1.0, color='grey', linestyle='dashed', zorder=0)
                ylimRes = list(axRes.get_ylim())
                ylimRes[0] = 0.0
                if ylimRes[1] <= 5.0: ylimRes[1] = 5.0
                else:                 ylimRes[1] = roundSig(ylimRes[1], 2)
                axRes.set_ylim(*ylimRes)  
            for peakName in lightYieldCorrCountRanges:
                countRangeNoFit = np.array(lightYieldCorrCountRanges[peakName])*relative_light_yields[-1].n
                rangeMaskNoFit = (countRangeNoFit[0] <= eventDataS1ampFit) & (eventDataS1ampFit < countRangeNoFit[1])
                ax[axIdx].hist(eventDataS1ampFit[rangeMaskNoFit], bins=s1_ampBinN, range=s1_ampRange,\
                               color=ALPHA_LINE_C[peakName], linewidth=1)
                axRes.hlines(y=1.0, xmin=countRangeNoFit[0], xmax=countRangeNoFit[1], color=ALPHA_LINE_C[peakName]) 
            plt.setp(legObj.get_title(), fontsize=16)
            if plot_extra is not None: plot_extra(plot_name_chunk, ax, mask_dict) 
            plot_subdir = dateStr + output_name_tag + plot_dir_name
            pathlib.Path(os.path.join(plot_dir, plot_subdir)).mkdir(parents=True, exist_ok=True)
            plot_path = os.path.join(plot_dir, plot_subdir, dateStr+plot_name_chunk+\
                                     str(runChunks[chunkIdx][0])+'.'+plot_format)
            plt.savefig(plot_path, format=plot_format)
            if verbosity >= 1: print('    saving:', plot_path)
            plt.close('all')
    #####
    # Rate Evalution
    #####################################################################################################################
    plot_name = 's1only_peak_region_rate_estimation_4peak_rate'
    if verbosity >= 1: print('  rate evalution:')    
    ###########################################################
    countRanges_gainCorr = {}
    eventCounts_gainCorr = {peakName:[] for peakName in fitLightYieldCorrPeaks}
    eventRates_gainCorr  = {peakName:[] for peakName in fitLightYieldCorrPeaks}
    mask_isotope = {}
    chunkIdxOffset = 0
    for peakName in fitLightYieldCorrPeaks: mask_isotope[peakName] = np.copy(mask1)
    ###########################################################
    if RnRatePickle is not None:
        for peakName in RnRatePickle['countRanges_gainCorr']:
            countRanges_gainCorr[peakName] = copy.copy(RnRatePickle['countRanges_gainCorr'][peakName])
            if chunkIdxOffset == 0: chunkIdxOffset = len(countRanges_gainCorr[peakName]) 
        for peakName in RnRatePickle['eventCounts_gainCorr']:
            eventCounts_gainCorr[peakName] = copy.copy(RnRatePickle['eventCounts_gainCorr'][peakName]) 
            eventRates_gainCorr[peakName]  = copy.copy(RnRatePickle['eventRates_gainCorr'][peakName])
    for chunkIdx, eventDataS1amp in enumerate(eventDataS1ampArr):  
        if verbosity >= 1: print(chunkIdx, runChunks[chunkIdx])
        mask_run_sel = np.array([(run_number in runChunks[chunkIdx]) for run_number in procData['run_number']])
        eventDataS1ampC = np.array(eventDataS1amp, dtype='float64')
        peak_sel_masks = {}
        if chi2_ndofs[chunkIdxOffset+chunkIdx] != -1:
            countRangeSet, peakMuSigSet = {}, {}
            for peakName in fitLightYieldCorrPeaks:
                fitMu  = peakFitPars[chunkIdxOffset+chunkIdx][peakName][0].n
                fitSig = peakFitPars[chunkIdxOffset+chunkIdx][peakName][1].n
                countRangeSet[peakName] = genCountRange(peakName, fitMu, fitSig, countRangeSigN)
                peakMuSigSet[peakName] = np.array(peakFitPars[chunkIdxOffset+chunkIdx][peakName])
            for peakName in lightYieldCorrCountRanges:
                countRangeSet[peakName] = np.array(lightYieldCorrCountRanges[peakName])*relative_light_yields[chunkIdx].n
                peakMuSigSet[peakName]  = None
                if peakName not in eventCounts_gainCorr:
                    eventCounts_gainCorr[peakName] = []
                    eventRates_gainCorr [peakName] = []
                if peakName not in mask_isotope: mask_isotope[peakName] = copy.deepcopy(mask1)
            for peakName in countRangeSet:
                if peakName not in countRanges_gainCorr: countRanges_gainCorr[peakName] = []
                countRanges_gainCorr[peakName].append(np.array(countRangeSet[peakName]))
                peakData = eventDataS1ampC[(countRangeSet[peakName][0] <= eventDataS1ampC)&
                                           (eventDataS1ampC < countRangeSet[peakName][1])]
                eventCounts_gainCorr[peakName].append(len(peakData))
            countDataSet = {peakName:eventCounts_gainCorr[peakName][-1] for peakName in eventCounts_gainCorr}
            activityPars = [countDataSet, countRangeSet, peakMuSigSet, liveTimes[chunkIdx], fidMassPar['XeFidMass']]
            for peakName in countRangeSet: 
                activity  = get_activity(peakName, *activityPars)
                #NOTE: s1only rate correction from s1s2
                activity *= fidMassPar['validS1S2Ratio']
                eventRates_gainCorr[peakName].append(activity)
                peak_sel_masks[peakName] = mask_run_sel & (countRangeSet[peakName][0] <= procData['s1_corrF']) &\
                                                          (procData['s1_corrF'] < countRangeSet[peakName][1]) 
                mask_isotope[peakName] &= (np.logical_not(mask_run_sel) | peak_sel_masks[peakName])
            if verbosity >= 1:
                print('    liveTime  =', liveTimes[chunkIdx], 's')
                print('    XeFidMass =', fidMassPar['XeFidMass'], 'kg')
                for peakName in peakMuSigSet:
                    print('    peakMuSigSet['+peakName+'] =', peakMuSigSet[peakName])
                for peakName in countRanges_gainCorr:
                    print('    countRanges_gainCorr['+peakName+'] =', countRanges_gainCorr[peakName][-1])
                for peakName in eventCounts_gainCorr:
                    print('    eventCounts_gainCorr['+peakName+'] =', eventCounts_gainCorr[peakName][-1])
                for peakName in eventRates_gainCorr:
                    print('    eventRates_gainCorr['+peakName+'] =',  eventRates_gainCorr[peakName][-1])
        else:
            for peakName in countRanges_gainCorr:
                countRanges_gainCorr[peakName].append(-1)
            for peakName in fitLightYieldCorrPeaks:
                eventCounts_gainCorr[peakName].append(0)
                eventRates_gainCorr [peakName].append(ufloat(-1.0, 0.0))
                mask_isotope[peakName] &= np.logical_not(mask_run_sel)
            for peakName in lightYieldCorrCountRanges:
                if peakName not in eventCounts_gainCorr:
                    eventCounts_gainCorr[peakName] = []
                    eventRates_gainCorr [peakName] = []
                    mask_isotope        [peakName] = copy.deepcopy(mask1)
                eventCounts_gainCorr[peakName].append(0)
                eventRates_gainCorr [peakName].append(ufloat(-1.0, 0.0))
                mask_isotope[peakName] &= np.logical_not(mask_run_sel)
    if RnRatePickle is not None:
        time_init = copy.copy(RnRatePickle['time_init'])
        runChunks       = copy.copy(RnRatePickle['runChunks'])
        runModes        = copy.copy(RnRatePickle['runModes'])
        timestampRanges = copy.copy(RnRatePickle['timestampRanges'])
        liveTimes       = copy.copy(RnRatePickle['liveTimes'])
        dayArr          = copy.copy(RnRatePickle['dayArr'])
        eventCounts = copy.copy(RnRatePickle['eventCounts'])
        eventRates  = copy.copy(RnRatePickle['eventRates'])
    ###########################################################
    timeDisplay = genTimeDisplay(dateDisplayRange, time_init, time_fin)
    eventDisplayMask  = (chi2DisplayRange[0] <= eventData['multipeak_chi2_ndof'])
    eventDisplayMask &= (eventData['multipeak_chi2_ndof'] < chi2DisplayRange[1])
    eventDisplayMask &= (timeDisplay[0] <= eventData['time']*1.0E-9)
    eventDisplayMask &= (eventData['endtime']*1.0E-9 < timeDisplay[1])
    if verbosity >= 1:  
        print('    [time_init, time_fin] =', [time_init, time_fin], ', timeDisplay =', timeDisplay)
        for peakName in mask_isotope: 
            print('    mask_isotope[\''+peakName+'\'],', np.sum(mask_isotope[peakName]))
    ###########################################################
    if save_plot == True:
        gridSpec = [1, 1]
        figSize, marginRatio = getSizeMargin(gridSpec, subplotSize=[18.0, 5.0], plotMargin=[0.06, 0.14, 0.02, 0.08])
        fig = plt.figure(figsize=figSize); fig.subplots_adjust(*marginRatio)
        gs = gridspec.GridSpec(*gridSpec)
        matplotlib.rc('xtick', labelsize=18)
        matplotlib.rc('ytick', labelsize=18)
        ax = []
        for axIdx in range(gs.nrows*gs.ncols):
            ax.append(fig.add_subplot(gs[axIdx]));
            ax[-1].ticklabel_format(style='sci', scilimits=(-2, 2))
            ax[-1].xaxis.set_tick_params(width=1.5, length=4)
            ax[-1].xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
            ax[-1].yaxis.set_tick_params(width=1.5, length=4)
            ax[-1].yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
        ######
        axIdx = 0
        eventDataToDisplay = eventData['s1_corrF_LYcorr'][eventDisplayMask]
        hist = ax[axIdx].hist(eventDataToDisplay, bins=s1_ampBinN, range=s1_ampRange, color='steelblue', linewidth=1)
        for peakName in lightYieldCorrCountRanges:
            countRangeNoFit = np.array(lightYieldCorrCountRanges[peakName])*1.0
            rangeMaskNoFit = (countRangeNoFit[0] <= eventDataToDisplay) & (eventDataToDisplay < countRangeNoFit[1])
            ax[axIdx].hist(eventDataToDisplay[rangeMaskNoFit], bins=s1_ampBinN, range=s1_ampRange,\
                           color=ALPHA_LINE_C[peakName], linewidth=1)
        ax[axIdx].set_xlabel('s1_corr3', fontsize=18)
        ax[axIdx].set_ylabel('counts', fontsize=18)
        ax[axIdx].set_xlim(s1_ampRange)
        ax[axIdx].set_ylim(bottom=0.5)
        ax[axIdx].set_yscale('log')

        linePlot, linePlotText= {}, {}
        for peakName in ALPHA_LINE_E:
            linePlot[peakName] = ax[axIdx].axvline(x=cali_factor*ALPHA_LINE_E[peakName], color=ALPHA_LINE_C[peakName], 
                                                   linestyle='dashed', linewidth=1.5)
            linePlotText[peakName] = peakName + ''
        legObj = ax[axIdx].legend([linePlot[peakName]     for peakName in ALPHA_LINE_E],\
                                  [linePlotText[peakName] for peakName in ALPHA_LINE_E], loc='upper left', fontsize=16)
        if plot_extra is not None: plot_extra(plot_name, ax, mask_dict) 
        plot_path = os.path.join(plot_dir, dateStr+output_name_tag+plot_name+'_.'+plot_format)
        plt.savefig(plot_path, format=plot_format)
        if verbosity >= 1: print('    saving:', plot_path)
        plt.close('all')
        ###########################################################
        gridSpec = [4, 1]
        figSize, marginRatio = getSizeMargin(gridSpec, subplotSize=[18.0, 4.0], plotMargin=[0.08, 0.16, 0.06, 0.06])
        fig = plt.figure(figsize=figSize); fig.subplots_adjust(*marginRatio)
        gs = gridspec.GridSpec(*gridSpec)
        matplotlib.rc('xtick', labelsize=18)
        matplotlib.rc('ytick', labelsize=18)
        ax = []
        for axIdx in range(gs.nrows*gs.ncols):
            ax.append(fig.add_subplot(gs[axIdx]));
            ax[-1].ticklabel_format(style='sci', scilimits=(-2, 2))
            ax[-1].xaxis.set_tick_params(which='major', width=2.0, length=4.0)
            ax[-1].xaxis.set_tick_params(which='minor', width=1.0, length=2.0)
            ax[-1].xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
            ax[-1].yaxis.set_tick_params(which='major', width=2.0, length=4.0)
            ax[-1].yaxis.set_tick_params(which='minor', width=1.0, length=2.0)
            ax[-1].yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
        ######
        validMask = (chi2DisplayRange[0] <= np.array(chi2_ndofs)) & (np.array(chi2_ndofs) < chi2DisplayRange[1])
        ######
        axIdx = 0
        tracePlot, tracePlotText = {}, {}
        for peakName in eventRates_gainCorr:
            tracePlot[peakName] = \
              ax[axIdx].errorbar(np.array(dayArr)[validMask],\
                                 [item.n for item in np.array(eventRates_gainCorr[peakName])[validMask]],\
                                 [item.s for item in np.array(eventRates_gainCorr[peakName])[validMask]],\
                                 color=ALPHA_LINE_C[peakName], fmt='s', **ERROR_PLOT_PARS)[0]
            tracePlotText[peakName] = peakName
        setXaxisDateLabel(ax[axIdx], time_init, timeDisplay, displayOff=True)
        ax[axIdx].set_ylabel('activity concentration\n(uBq/kg_Xe)', fontsize=18)
        ax[axIdx].set_yscale('log')
        ax[axIdx].secondary_yaxis('right')
        if 'Po218' in tracePlotText: tracePlotText['Po218'] = tracePlotText['Po218'].replace('Po218', 'Po218/Bi212')
        legObj = ax[axIdx].legend([tracePlot[peakName]     for peakName in eventRates_gainCorr\
                                   if (peakName not in peakDisplaySubset)],\
                                  [tracePlotText[peakName] for peakName in eventRates_gainCorr\
                                   if (peakName not in peakDisplaySubset)],\
                                  loc='upper left', fontsize=16, markerscale=1.5)
        ######
        axIdx = 1
        subTracePlot, subTracePlotText = {}, {}
        for peakName in peakDisplaySubset:
            subTracePlot[peakName] = \
              ax[axIdx].errorbar(np.array(dayArr)[validMask],\
                                 [item.n for item in np.array(eventRates_gainCorr[peakName])[validMask]],
                                 [item.s for item in np.array(eventRates_gainCorr[peakName])[validMask]],
                                 color=ALPHA_LINE_C[peakName], fmt='s', **ERROR_PLOT_PARS)[0]
            subTracePlotText[peakName] = peakName
        setXaxisDateLabel(ax[axIdx], time_init, timeDisplay, displayOff=True)
        ax[axIdx].set_ylabel('activity concentration\n(uBq/kg_Xe)', fontsize=18)
        ax[axIdx].set_yscale('log')
        ax[axIdx].secondary_yaxis('right')
        legObj = ax[axIdx].legend([subTracePlot[peakName]     for peakName in peakDisplaySubset],\
                                  [subTracePlotText[peakName] for peakName in peakDisplaySubset],\
                                  loc='upper right', fontsize=14, markerscale=1.5)
        ######
        axIdx = 2
        ax[axIdx].errorbar(np.array(dayArr)[validMask], [item.n for item in np.array(relative_light_yields)[validMask]], 
                                                        [item.s for item in np.array(relative_light_yields)[validMask]], 
                           color='blue', fmt='s', **ERROR_PLOT_PARS)[0]
        setXaxisDateLabel(ax[axIdx], time_init, timeDisplay, displayOff=True)
        ax[axIdx].set_ylabel('relative light yield', fontsize=18)
        ax[axIdx].secondary_yaxis('right')
        ######
        axIdx = 3
        ax[axIdx].scatter(np.array(dayArr)[validMask], np.array(chi2_ndofs)[validMask], color='blue')
        setXaxisDateLabel(ax[axIdx], time_init, timeDisplay)
        ax[axIdx].set_ylabel('multi-peak fit chi2/ndof', fontsize=18)
        ax[axIdx].set_ylim(bottom=0.0)
        ax[axIdx].secondary_yaxis('right')
        ax[axIdx].axhline(y=1.0, color='grey', linestyle='dashed', zorder=0)
        if plot_extra is not None: plot_extra(plot_name, ax, mask_dict, time_init=time_init) 
        plot_path = os.path.join(plot_dir, dateStr+output_name_tag+plot_name+'.'+plot_format)
        plt.savefig(plot_path, format=plot_format)
        if verbosity >= 1: print('    saving:', plot_path)
        plt.close('all')
    #####################################################################################################################
    rateDict = {}
    if time_init is None:
        time_init = 1.0E-9*procData['time'][mask_final_FV][0]
        time_init = toBeginningOfDay(time_init)
    rateDict['time_init'] = time_init
    rateDict['time_fin']  = time_fin
    rateDict['runChunks']       = runChunks
    rateDict['runModes']        = runModes
    rateDict['timestampRanges'] = timestampRanges
    rateDict['liveTimes']       = liveTimes
    rateDict['dayArr']          = dayArr
    rateDict['dayErr']          = []
    for timestampRangeChunk in timestampRanges:
        rateDict['dayErr'].append(0.5*(timestampRangeChunk[-1][1] - timestampRangeChunk[0][0])/(24*3600.0))
    rateDict['eventData'] = eventData
    rateDict['eventCounts'] = eventCounts
    rateDict['eventRates']  = eventRates
    rateDict['cali_factor_ref']        = cali_factor
    rateDict['fitLightYieldCorrPeaks'] = fitLightYieldCorrPeaks
    rateDict['fitRanges_gainCorr']     = fitRanges_gainCorr
    rateDict['optFitPars']             = optFitPars
    rateDict['peakFitPars']            = peakFitPars
    rateDict['relative_light_yields']  = relative_light_yields
    rateDict['corr_factors']           = corr_factors
    rateDict['chi2_ndofs']             = chi2_ndofs
    rateDict['countRanges_gainCorr'] = countRanges_gainCorr
    rateDict['eventCounts_gainCorr'] = eventCounts_gainCorr
    rateDict['eventRates_gainCorr']  = eventRates_gainCorr

    correlationFactor = fidMassPar['validS1S2Ratio']/fidMassPar['XeFidMass']
    rateDict['eventRates_gainCorr_noCorrErr'] = {}
    for peakName in rateDict['eventRates_gainCorr']:
        rateDict['eventRates_gainCorr_noCorrErr'][peakName] = []
        for eventRate in rateDict['eventRates_gainCorr'][peakName]:
            eventRate_noCorrErr = ufloat(-1, 0)
            if eventRate > 0: 
                eventRate_noCorrErr = removeMultiplicationErr(eventRate, correlationFactor)*correlationFactor.n
            rateDict['eventRates_gainCorr_noCorrErr'][peakName].append(eventRate_noCorrErr)
    rateDict['correction_derivation_parameters']  = corrParPickle
    rateDict['estimate_fiducial_mass_parameters'] = fidMassPar
    rateDict['parameters'] = {}
    for key in inputPars: rateDict['parameters'][key] = inputPars[key]
    for key in mask_dict:
        if ('mask' not in key) and ('Xe' not in key): rateDict['parameters'][key] = mask_dict[key]
    rateDict['parameters']['XeDensity_kg/cm^3'] = fidMassPar['XeDensity']
    rateDict['parameters']['XeFidVol_ref_cm^3'] = fidMassPar['XeFidVol_ref']
    rateDict['parameters']['XeFidMass_ref_kg']  = fidMassPar['XeFidMass_ref']
    rateDict['parameters']['XeFidVol_cm^3']     = fidMassPar['XeFidVol']
    rateDict['parameters']['XeFidMass_kg']      = fidMassPar['XeFidMass']
    with open(output_path, 'wb') as handle:
        pickle.dump(rateDict, handle, protocol=pickle.HIGHEST_PROTOCOL)
    if verbosity >=1: print('    saving:', output_path)

    rateDictLite = ufloatToErrDict(rateDict)
    del rateDictLite['eventData']
    with open(output_path.replace('.pkl', '_lite.pkl'), 'wb') as handle:
        pickle.dump(rateDictLite, handle, protocol=pickle.HIGHEST_PROTOCOL)
    if verbosity >=1: print('    saving:', output_path.replace('.pkl', '_lite.pkl'))
    #####################################################################################################################
    ### fiducial volume consistency check
    if (fidMassPar['S1AFT_cutRange']   != mask_dict['S1AFT_cutRange']) or\
       (fidMassPar['r2_corr_cutRange'] != mask_dict['r2_corr_cutRange']):
        warnings.warn('peak_region_rate_estimation(): s1only FV does NOT match that from fidMassPar', Warning)
        print('    fidMassPar[\'S1AFT_cutRange\'] =', fidMassPar['S1AFT_cutRange'])
        print('    fidMassPar[\'r2_corr_cutRange\'] =', fidMassPar['r2_corr_cutRange'])
        print('    mask_dict[\'S1AFT_cutRange\'] =', mask_dict['S1AFT_cutRange'])
        print('    mask_dict[\'r2_corr_cutRange\'] =', mask_dict['r2_corr_cutRange'])

#########################################################################################################################
def get_activity_surface_count(peakName, countDataSet, countRangeSet, peakMuSigSet, dataTime, XeFidArea):
    activity = ufloat(countDataSet[peakName], np.sqrt(countDataSet[peakName]))
    if peakName in peakMuSigSet:
        activity /= (gausCDF_err(countRangeSet[peakName][1], 1.0, *peakMuSigSet[peakName])
                    -gausCDF_err(countRangeSet[peakName][0], 1.0, *peakMuSigSet[peakName]))
    activity *= 1.0E3/((XeFidArea*1.0E-4)*dataTime)     #mBq/m^2
    return activity
def get_activity_surface_fit(peakName, peakFitPars, xBinSize, dataTime, XeFidArea):
    activity  = peakFitPars[peakName][0]/xBinSize
    activity *= 1.0E3/((XeFidArea*1.0E-4)*dataTime)     #mBq/m^2
    return activity
def surface_rate_estimation(runInfoData, procData, mask_dict, s1corr_paths, fidMass_path, output_dir=None,plot_dir=None,\
                            save_plot=None, plot_format=None, output_name_tag=None, config_dict=None, plot_extra=None,\
                            verbosity=None):
    print('--------------------------------passing surface_rate_estimation--------------------------------')
    pass
    '''
    if isinstance(runInfoData, pd.DataFrame) == True:
        runInfoData = runInfoData.to_dict('series')
        for key in runInfoData: runInfoData[key] = np.array(runInfoData[key])
    if isinstance(procData, pd.DataFrame) == True:
        procData = procData.to_dict('series')
        for key in procData: procData[key] = np.array(procData[key])
    if output_dir      is None: output_dir      = './'
    if plot_dir        is None: plot_dir        = './'
    if save_plot       is None: save_plot       = True
    if plot_format     is None: plot_format     = 'png'
    if output_name_tag is None: output_name_tag = ''
    if verbosity       is None: verbosity       = 1
    dateStr = date.today().strftime('%y%m%d')

    lightYieldRefPath = s1corr_paths[0]
    deadtime_type  = None
    generalBinN    = 400    
    countRangeSigN = 3
    chunkCheckTags = False
    combTimeThres  = 6*3600.0
    fix_crystalBall_shape = False
    crystalBall_alpha     = 1.0
    crystalBall_n         = 2.0
    s1_ampBinN       = 800
    s1_ampRange      = [1.0E4, 1.0E5]
    chi2DisplayRange = [0, 'inf']
    dateDisplayRange = ['-inf', 'inf']
    s1_corrFDisplay_forTheta        = [2.6E4, 4.0E4]
    s1_corrF_LYcorrDisplay_forTheta = [2.6E4, 4.0E4]
    peakDisplaySubset = ['Po210']
    if config_dict is not None:
        if 'deadtime_type'    in config_dict: deadtime_type  = config_dict['deadtime_type']
        if 'generalBinN'      in config_dict: generalBinN    = config_dict['generalBinN']
        if 'countRangeSigN'   in config_dict: countRangeSigN = config_dict['countRangeSigN']
        if 'chunkCheckTags'   in config_dict: chunkCheckTags = config_dict['chunkCheckTags']
        if 'combTimeThres'    in config_dict: combTimeThres  = config_dict['combTimeThres']
        if 'fix_crystalBall_shape' in config_dict: fix_crystalBall_shape = config_dict['fix_crystalBall_shape']
        if 'crystalBall_alpha'     in config_dict: crystalBall_alpha     = config_dict['crystalBall_alpha']
        if 'crystalBall_n'         in config_dict: crystalBall_n         = config_dict['crystalBall_n']
        if 's1_ampBinN'       in config_dict: s1_ampBinN       = config_dict['s1_ampBinN']
        if 's1_ampRange'      in config_dict: s1_ampRange      = config_dict['s1_ampRange']
        if 'chi2DisplayRange' in config_dict: chi2DisplayRange = config_dict['chi2DisplayRange'] 
        if 'dateDisplayRange' in config_dict: dateDisplayRange = config_dict['dateDisplayRange']
        if 's1_corrFDisplay_forTheta' in config_dict: s1_corrFDisplay_forTheta = config_dict['s1_corrFDisplay_forTheta']
        if 's1_corrF_LYcorrDisplay_forTheta' in config_dict: 
            s1_corrF_LYcorrDisplay_forTheta = config_dict['s1_corrF_LYcorrDisplay_forTheta']
        if 'peakDisplaySubset' in config_dict: peakDisplaySubset = config_dict['peakDisplaySubset']
    inputPars = {'s1corr_paths':         s1corr_paths,
                 'light_yield_ref_path': lightYieldRefPath,
                 'deadtime_type':     deadtime_type,
                 'generalBinN':       generalBinN,
                 'countRangeSigN':    countRangeSigN,
                 'chunkCheckTags':    chunkCheckTags,
                 'combTimeThres':     combTimeThres,
                 'crystalBall_alpha': crystalBall_alpha,
                 'crystalBall_n':     crystalBall_n,
                 's1_ampBinN':        s1_ampBinN,
                 's1_ampRange':       s1_ampRange,
                 'chi2DisplayRange':  chi2DisplayRange,
                 'dateDisplayRange':  dateDisplayRange,
                 's1_corrFDisplay_forTheta': s1_corrFDisplay_forTheta,
                 's1_corrF_LYcorrDisplay_forTheta': s1_corrF_LYcorrDisplay_forTheta,
                 'peakDisplaySubset': peakDisplaySubset}
    for chi2DispIdx in range(len(chi2DisplayRange)):
        if chi2DisplayRange[chi2DispIdx] =='-inf': chi2DisplayRange[chi2DispIdx] =-np.inf
        if chi2DisplayRange[chi2DispIdx] == 'inf': chi2DisplayRange[chi2DispIdx] = np.inf
    countRangeSet = config_dict['countRangeSet'] 

    if verbosity >= 1: print('  loading light yield reference file:', lightYieldRefPath)
    lightYieldRefPickle, cali_factor = None, None
    with open(lightYieldRefPath, 'rb') as handle:
        lightYieldRefPickle = pickle.load(handle)
    cali_factor = copy.deepcopy(lightYieldRefPickle['cali_factor'])
    if verbosity >= 1: print('    cali_factor =', cali_factor, '/MeV')
 
    if verbosity >= 1: print('  loading correction_derivation() files:', s1corr_paths) 
    corrPar, corrParPickle = {}, None
    fitRangeSetPickle, peak_normParPickle, peak_normErrPickle = {}, {}, {}
    for s1corrPath in s1corr_paths:
        with open(s1corrPath, 'rb') as handle:
            corrParPickle = pickle.load(handle)
        norm_factor = cali_factor/corrParPickle['cali_factor']
        for key in corrParPickle:
            if key not in corrPar:
                if isinstance(corrParPickle[key], dict) == False: 
                    corrPar[key] = np.array(corrParPickle[key])
                    if verbosity >= 1: print('   ', key, '=', corrPar[key])
        for peakName in corrParPickle['fitRangeSet']:
            if peakName not in fitRangeSetPickle:
                if (peakName not in ['Po218', 'Bi212']) or\
                   ((peakName == 'Po218') and ('Bi212' not in fitRangeSetPickle)) or\
                   ((peakName == 'Bi212') and ('Po218' not in fitRangeSetPickle)):
                    fitRangeSetPickle [peakName] = np.array(corrParPickle['fitRangeSet'][peakName]) *norm_factor
                    peak_normParPickle[peakName] = np.array(corrParPickle['peak_normPar'][peakName])*norm_factor
                    peak_normErrPickle[peakName] = np.array(corrParPickle['peak_normErr'][peakName])*norm_factor
                    if verbosity >= 1: 
                        print('    fitRangeSet['+peakName+']  =', fitRangeSetPickle[peakName])
                        print('    peak_normPar['+peakName+'] =', peak_normParPickle[peakName])
                        print('    peak_normErr['+peakName+'] =', peak_normErrPickle[peakName])
   
    if verbosity >= 1: print('  loading estimate_fiducial_mass() file:', fidMass_path)
    fidMassPar = None
    with open(fidMass_path, 'rb') as handle:
        fidMassPar = pickle.load(handle)
    if verbosity >= 1:
        for key in fidMassPar: print('   ', key, '=', fidMassPar[key])

    fitLightYieldCorrPeaks = ['Po210']
    fitLightYieldCorrCountThres = 20
    if config_dict is not None:      
        if 'fitLightYieldCorrPeaks' in config_dict: fitLightYieldCorrPeaks = config_dict['fitLightYieldCorrPeaks']
        if 'fitLightYieldCorrCountThres' in config_dict:      
            fitLightYieldCorrCountThres = config_dict['fitLightYieldCorrCountThres']
    #####################################################################################################################
    zName, r2name, r2Range = None, None, None
    if 'r2_naive'  in procData:
        zName   = 'area_fraction_top'
        r2name  = 'r2_naive'
        zRange  = [0.0, 1.0]
        r2Range = [0.0, 3.0E2]
    if 'r2_s1_cnn' in procData:
        zName  = 'z_s1_cnn'
        r2name = 'r2_s1_cnn'
        zRange  = [-150.0, 0.0]
        r2Range = [0.0, 6.0E3]
    r2_corrRange = [0.0, 6.0E3]
    if r2name is None:
        raise AssertionError('surface_rate_estimation(): either r2_naive or r2_s1_cnn are required for s1only analysis')
    #####################################################################################################################
    mask0 = np.copy(mask_dict['mask_base'])
    mask1 = np.copy(mask_dict['mask_final_surface'])
    #####################################################################################################################
    if verbosity >= 1: print('  base data:')
    plot_name = 's1only_surface_rate_estimation_0base' 
    if save_plot == True:
        plt, ax = plotStandard_s1only_highA(procData, mask0, mask1, verbosity=verbosity,binN=generalBinN,\
                                            r2Key=r2name, r2Range=r2Range, zKey=zName, zRange=zRange)
        if plot_extra is not None: plot_extra(plot_name, ax, mask_dict) 
        plot_path = os.path.join(plot_dir, dateStr+output_name_tag+plot_name+'.'+plot_format)
        plt.savefig(plot_path, format=plot_format)
        if verbosity >= 1: print('    saving:', plot_path)
        plt.close('all')
    ###S1 Corr###############################################################################
    plot_name = 's1only_surface_rate_estimation_1s1corr'
    if save_plot == True:
        plt, ax = plotStandard_s1only_highA(procData, mask0, mask1, verbosity=verbosity, binN=generalBinN,\
                                            ampKey='s1_corrF', r2Key='r2_corr', r2Range=r2_corrRange, zKey=zName,\
                                            zRange=zRange)
        if plot_extra is not None: plot_extra(plot_name, ax, mask_dict)
        plot_path = os.path.join(plot_dir, dateStr+output_name_tag+plot_name+'.'+plot_format)
        plt.savefig(plot_path, format=plot_format)
        if verbosity >= 1: print('    saving:', plot_path)
        plt.close('all') 
    #####################################################################################################################
    ###Wall Spatial Distribution###############################################################################
    theta_cmin = 1.0        # NOTE: must be smaller than the maximum bin
    s1_corrF_wallMask = copy.deepcopy(mask_dict['mask_wall_save_data'])
    s1_corrF_wallMask &= (s1_corrFDisplay_forTheta[0] <= procData['s1_corrF'])
    s1_corrF_wallMask &= (procData['s1_corrF'] < s1_corrFDisplay_forTheta[1])
    if 'mask_wall_save_data' in mask_dict:
        plot_name = 's1only_surface_rate_estimation_2s1corr'
        if save_plot == True:
            plt, ax = plotStandard_s1only_highA(procData, mask0, s1_corrF_wallMask,verbosity=verbosity,binN=generalBinN,\
                                                ampKey='s1_corrF', r2Key='r2_corr', r2Range=r2_corrRange,\
                                                zKey=zName, zRange=zRange)
            if plot_extra is not None: plot_extra(plot_name, ax, mask_dict)
            plot_path = os.path.join(plot_dir, dateStr+output_name_tag+plot_name+'.'+plot_format)
            plt.savefig(plot_path, format=plot_format)
            if verbosity >= 1: print('    saving:', plot_path)
            plt.close('all')
        plot_name = 's1only_surface_rate_estimation_2spatial'
        if save_plot == True:
            gridSpec = [1, 1]
            figSize, marginRatio = getSizeMargin(gridSpec, subplotSize=[18.0, 8.0], plotMargin=[0.07, 0.1, 0.07, 0.02])
            fig = matplotlib.pyplot.figure(figsize=figSize); fig.subplots_adjust(*marginRatio)
            gs = gridspec.GridSpec(*gridSpec)
            matplotlib.rc('xtick', labelsize=18)
            matplotlib.rc('ytick', labelsize=18)
            ax = list()
            for axIdx in range(gs.nrows*gs.ncols):
                ax.append(fig.add_subplot(gs[axIdx]));
                ax[-1].ticklabel_format(style='sci', scilimits=(-3, 3))
                ax[-1].xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
                ax[-1].yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())    
            ######
            axIdx = 0
            thetaRange = [-np.pi, np.pi]
            rangeXYZ = [thetaRange, zRange]
            hist = ax[axIdx].hist2d(procData[thetaName][s1_corrF_wallMask],procData[zName][s1_corrF_wallMask], 
                                    bins=[3*generalBinN, generalBinN], range=[rangeXYZ[0], rangeXYZ[1]],
                                    cmin=theta_cmin, norm=LogNorm(), cmap=CMAP)
            divider = make_axes_locatable(ax[axIdx])
            cax = divider.append_axes('right', size='5%', pad=0.1)
            cbar = fig.colorbar(hist[3], ax=ax[axIdx], cax=cax)
            ax[axIdx].get_yaxis().set_visible(False)
            ax[axIdx].set_xlabel(thetaName, fontsize=18)
            ax[axIdx].set_xlim(*rangeXYZ[0])
            ax[axIdx].set_ylim(*rangeXYZ[1])
            plot_pillars_theta(ax[axIdx], zorder=10)
            ###
            axUpper = divider.append_axes('top', size='60%', pad='6%')
            ax[axIdx].figure.add_axes(axUpper)
            hist = axUpper.hist(procData[thetaName][s1_corrF_wallMask], bins=3*generalBinN, range=rangeXYZ[0],\
                                color='steelblue', linewidth=1)
            axUpper.set_ylabel('counts', fontsize=18)
            axUpper.get_xaxis().set_visible(False)
            axUpper.set_xlim(*rangeXYZ[0])
            plot_pillars_theta(axUpper, zorder=10)
            ###
            axSide = divider.append_axes('left', size='15%', pad='2%')
            ax[axIdx].figure.add_axes(axSide)
            axSide.set_xlabel('counts', fontsize=18)
            axSide.set_ylabel(zName, fontsize=18)
            axSide.ticklabel_format(style='sci', scilimits=(-2, 2))
            axSide.xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
            axSide.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
            hist = axSide.hist(procData[zName][s1_corrF_wallMask], bins=generalBinN,range=rangeXYZ[1],\
                               color='steelblue', linewidth=1, orientation='horizontal')
            axSide.set_ylim(*rangeXYZ[1])
            xRange = axSide.get_xlim()
            axSide.set_xlim(xRange[1], xRange[0])
           
            if plot_extra is not None: plot_extra(plot_name, ax, mask_dict)
            plot_path = os.path.join(plot_dir, dateStr+output_name_tag+plot_name+'.'+plot_format)
            matplotlib.pyplot.savefig(plot_path, format=plot_format)
            if verbosity >= 1: print('    saving:', plot_path)
            matplotlib.pyplot.close('all')
    #####
    # Chunk Multi-Peak Rate with No Light Yeild Correction
    #####################################################################################################################
    plot_name = 's1only_surface_rate_estimation_3peak_fit0'
    if verbosity >= 1: print('  multi-peak fit, no LY corr:')
    if len(procData['run_number'][mask0]) == 0: 
        raise AssertionError('surface_rate_estimation(): mask0 cuts out all events')
    time_init = toBeginningOfDay(1.0E-9*procData['time'][mask0][0])
    time_fin  = 1.0E-9*procData['endtime'][mask0][-1]

    peak_amps, peak_width_percent = {}, {}
    for peakIdx, peakName in enumerate(fitRangeSetPickle):
        peak_amps[peakName]          = peak_normParPickle[peakName][0]
        peak_width_percent[peakName] = 100.0*peak_normParPickle[peakName][1]/peak_normParPickle[peakName][0]

    if verbosity >= 1:
        for peakName in countRangeSet: print('    countRangeSet['+peakName+'] =',  countRangeSet[peakName])
    ###########################################################
    runChunks, runModes = genRunChunks(np.unique(procData['run_number']), runInfoData, combTimeThres,\
                                       checkTags=chunkCheckTags, verbosity=(verbosity-1))
    timestampRanges, liveTimes, dayArr = [], [], []
    eventDataS1ampArr, eventData_chunk = [], []
    eventCounts = {peakName:[] for peakName in countRangeSet}
    eventRates  = {peakName:[] for peakName in countRangeSet}
    if (deadtime_type is not None) and ((deadtime_type in runInfoData) == False): deadtime_type == None
    for runChunk in tqdm(runChunks, disable=(verbosity<=0)):
        timestampRangeChunk, liveTimeChunk = [], ufloat(0, 0)
        mask_run_sel = (procData['run_number'] == -1)
        for run_number in runChunk:      
            timestampRangeChunk.append([runInfoData['start_timestamp'][runInfoData['name']==run_number][0],
                                        runInfoData['end_timestamp'][runInfoData['name']==run_number][0]])
            liveTime = pd.Timedelta(runInfoData['livetime'][runInfoData['name']==run_number][0]).total_seconds()
            if deadtime_type is not None: liveTime-=runInfoData[deadtime_type][runInfoData['name']==run_number][0]*1.0E-9
            else: warnings.warn('surface_rate_estimation(): deadtime plugin NOT implemented', Warning)
            liveTimeChunk += ufloat(liveTime, 0.5)          #NOTE: 0.5 gps uncertainty for livetime
            mask_run_sel |= (procData['run_number'] == run_number)
        timestampRanges.append(timestampRangeChunk)
        liveTimes      .append(liveTimeChunk)
        dayArr         .append((timestampRangeChunk[0][0] - time_init)/(3600.0*24.0))

        eventDataS1amp = np.array(procData['s1_corrF'][mask1 & mask_run_sel], dtype='float64')
        eventDataS1ampArr.append(eventDataS1amp)
        eventData_chunk.append({key:procData[key][mask_dict['mask_wall_save_data'] & mask_run_sel]\
                                for key in procData})
        countDataSet = {}
        for peakName in countRangeSet:
            peakData = eventDataS1amp[(countRangeSet[peakName][0] <= eventDataS1amp) &
                                      (eventDataS1amp < countRangeSet[peakName][1])]
            countDataSet[peakName] = len(peakData)
            eventCounts[peakName].append(len(peakData))
        activityPars = [countDataSet, countRangeSet, {}, liveTimeChunk,\
                        2*np.pi*mask_dict['tpc_r_cm']*fidMassPar['XeFidZ_wall']]
        activity  = get_activity_surface_count(peakName, *activityPars)
        #NOTE: s1only surface rate correction from s1s2
        activity *= fidMassPar['validS1S2ZRatio']
        for peakName in countRangeSet: eventRates[peakName].append(activity)
    ###########################################################
    if save_plot == True:
        gridSpec = [2, 1]
        figSize, marginRatio = getSizeMargin(gridSpec, subplotSize=[18.0, 5.0], plotMargin=[0.08, 0.14, 0.04, 0.06])
        fig = plt.figure(figsize=figSize); fig.subplots_adjust(*marginRatio)
        gs = gridspec.GridSpec(*gridSpec)
        matplotlib.rc('xtick', labelsize=18)
        matplotlib.rc('ytick', labelsize=18)
        ax = []
        for axIdx in range(gs.nrows*gs.ncols):
            ax.append(fig.add_subplot(gs[axIdx]));
            ax[-1].ticklabel_format(style='sci', scilimits=(-2, 2))
            ax[-1].xaxis.set_tick_params(which='major', width=2.0, length=4.0)
            ax[-1].xaxis.set_tick_params(which='minor', width=1.0, length=2.0)
            ax[-1].xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
            ax[-1].yaxis.set_tick_params(which='major', width=2.0, length=4.0)
            ax[-1].yaxis.set_tick_params(which='minor', width=1.0, length=2.0)
            ax[-1].yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
        ######
        axIdx = 0
        hist = ax[axIdx].hist(procData['s1_corrF'][mask1], bins=s1_ampBinN, range=s1_ampRange, color='steelblue',\
                              linewidth=1)
        ax[axIdx].set_xlabel('s1_corrF', fontsize=18)
        ax[axIdx].set_ylabel('counts', fontsize=18)
        ax[axIdx].set_xlim(s1_ampRange)
        ax[axIdx].set_ylim(bottom=0.5)
        ax[axIdx].set_yscale('log')

        for peakName in countRangeSet:
            mask_isotope  = (countRangeSet[peakName][0] <= procData['s1_corrF']) 
            mask_isotope &= (procData['s1_corrF'] < countRangeSet[peakName][1])
            mask_isotope &= mask1
            hist = ax[axIdx].hist(procData['s1_corrF'][mask_isotope], bins=s1_ampBinN, range=s1_ampRange,
                                  color='gainsboro', linewidth=1, alpha=0.3)
        ######
        axIdx = 1
        tracePlot, tracePlotText = {}, {}
        for peakName in countRangeSet:
            tracePlot[peakName] = ax[axIdx].errorbar(dayArr, [item.n for item in eventRates[peakName]], 
                                                              [item.s for item in eventRates[peakName]], 
                                                      color=ALPHA_LINE_C[peakName], fmt='s', **ERROR_PLOT_PARS)[0]
            tracePlotText[peakName] = peakName
        setXaxisDateLabel(ax[axIdx], time_init, [time_init, time_fin])
        ax[axIdx].set_ylabel('activity concentration (mBq/m^2)', fontsize=18)
        ax[axIdx].set_yscale('log')
        ax[axIdx].secondary_yaxis('right')
        if plot_extra is not None: plot_extra(plot_name, ax, mask_dict) 
        plot_path = os.path.join(plot_dir, dateStr+output_name_tag+plot_name+'.'+plot_format)
        plt.savefig(plot_path, format=plot_format)
        if verbosity >= 1: print('    saving:', plot_path)
        plt.close('all')
    #####
    # Chunk Multi-peak Fit, Light Yeild Correction
    #####################################################################################################################
    plot_dir_name   = 's1only_surface_rate_estimation_4chunk_fit'
    plot_name_chunk = 's1only_surface_rate_estimation_4chunk_multipeak_fit'

    peakFitFunc = multiSigLinkedGausPDF_withOneCrystalBall_scipy
    output_path = os.path.join(output_dir, output_name_tag+'s1only_surface_rate_estimation.pkl')
    ###########################################################
    fitLightYieldCorrPeaks = [fitLightYieldCorrPeak for fitLightYieldCorrPeak in fitLightYieldCorrPeaks
                              if fitLightYieldCorrPeak in countRangeSet]
    fitRanges_gainCorr, peakFitPars, relative_light_yields, corr_factors, chi2_ndofs = {}, [], [], [], []
    eventData = {}
    ###########################################################
    RnRatePickle, inputDictToExtend = None, None    
    inputDictToExtend = {'runChunks': runChunks, 'runModes': runModes, 'timestampRanges': timestampRanges,\
                         'liveTimes': liveTimes, 'eventCounts': eventCounts, 'eventRates': eventRates}
    RnRatePickle = loadPickleRateData(output_path, dayArr, time_init, inputDictToExtend=inputDictToExtend)
    if RnRatePickle is not None:
        if fitLightYieldCorrPeaks != RnRatePickle['fitLightYieldCorrPeaks']:
            raise AssertionError('surface_rate_estimation(): fitLightYieldCorrPeaks =',fitLightYieldCorrPeaks,'\nwhile',\
                                 output_path, 'provides fitLightYieldCorrPeaks =',RnRatePickle['fitLightYieldCorrPeaks'])
        for peakName in fitLightYieldCorrPeaks:
            fitRanges_gainCorr[peakName] = copy.deepcopy(RnRatePickle['fitRanges_gainCorr'][peakName])
        peakFitPars           = copy.deepcopy(RnRatePickle['peakFitPars'])
        relative_light_yields = copy.deepcopy(RnRatePickle['relative_light_yields'])
        corr_factors          = copy.deepcopy(RnRatePickle['corr_factors'])
        chi2_ndofs            = copy.deepcopy(RnRatePickle['chi2_ndofs'])
        if 'multipeak_chi2_ndof' in RnRatePickle['eventData']:
            eventData = copy.deepcopy(RnRatePickle['eventData'])
            timeMask = (eventData['run_number'] < runChunks[0][0])
            eventData = {key:eventData[key][timeMask] for key in eventData}
    ###########################################################
    if verbosity >= 1: 
        print('  fitLightYieldCorrPeaks:', fitLightYieldCorrPeaks)
        print('  chunk fits to process:', len(eventDataS1ampArr), '\n')
    if len(eventDataS1ampArr) == 0: raise AssertionError('surface_rate_estimation(): no run chunks to process')
    fitLightYieldCorrPeaks_noPo210 = copy.deepcopy(fitLightYieldCorrPeaks)
    fitLightYieldCorrPeaks_noPo210.remove('Po210')
    fitX = np.linspace(*s1_ampRange, 10*s1_ampBinN+1)[:-1]
    s1ampXbinSize, fitY, fitYInit = None, None, None
    for chunkIdx, eventDataS1amp in enumerate(eventDataS1ampArr):
        eventDataS1ampC = np.array(eventDataS1amp, dtype='float64')
        s1ampHist, s1ampX = np.histogram(eventDataS1ampC, bins=s1_ampBinN, range=s1_ampRange)
        s1ampXcenter = 0.5*(s1ampX[1:] + s1ampX[:-1])
        s1ampX = s1ampX[:-1]
        s1ampXbinSize = s1ampX[1] - s1ampX[0]
        s1ampErr = np.array([np.sqrt(count) if count != 0 else 1.0 for count in s1ampHist])
  
        ###########################scaling for the best fit range 
        scaleStep   = 0.02
        segIdxRange = [-25, 25] 
        if config_dict is not None:
            if 'scaleStep'   in config_dict: scaleStep   = config_dict['scaleStep']
            if 'segIdxRange' in config_dict: segIdxRange = config_dict['segIdxRange']
        fitRangeChunk = {}
        optCount, optIdx, optIdx_sameN = -1, 0, 0
        if verbosity >= 3: print('  find peak fit range shift:')
        for segIdx in range(*segIdxRange, 1):
            countInRange = 0
            for peakName in fitLightYieldCorrPeaks:
                fitRangeChunk[peakName] = np.array(countRangeSet[peakName])*(1.0 + segIdx*scaleStep)
                countInRange += len(eventDataS1ampC[(fitRangeChunk[peakName][0] <= eventDataS1ampC) &
                                                    (eventDataS1ampC < fitRangeChunk[peakName][1])])
            if countInRange >= optCount:
                if countInRange == optCount: optIdx_sameN += 1
                else:                        optIdx_sameN  = 0
                optCount = np.copy(countInRange)
                optIdx   = np.copy(segIdx)
            if verbosity >= 3: print('   ', segIdx, segIdx*scaleStep, countInRange, optCount)
        optIdx -= int(optIdx_sameN/2)
        for peakName in fitLightYieldCorrPeaks:
            fitRangeChunk[peakName] = np.array(countRangeSet[peakName])*(1.0 + optIdx*scaleStep)
            if peakName not in fitRanges_gainCorr: fitRanges_gainCorr[peakName] = []
            fitRanges_gainCorr[peakName].append(fitRangeChunk[peakName])
        ##################################################################
        if verbosity >= 1: 
            print(chunkIdx, '/', len(eventDataS1ampArr), runChunks[chunkIdx])
            print('  day:     ', dayArr[chunkIdx])
            print('  eventN:  ', len(eventDataS1ampC))
            print('  livetime:', liveTimes[chunkIdx])
        countInRange = 0
        for peakName in fitLightYieldCorrPeaks:
            if verbosity >= 1: print('   ', peakName, countRangeSet[peakName], fitRangeChunk[peakName])
            countInRange += len(eventDataS1ampC[(fitRangeChunk[peakName][0] <= eventDataS1ampC) &
                                                (eventDataS1ampC < fitRangeChunk[peakName][1])])
        if verbosity >= 2: print('    optIdx =', optIdx, ', countInRange =', countInRange)
        ##################################################################
        fitPars, fitErrs, fitCov = None, None, None
        if countInRange >= fitLightYieldCorrCountThres:
            fitDataSet, initPars, initPars_crystalBall = {}, [], []
            fitRangeMask = np.full(len(s1ampX), False)
            for peakName in fitLightYieldCorrPeaks:
                fitDataSet[peakName] = eventDataS1ampC[(fitRangeChunk[peakName][0] <= eventDataS1ampC)&
                                                       (eventDataS1ampC < fitRangeChunk[peakName][1])] 
                fitRangeMask |= ((fitRangeChunk[peakName][0] < s1ampX) & (s1ampX < fitRangeChunk[peakName][1]))
                if peakName == 'Po210':
                    initPars_crystalBall.append(len(fitDataSet[peakName])*(s1ampX[1]-s1ampX[0]))            #A
                    initPars_crystalBall.append(crystalBall_alpha)                                          #alpha / beta
                    initPars_crystalBall.append(crystalBall_n)                                              #n     / m
                    if len(fitDataSet[peakName]) != 0: initPars_crystalBall.append(np.median(fitDataSet[peakName]))
                    else:                              initPars_crystalBall.append(np.median(fitRangeChunk[peakName]))
                    initPars_crystalBall.append(np.std(fitDataSet['Po210'])/2.0)
                else:
                    initPars.append(len(fitDataSet[peakName])*(s1ampX[1]-s1ampX[0]))
                    if len(fitDataSet[peakName]) != 0: initPars.append(np.mean(fitDataSet[peakName]))
                    else:                              initPars.append(np.mean(fitRangeChunk[peakName]))
            if len(initPars) != 0: initPars.append(np.std(fitDataSet[fitLightYieldCorrPeaks[-1]]))
            initPars.extend(initPars_crystalBall)

            costFunc = iminuit.cost.LeastSquares(s1ampX, s1ampHist, s1ampErr, peakFitFunc)
            costFunc.mask = fitRangeMask
            objMinuit = iminuit.Minuit(costFunc, initPars)
            for peakIdx, peakName in enumerate(fitLightYieldCorrPeaks_noPo210):
                peakLocRange = np.array(fitRangeChunk[peakName])
                peakLocRangeLen = peakLocRange[1] - peakLocRange[0]
                peakLocRange[0] += peakLocRangeLen/4.0
                peakLocRange[1] -= peakLocRangeLen/4.0
                objMinuit.limits[2*peakIdx]   = [0.1*initPars[2*peakIdx], 10.0*initPars[2*peakIdx]]
                objMinuit.limits[2*peakIdx+1] = peakLocRange
                if initPars[2*peakIdx] == 0:
                    objMinuit.fixed[2*peakIdx]   = True
                    objMinuit.fixed[2*peakIdx+1] = True
            if len(initPars_crystalBall) == 0:
                objMinuit.limits[-1] = [0.8*initPars[-1], 1.2*initPars[-1]]     #NOTE the fit limit for the peak width
            else:
                if len(initPars) != len(initPars_crystalBall): objMinuit.limits[-6] = [0.8*initPars[-6],1.2*initPars[-6]]
                peakLocRange = np.array(fitRangeChunk['Po210'])
                peakLocRangeLen = peakLocRange[1] - peakLocRange[0]
                peakLocRange[0] += peakLocRangeLen/2.0
                objMinuit.limits[-5] = [0.1*initPars[-5], 10.0*initPars[-5]]
                objMinuit.limits[-2] = fitRangeChunk['Po210']  
                objMinuit.limits[-1] = [0.8*initPars[-1], 1.2*initPars[-1]]
                if fix_crystalBall_shape == True:
                    objMinuit.fixed[-4] = True
                    objMinuit.fixed[-5] = True
  
            if verbosity >= 2: print('    initial values:', initPars, '\n    fit limits:', objMinuit.limits)
            optResult = objMinuit.migrad()
            try:
                fitPars = [val for val in optResult.values]
                fitErrs = [err for err in optResult.errors]
                fitCov  = [row for row in optResult.covariance]
                fitPars_gaus = [val for val in fitPars[:-5]]
                fitErrs_gaus = [val for val in fitErrs[:-5]]
                fitPars_crystalBall = [val for val in fitPars[-5:]]
                fitErrs_crystalBall = [val for val in fitErrs[-5:]]
            except:
                pass
            if verbosity >= 1: print('    multi-peak fit:\n', optResult)
            ############################
            fitYinit = peakFitFunc(fitX, initPars)    
            if fitCov is not None:
                peak_normPar_chunk, peak_locs = {}, []
                for peakIdx, peakName in enumerate(fitLightYieldCorrPeaks_noPo210):
                    peak_normA    = ufloat(fitPars_gaus[2*peakIdx],   fitErrs_gaus[2*peakIdx])
                    peak_normMu   = ufloat(fitPars_gaus[2*peakIdx+1], fitErrs_gaus[2*peakIdx+1])
                    peak_normSig  = ufloat(fitPars_gaus[-1],          fitErrs_gaus[-1])
                    peak_normSig *= peak_normMu/ufloat(fitPars_gaus[1], fitErrs_gaus[1])
                    peak_normPar_chunk[peakName] = [peak_normA, peak_normMu, peak_normSig]
                    if peakName != 'Rn222': peak_locs.append(peak_normMu/ALPHA_LINE_E[peakName])
                if 'Po210' in fitLightYieldCorrPeaks:
                    peak_normA   = ufloat(fitPars_crystalBall[0], fitErrs_crystalBall[0])
                    peak_normMu  = ufloat(fitPars_crystalBall[3], fitErrs_crystalBall[3])
                    peak_normSig = ufloat(fitPars_crystalBall[4], fitErrs_crystalBall[4])
                    peak_norm_a  = ufloat(fitPars_crystalBall[1], fitErrs_crystalBall[1])
                    peak_norm_n  = ufloat(fitPars_crystalBall[2], fitErrs_crystalBall[2])
                    peak_normPar_chunk['Po210'] = [peak_normA, peak_normMu, peak_normSig, peak_norm_a, peak_norm_n]
                    #peak_locs.append(peak_normMu/ALPHA_LINE_E['Po210'])
                peak_loc_ave = weightedAverage(peak_locs)
                peakFitPars          .append(peak_normPar_chunk)
                relative_light_yields.append(peak_loc_ave/cali_factor)
                corr_factors         .append(peak_loc_ave)
                chi2_ndofs           .append(optResult.fval/optResult.ndof)
                if verbosity >= 1: 
                    print('    relative_light_yield =', relative_light_yields[-1])
                    print('    chi2_ndof =', chi2_ndofs[-1])
                fitY= peakFitFunc(fitX, fitPars)
                fitResidual = np.array([np.sqrt(iminuit.cost.poisson_chi2(np.array([1.0*histYval]),\
                                                                          np.array([1.0*fitYval])))
                                        for histYval, fitYval in zip(s1ampHist, peakFitFunc(s1ampX, fitPars))])
            else:
                peak_normPar_chunk = {}
                for peakIdx, peakName in enumerate(fitLightYieldCorrPeaks): peak_normPar_chunk[peakName] = None
                peakFitPars          .append(peak_normPar_chunk)
                relative_light_yields.append(ufloat(-1, 0))
                corr_factors         .append(ufloat(-1, 0))
                chi2_ndofs           .append(-1)
                fitY = None
        else:
            peak_normPar_chunk = {}
            for peakIdx, peakName in enumerate(fitLightYieldCorrPeaks): peak_normPar_chunk[peakName] = None
            peakFitPars          .append(peak_normPar_chunk)
            relative_light_yields.append(ufloat(-1, 0))
            corr_factors         .append(ufloat(-1, 0))
            chi2_ndofs           .append(-1)
            fitYinit = None
            fitY     = None
        ###########################################################
        eventData_chunk[chunkIdx]['relative_light_yield'] = 0*eventData_chunk[chunkIdx]['s1_corrF'] +\
                                                            relative_light_yields[-1].n
        eventData_chunk[chunkIdx]['multipeak_chi2_ndof']  = 0*eventData_chunk[chunkIdx]['s1_corrF'] + chi2_ndofs[-1]
        eventData_chunk[chunkIdx]['s1_corrF_LYcorr']      = eventData_chunk[chunkIdx]['s1_corrF']/\
                                                            relative_light_yields[-1].n
        if len(eventData) == 0: eventData = copy.deepcopy(eventData_chunk[chunkIdx])
        else:
            for key in eventData: 
                if key in eventData_chunk[chunkIdx]:
                    eventData[key] = np.append(eventData[key], eventData_chunk[chunkIdx][key]) 
                else:
                    eventData[key] =np.append(eventData[key],np.nan*np.empty(len(eventData_chunk[chunkIdx]['s1_corrF'])))
        ###########################################################
        if save_plot == True:
            import matplotlib.pyplot as plt
            gridSpec = [1, 1]
            figSize, marginRatio = getSizeMargin(gridSpec, subplotSize=[18.0, 5.0], plotMargin=[0.06, 0.16, 0.02, 0.12])
            fig = plt.figure(figsize=figSize); fig.subplots_adjust(*marginRatio)
            gs = gridspec.GridSpec(*gridSpec)
            matplotlib.rc('xtick', labelsize=18)
            matplotlib.rc('ytick', labelsize=18)
            ax = []
            for axIdx in range(gs.nrows*gs.ncols):
                ax.append(fig.add_subplot(gs[axIdx]));
                ax[-1].ticklabel_format(style='sci', scilimits=(-2, 2))
                ax[-1].xaxis.set_tick_params(width=1.5, length=4)
                ax[-1].xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
                ax[-1].yaxis.set_tick_params(width=1.5, length=4)
                ax[-1].yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
            ######
            axIdx = 0
            hist = ax[axIdx].hist(eventDataS1ampC, bins=s1_ampBinN, range=s1_ampRange, color='steelblue', linewidth=1)
            runChunkTitle = [str(runChunk[runIdx:runIdx          + 12])
                             for runIdx in range(0, len(runChunk), 12)]
            dateInit = timestampToDate(timestampRanges[chunkIdx][0][0]).split('\n')
            dateInit = dateInit[0] + '_' + dateInit[1]
            dateFin  = timestampToDate(timestampRanges[chunkIdx][-1][1]).split('\n')
            dateFin  = dateFin[0]  + '_' + dateFin[1]
            titleStr = ''
            titleStr += dateInit + ' to ' + dateFin + '\n'
            # titleStr += 'Run: ' + '\n'.join(runChunkTitle) + '\n'
            titleStr += 'Run: ' + str(runChunks[chunkIdx][0]) + '-' + str(runChunks[chunkIdx][-1]) 
            titleStr += ' ' + runModes[chunkIdx][0] + ', '
            titleStr += 'Count: ' + str(len(eventDataS1ampC)) + ', '
            titleStr += 'Duration: ' + str(liveTimes[chunkIdx]) + ' s'
            ax[axIdx].set_title(titleStr, fontsize=16, y=1.03)
            ax[axIdx].set_xlabel('s1_corrF', fontsize=18)
            ax[axIdx].set_ylabel('counts', fontsize=18)
            ax[axIdx].set_xlim(s1_ampRange) 
            ax[axIdx].set_ylim(bottom=0.5)
            ax[axIdx].set_yscale('log')
            if fitY is not None:
                ax[axIdx].get_xaxis().set_visible(False)
                divider = make_axes_locatable(ax[axIdx])
                axRes = divider.append_axes('bottom', size='20%', pad='2%')
                ax[axIdx].figure.add_axes(axRes)
                axRes.set_xlim(s1_ampRange)
                axRes.set_xlabel('s1_corrF', fontsize=18)
                axRes.set_ylabel('fit\nresidual', fontsize=18)
                axRes.ticklabel_format(style='sci', scilimits=(-2, 2))
                axRes.xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
                axRes.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())    
                
                linePlot, linePlotText= {}, {}
                for peakName in fitLightYieldCorrPeaks:
                    fitMask = ((fitRangeChunk[peakName][0] < fitX) & (fitX < fitRangeChunk[peakName][1]))
                    ax[axIdx].plot(fitX[fitMask], fitY[fitMask], color='red', linewidth=1.5, alpha=0.8, zorder=3)
                    linePlot[peakName] = ax[axIdx].axvline(x=peak_normPar_chunk[peakName][1].n, 
                                                           color=ALPHA_LINE_C[peakName],linestyle='dashed',linewidth=1.5)
                    linePlotText[peakName] = peakName + ''
                for peakName in ALPHA_LINE_E:
                    if peakName not in fitLightYieldCorrPeaks:
                        linePlot[peakName] = ax[axIdx].axvline(x=corr_factors[-1].n*ALPHA_LINE_E[peakName],
                                                               color=ALPHA_LINE_C[peakName], linestyle='dashed', 
                                                               linewidth=1.5)
                        linePlotText[peakName] = peakName + ''        
                legObj = ax[axIdx].legend([linePlot[peakName] for peakName in ALPHA_LINE_E],\
                                          [linePlotText[peakName] for peakName in ALPHA_LINE_E],\
                                          title=('chi2/ndof=' + scientificStr(chi2_ndofs[-1])),\
                                          loc='upper left', fontsize=16)
                for peakName in fitLightYieldCorrPeaks:
                    fitMask = ((fitRangeChunk[peakName][0] <= s1ampX) & (s1ampX < fitRangeChunk[peakName][1]))
                    axRes.plot(s1ampX[fitMask], fitResidual[fitMask], color='black')
                axRes.axhline(y=1.0, color='grey', linestyle='dashed', zorder=0)
                ylimRes = list(axRes.get_ylim())
                ylimRes[0] = 0.0
                if ylimRes[1] <= 5.0: ylimRes[1] = 5.0
                else:                 ylimRes[1] = roundSig(ylimRes[1], 2)
                axRes.set_ylim(*ylimRes)
                plt.setp(legObj.get_title(), fontsize=16) 
            elif fitYinit is not None:
                fitMask = ((fitRangeChunk[peakName][0] < fitX) & (fitX < fitRangeChunk[peakName][1]))
                ax[axIdx].plot(fitX[fitMask], fitYinit[fitMask],color='purple',linewidth=2.0,linestyle='dashed',zorder=2)
            if plot_extra is not None: plot_extra(plot_name_chunk, ax, mask_dict) 
            plot_subdir = dateStr + output_name_tag + plot_dir_name
            pathlib.Path(os.path.join(plot_dir, plot_subdir)).mkdir(parents=True, exist_ok=True)
            plot_path = os.path.join(plot_dir, plot_subdir, dateStr+plot_name_chunk+\
                                     str(runChunks[chunkIdx][0])+'.'+plot_format)
            plt.savefig(plot_path, format=plot_format)
            if verbosity >= 1: print('    saving:', plot_path)
            plt.close('all')
    #####
    # Rate Evalution
    #####################################################################################################################
    plot_name = 's1only_surface_estimation_4peak_rate'
    if verbosity >= 1: print('  rate evalution:')    
    ###########################################################
    eventRates_fitA = {peakName:[] for peakName in fitLightYieldCorrPeaks}
    mask_isotope = {}
    chunkIdxOffset = 0
    for peakName in fitLightYieldCorrPeaks: mask_isotope[peakName] = np.copy(mask1)
    ###########################################################
    if RnRatePickle is not None:
        for peakName in RnRatePickle['eventRates_fitA']:
            eventRates_fitA[peakName]  = copy.copy(RnRatePickle['eventRates_fitA'][peakName])
    for chunkIdx, eventDataS1amp in enumerate(eventDataS1ampArr):  
        if verbosity >= 1: print(chunkIdx, runChunks[chunkIdx])
        if chi2_ndofs[chunkIdxOffset+chunkIdx] != -1:
            activityPars = [peakFitPars[chunkIdxOffset+chunkIdx], s1ampXbinSize, liveTimeChunk,\
                            2*np.pi*mask_dict['tpc_r_cm']*fidMassPar['XeFidZ_wall']]
            for peakName in fitLightYieldCorrPeaks:
                activity  = get_activity_surface_fit(peakName, *activityPars)
                activity *= fidMassPar['validS1S2ZRatio']               #NOTE: s1only surface rate correction from s1s2
                eventRates_fitA[peakName].append(activity)
            if verbosity >= 1:
                print('    liveTime    =', liveTimes[chunkIdx],   's')
                print('    XeFidZ_wall =', fidMassPar['XeFidZ_wall'],  'cm')
                for peakName in eventRates_fitA:
                    print('    eventRates_fitA['+peakName+'] =',  eventRates_fitA[peakName][-1])
        else:
            for peakName in fitLightYieldCorrPeaks: eventRates_fitA[peakName].append(ufloat(-1.0, 0.0))
    if RnRatePickle is not None:
        time_init = copy.copy(RnRatePickle['time_init'])
        runChunks       = copy.copy(RnRatePickle['runChunks'])
        runModes        = copy.copy(RnRatePickle['runModes'])
        timestampRanges = copy.copy(RnRatePickle['timestampRanges'])
        liveTimes       = copy.copy(RnRatePickle['liveTimes'])
        dayArr          = copy.copy(RnRatePickle['dayArr'])
        eventCounts       = copy.copy(RnRatePickle['eventCounts'])
        eventRates        = copy.copy(RnRatePickle['eventRates'])
    ###########################################################
    timeDisplay = genTimeDisplay(dateDisplayRange, time_init, time_fin)
    eventDisplayMask  = (chi2DisplayRange[0] <= eventData['multipeak_chi2_ndof'])
    eventDisplayMask &= (eventData['multipeak_chi2_ndof'] < chi2DisplayRange[1])
    eventDisplayMask &= (timeDisplay[0] <= eventData['time']*1.0E-9)
    eventDisplayMask &= (eventData['endtime']*1.0E-9 < timeDisplay[1])
    if 'z_s1_cnn' in procData:
        eventDisplayMask &= (mask_dict['z_s1_cnn_cutRange'][0] <= eventData['z_s1_cnn'])
        eventDisplayMask &= (eventData['z_s1_cnn'] < mask_dict['z_s1_cnn_cutRange'][1])
    else:
        eventDisplayMask &= (mask_dict['S1AFT_cutRange'][0] <= eventData['area_fraction_top'])
        eventDisplayMask &= (eventData['area_fraction_top'] < mask_dict['S1AFT_cutRange'][1])
    if verbosity >= 1:  
        print('    [time_init, time_fin] =', [time_init, time_fin], ', timeDisplay =', timeDisplay)
        for peakName in mask_isotope: 
            print('    mask_isotope[\''+peakName+'\'],', np.sum(mask_isotope[peakName]))
    ###########################################################
    if save_plot == True:
        gridSpec = [1, 1]
        figSize, marginRatio = getSizeMargin(gridSpec, subplotSize=[18.0, 5.0], plotMargin=[0.06, 0.14, 0.02, 0.08])
        fig = plt.figure(figsize=figSize); fig.subplots_adjust(*marginRatio)
        gs = gridspec.GridSpec(*gridSpec)
        matplotlib.rc('xtick', labelsize=18)
        matplotlib.rc('ytick', labelsize=18)
        ax = []
        for axIdx in range(gs.nrows*gs.ncols):
            ax.append(fig.add_subplot(gs[axIdx]));
            ax[-1].ticklabel_format(style='sci', scilimits=(-2, 2))
            ax[-1].xaxis.set_tick_params(width=1.5, length=4)
            ax[-1].xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
            ax[-1].yaxis.set_tick_params(width=1.5, length=4)
            ax[-1].yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
        ######
        axIdx = 0
        hist = ax[axIdx].hist(eventData['s1_corrF_LYcorr'][eventDisplayMask], bins=s1_ampBinN, range=s1_ampRange,\
                              color='steelblue', linewidth=1)
        ax[axIdx].set_xlabel('s1_corr3', fontsize=18)
        ax[axIdx].set_ylabel('counts', fontsize=18)
        ax[axIdx].set_xlim(s1_ampRange)
        ax[axIdx].set_ylim(bottom=0.5)
        ax[axIdx].set_yscale('log')

        linePlot, linePlotText= {}, {}
        for peakName in ALPHA_LINE_E:
            linePlot[peakName] = ax[axIdx].axvline(x=cali_factor*ALPHA_LINE_E[peakName], color=ALPHA_LINE_C[peakName], 
                                                   linestyle='dashed', linewidth=1.5)
            linePlotText[peakName] = peakName + ''
        legObj = ax[axIdx].legend([linePlot[peakName]     for peakName in ALPHA_LINE_E],\
                                  [linePlotText[peakName] for peakName in ALPHA_LINE_E], loc='upper left', fontsize=16)
        if plot_extra is not None: plot_extra(plot_name, ax, mask_dict) 
        plot_path = os.path.join(plot_dir, dateStr+output_name_tag+plot_name+'_.'+plot_format)
        plt.savefig(plot_path, format=plot_format)
        if verbosity >= 1: print('    saving:', plot_path)
        plt.close('all')
        ###########################################################
        gridSpec = [4, 1]
        figSize, marginRatio = getSizeMargin(gridSpec, subplotSize=[18.0, 4.0], plotMargin=[0.08, 0.16, 0.06, 0.06])
        fig = plt.figure(figsize=figSize); fig.subplots_adjust(*marginRatio)
        gs = gridspec.GridSpec(*gridSpec)
        matplotlib.rc('xtick', labelsize=18)
        matplotlib.rc('ytick', labelsize=18)
        ax = []
        for axIdx in range(gs.nrows*gs.ncols):
            ax.append(fig.add_subplot(gs[axIdx]));
            ax[-1].ticklabel_format(style='sci', scilimits=(-2, 2))
            ax[-1].xaxis.set_tick_params(which='major', width=2.0, length=4.0)
            ax[-1].xaxis.set_tick_params(which='minor', width=1.0, length=2.0)
            ax[-1].xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
            ax[-1].yaxis.set_tick_params(which='major', width=2.0, length=4.0)
            ax[-1].yaxis.set_tick_params(which='minor', width=1.0, length=2.0)
            ax[-1].yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
        ######
        validMask = (chi2DisplayRange[0] <= np.array(chi2_ndofs)) & (np.array(chi2_ndofs) < chi2DisplayRange[1])
        ######
        axIdx = 0
        tracePlot, tracePlotText = {}, {}
        for peakName in eventRates_fitA:
            tracePlot[peakName] = \
              ax[axIdx].errorbar(np.array(dayArr)[validMask],\
                                 [item.n for item in np.array(eventRates_fitA[peakName])[validMask]],\
                                 [item.s for item in np.array(eventRates_fitA[peakName])[validMask]],\
                                 color=ALPHA_LINE_C[peakName], fmt='s', **ERROR_PLOT_PARS)[0]
            tracePlotText[peakName] = peakName
        setXaxisDateLabel(ax[axIdx], time_init, timeDisplay, displayOff=True)
        ax[axIdx].set_ylabel('activity concentration\n(uBq/cm^2)', fontsize=18)
        ax[axIdx].set_yscale('log')
        ax[axIdx].secondary_yaxis('right')
        if 'Po218' in tracePlotText: tracePlotText['Po218'] = tracePlotText['Po218'].replace('Po218', 'Po218/Bi212')
        legObj = ax[axIdx].legend([tracePlot[peakName]     for peakName in eventRates_fitA\
                                   if (peakName in ALPHA_LINE_E) and (peakName not in peakDisplaySubset)],\
                                  [tracePlotText[peakName] for peakName in eventRates_fitA\
                                   if (peakName in ALPHA_LINE_E) and (peakName not in peakDisplaySubset)],\
                                  loc='upper left', fontsize=16, markerscale=1.5)
        ######
        axIdx = 1
        subTracePlot, subTracePlotText = {}, {}
        for peakName in peakDisplaySubset:
            subTracePlot[peakName] = \
              ax[axIdx].errorbar(np.array(dayArr)[validMask],\
                                 [item.n for item in np.array(eventRates_fitA[peakName])[validMask]],
                                 [item.s for item in np.array(eventRates_fitA[peakName])[validMask]],
                                 color=ALPHA_LINE_C[peakName], fmt='s', **ERROR_PLOT_PARS)[0]
            subTracePlotText[peakName] = peakName
        setXaxisDateLabel(ax[axIdx], time_init, timeDisplay, displayOff=True)
        ax[axIdx].set_ylabel('activity concentration\n(Bq/m^2)', fontsize=18)
        ax[axIdx].set_ylim(bottom=1.0)
        ax[axIdx].set_yscale('log')
        ax[axIdx].secondary_yaxis('right')
        legObj = ax[axIdx].legend([subTracePlot[peakName]     for peakName in peakDisplaySubset\
                                   if peakName in ALPHA_LINE_E],\
                                  [subTracePlotText[peakName] for peakName in peakDisplaySubset\
                                   if peakName in ALPHA_LINE_E],\
                                  loc='lower right', fontsize=14, markerscale=1.5)
        ######
        axIdx = 2
        ax[axIdx].errorbar(np.array(dayArr)[validMask], [item.n for item in np.array(relative_light_yields)[validMask]], 
                                                        [item.s for item in np.array(relative_light_yields)[validMask]], 
                           color='blue', fmt='s', **ERROR_PLOT_PARS)[0]
        setXaxisDateLabel(ax[axIdx], time_init, timeDisplay, displayOff=True)
        ax[axIdx].set_ylabel('relative light yield', fontsize=18)
        ax[axIdx].secondary_yaxis('right')
        ######
        axIdx = 3
        ax[axIdx].scatter(np.array(dayArr)[validMask], np.array(chi2_ndofs)[validMask], color='blue')
        setXaxisDateLabel(ax[axIdx], time_init, timeDisplay)
        ax[axIdx].set_ylabel('multi-peak fit chi2/ndof', fontsize=18)
        ax[axIdx].set_ylim(bottom=0.0)
        ax[axIdx].secondary_yaxis('right')
        ax[axIdx].axhline(y=1.0, color='grey', linestyle='dashed', zorder=0)
        if plot_extra is not None: plot_extra(plot_name, ax, mask_dict, time_init=time_init) 
        plot_path = os.path.join(plot_dir, dateStr+output_name_tag+plot_name+'.'+plot_format)
        plt.savefig(plot_path, format=plot_format)
        if verbosity >= 1: print('    saving:', plot_path)
        plt.close('all')
    #####################################################################################################################
    ###Wall Spatial Distribution###############################################################################
    s1_corrF_LYcorr_wallMask  = (chi2DisplayRange[0] <= eventData['multipeak_chi2_ndof'])
    s1_corrF_LYcorr_wallMask &= (eventData['multipeak_chi2_ndof'] < chi2DisplayRange[1])
    s1_corrF_LYcorr_wallMask &= (timeDisplay[0] <= eventData['time']*1.0E-9)
    s1_corrF_LYcorr_wallMask &= (eventData['endtime']*1.0E-9 < timeDisplay[1])
    s1_corrF_LYcorr_wallMask &= (s1_corrF_LYcorrDisplay_forTheta[0] <= eventData['s1_corrF_LYcorr'])
    s1_corrF_LYcorr_wallMask &= (eventData['s1_corrF_LYcorr'] < s1_corrF_LYcorrDisplay_forTheta[1])
    nullMask = np.ones(len(s1_corrF_LYcorr_wallMask), dtype=bool)
    if 'mask_wall_save_data' in mask_dict:
        plot_name = 's1only_surface_rate_estimation_5s1corr'
        if save_plot == True:
            plt, ax = plotStandard_s1only_highA(eventData, nullMask, s1_corrF_LYcorr_wallMask,\
                                                binN=generalBinN, ampKey='s1_corrF_LYcorr', r2Key='r2_corr',\
                                                r2Range=r2_corrRange, zKey=zName, zRange=zRange, verbosity=verbosity)
            if plot_extra is not None: plot_extra(plot_name, ax, mask_dict)
            plot_path = os.path.join(plot_dir, dateStr+output_name_tag+plot_name+'.'+plot_format)
            plt.savefig(plot_path, format=plot_format)
            if verbosity >= 1: print('    saving:', plot_path)
            plt.close('all')
        plot_name = 's1only_surface_rate_estimation_5spatial'
        if save_plot == True:
            gridSpec = [1, 1]
            figSize, marginRatio = getSizeMargin(gridSpec, subplotSize=[18.0, 8.0], plotMargin=[0.07, 0.1, 0.07, 0.02])
            fig = matplotlib.pyplot.figure(figsize=figSize); fig.subplots_adjust(*marginRatio)
            gs = gridspec.GridSpec(*gridSpec)
            matplotlib.rc('xtick', labelsize=18)
            matplotlib.rc('ytick', labelsize=18)
            ax = list()
            for axIdx in range(gs.nrows*gs.ncols):
                ax.append(fig.add_subplot(gs[axIdx]));
                ax[-1].ticklabel_format(style='sci', scilimits=(-3, 3))
                ax[-1].xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
                ax[-1].yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())    
            ######
            axIdx = 0
            thetaRange = [-np.pi, np.pi]
            rangeXYZ = [thetaRange, zRange]
            hist = ax[axIdx].hist2d(eventData[thetaName][s1_corrF_LYcorr_wallMask],\
                                    eventData[zName][s1_corrF_LYcorr_wallMask], bins=[3*generalBinN, generalBinN],\
                                    range=[rangeXYZ[0], rangeXYZ[1]], cmin=theta_cmin, norm=LogNorm(), cmap=CMAP)
            divider = make_axes_locatable(ax[axIdx])
            cax = divider.append_axes('right', size='5%', pad=0.1)
            cbar = fig.colorbar(hist[3], ax=ax[axIdx], cax=cax)
            ax[axIdx].get_yaxis().set_visible(False)
            ax[axIdx].set_xlabel(thetaName, fontsize=18)
            ax[axIdx].set_xlim(*rangeXYZ[0])
            ax[axIdx].set_ylim(*rangeXYZ[1])
            plot_pillars_theta(ax[axIdx], zorder=10)
            ###
            axUpper = divider.append_axes('top', size='60%', pad='6%')
            ax[axIdx].figure.add_axes(axUpper)
            hist = axUpper.hist(eventData[thetaName][s1_corrF_LYcorr_wallMask], bins=3*generalBinN,\
                                range=rangeXYZ[0], color='steelblue', linewidth=1)
            axUpper.set_ylabel('counts', fontsize=18)
            axUpper.get_xaxis().set_visible(False)
            axUpper.set_xlim(*rangeXYZ[0])
            plot_pillars_theta(axUpper, zorder=10)
            ###
            axSide = divider.append_axes('left', size='15%', pad='2%')
            ax[axIdx].figure.add_axes(axSide)
            axSide.set_xlabel('counts', fontsize=18)
            axSide.set_ylabel(zName, fontsize=18)
            axSide.ticklabel_format(style='sci', scilimits=(-2, 2))
            axSide.xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
            axSide.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
            hist = axSide.hist(eventData[zName][s1_corrF_LYcorr_wallMask], bins=generalBinN,\
                               range=rangeXYZ[1], color='steelblue', linewidth=1, orientation='horizontal')
            axSide.set_ylim(*rangeXYZ[1])
            xRange = axSide.get_xlim()
            axSide.set_xlim(xRange[1], xRange[0])
           
            if plot_extra is not None: plot_extra(plot_name, ax, mask_dict)
            plot_path = os.path.join(plot_dir, dateStr+output_name_tag+plot_name+'.'+plot_format)
            matplotlib.pyplot.savefig(plot_path, format=plot_format)
            if verbosity >= 1: print('    saving:', plot_path)
            matplotlib.pyplot.close('all')
    #####################################################################################################################
    rateDict = {}
    if time_init is None:
        time_init = 1.0E-9*procData['time'][mask_final_FV][0]
        time_init = toBeginningOfDay(time_init)
    rateDict['time_init'] = time_init
    rateDict['time_fin']  = time_fin
    rateDict['runChunks']       = runChunks
    rateDict['runModes']        = runModes
    rateDict['timestampRanges'] = timestampRanges
    rateDict['liveTimes']       = liveTimes
    rateDict['dayArr']          = dayArr
    rateDict['dayErr']          = []
    for timestampRangeChunk in timestampRanges:
        rateDict['dayErr'].append(0.5*(timestampRangeChunk[-1][1] - timestampRangeChunk[0][0])/(24*3600.0))
    rateDict['eventData']   = eventData
    rateDict['eventCounts'] = eventCounts 
    rateDict['eventRates']  = eventRates
    rateDict['cali_factor_ref']        = cali_factor
    rateDict['fitLightYieldCorrPeaks'] = fitLightYieldCorrPeaks
    rateDict['fitRanges_gainCorr']     = fitRanges_gainCorr
    rateDict['peakFitPars']            = peakFitPars
    rateDict['relative_light_yields'] = relative_light_yields
    rateDict['corr_factors']          = corr_factors
    rateDict['chi2_ndofs']            = chi2_ndofs
    rateDict['eventRates_fitA'] = eventRates_fitA
    rateDict['eventRates_fullFid_fitA'] = {}
    for peakName in rateDict['eventRates_fitA']:
        rateDict['eventRates_fullFid_fitA'][peakName] = []
        for eventRate in rateDict['eventRates_fitA'][peakName]:
            eventRate_fullFid = ufloat(-1, 0)
            if eventRate > 0: eventRate_fullFid = removeDivisionErr(eventRate, fidMassPar['XeFidZ_wall'])
            rateDict['eventRates_fullFid_fitA'][peakName].append(eventRate_fullFid)
    rateDict['correction_derivation_parameters']  = corrParPickle
    rateDict['estimate_fiducial_mass_parameters'] = fidMassPar
    rateDict['parameters'] = {}
    for key in inputPars: rateDict['parameters'][key] = inputPars[key]
    for key in mask_dict:
        if ('mask' not in key) and ('Xe' not in key): rateDict['parameters'][key] = mask_dict[key]
    rateDict['parameters']['XeFidZ_ref_cm']    = fidMassPar['XeFidZ_ref']
    rateDict['parameters']['XeFidR2_ref_cm^2'] = fidMassPar['XeFidR2_ref']
    rateDict['parameters']['XeFidZ_cm']        = fidMassPar['XeFidZ']
    rateDict['parameters']['XeFidZ_wall_cm']   = fidMassPar['XeFidZ_wall']
    rateDict['parameters']['XeFidR2_cm^2']     = fidMassPar['XeFidR2']
    with open(output_path, 'wb') as handle:
        pickle.dump(rateDict, handle, protocol=pickle.HIGHEST_PROTOCOL)
    if verbosity >=1: print('    saving:', output_path)

    rateDictLite = ufloatToErrDict(rateDict)
    del rateDictLite['eventData']
    with open(output_path.replace('.pkl', '_lite.pkl'), 'wb') as handle:
        pickle.dump(rateDictLite, handle, protocol=pickle.HIGHEST_PROTOCOL)
    if verbosity >=1: print('    saving:', output_path.replace('.pkl', '_lite.pkl'))
    '''
#########################################################################################################################
# helper functions
def load_context_with_run_info_s1only(self):
    import strax, straxen, cutax
    if self.verbosity >= 1: straxen.print_versions(('strax', 'straxen', 'cutax', 'xedocs', 'rframe'))
    #####################################################################################################################
    @strax.takes_config(strax.Option('area_hard_high_thres', default=2.0E4,\
                        help='To reveal only high energy events'),)
    class cut_area_hard_high_thres(strax.CutPlugin):
        depends_on = ('peak_basics',)
        provides = 'cut_area_hard_high_thres'
        cut_name = 'cut_area_hard_high_thres'
        cut_description = 'To reveal only high energy events'
        __version__ = '0.0.1'
        def cut_by(self, peaks):
            greater_than_thres = peaks['area'] > self.config['area_hard_high_thres']
            return greater_than_thres
    class naive_peak_positions(strax.Plugin):
        depends_on = ('peaks',)
        provides = 'naive_peak_positions'
        dtype = [('x_naive', np.float32, 'Naive X position (cm)'),
                 ('y_naive', np.float32, 'Naive Y position (cm)')]
        dtype.extend(strax.time_fields)
        __version__ = '0.0.1'
        def compute(self, peaks):
            weight_map  = np.abs(peaks['area_per_channel'])
            weight_map /= np.sum(weight_map, axis=1)[:, None]
            result = np.ones(len(peaks), self.dtype)   
            result['time']    = peaks['time'] 
            result['endtime'] = strax.endtime(peaks)
            result['x_naive'] = np.sum(weight_map*straxen.pmt_positions().x.values, axis=1)
            result['y_naive'] = np.sum(weight_map*straxen.pmt_positions().y.values, axis=1)
            return result

    self.register_lists = []                                                                        #plugin registers
    try:
        self.register_lists.append(cutax.plugins.deadtime)
    except:
        pass
    self.register_singles    = [naive_peak_positions, cut_area_hard_high_thres]
    self.targets_fixed       = ['cut_area_hard_high_thres']
    self.selection_str_fixed = ['cut_area_hard_high_thres']
    #####################################################################################################################
    context_args = self.config['general']['context']
    context_func = getattr(cutax, context_args['context'])
    self.context = context_func(**{key:context_args[key] for key in context_args\
                                   if (key not in ['context', 'storage_append'])})
    if 'storage_append' in context_args:
        for storage_dir in context_args['storage_append']:
            self.context.storage.append(strax.DataDirectory(storage_dir, readonly=True)) 
    for register_list   in self.register_lists:   self.context.register_all(register_list)
    for register_single in self.register_singles: self.context.register(register_single)
    if self.verbosity >= 1:
        print('\nst.storage', self.context.storage)
        print('lineage:', self.context.key_for('00', self.config['general']['data_df']['targets'][0]))
    if self.verbosity >= 2: 
        pluginKeys = sorted(list(self.context._plugin_class_registry.keys()))
        for pluginKey in pluginKeys: print('   ', pluginKey, self.context._plugin_class_registry[pluginKey].__version__)

    self.run_info_available_target = None
    if 'run_info_available_target' in self.config['general']:
        self.run_info_available_target = self.config['general']['run_info_available_target']
    else: self.run_info_available_target = self.config['general']['data_df']['targets'][0] 
    run_info_args = self.config['general']['run_info']
    self.run_info = self.context.select_runs(available=self.run_info_available_target,\
                                             **{key:run_info_args[key] for key in run_info_args})
    lineage_hash = [self.context.key_for(run_id, self.run_info_available_target) for run_id in self.run_info['name']]
    self.run_info['lineage_hash'] = np.array(lineage_hash)
    if self.verbosity >= 2: print('  loading run_info:\n', self.run_info.head(10), '\n', self.run_info.tail(10))
 
    self.python_module = None
    if os.path.isfile(self.config['general']['python_config_file']) == True:
        python_config_file = self.config['general']['python_config_file'].replace('.py', '')
        python_config_path = os.path.join(*python_config_file.split('/')[:-1])
        if python_config_file[0] == '/': python_config_path = '/' + python_config_path
        sys.path.append(python_config_path)
        self.python_module = importlib.import_module(python_config_file.split('/')[-1])
        if hasattr(self.python_module, 'runInfoDF_selection'):
            self.run_info = self.python_module.runInfoDF_selection(self.run_info)
        else: print('  Note: runInfoDF_selection() does not exist in general:python_config_file')
    else: print('  Note: general:python_config_file not set')
    if len(self.run_info.index) == 0: 
        raise AssertionError('load_context_with_run_info_s1only(): no runs have been found in self.run_info, '+\
                             'check runInfoDF_selection function from:\n  '+self.config['general']['python_config_file'])
    self.run_info.reset_index(drop=True, inplace=True)
    if self.verbosity >= 1: 
        print('  available modes:', set(self.run_info['mode']), '\n  leftover tags:', set(self.run_info['tags']))

    self.targets = copy.copy(self.config['general']['data_df']['targets'])
    self.targets.extend([target for target in self.targets_fixed if target not in self.targets])
    self.selection_str = copy.copy(self.config['general']['data_df']['selection_str'])
    self.selection_str.extend([selection for selection in self.selection_str_fixed\
                               if selection not in self.selection_str]) 
def select_runs_anode_voltage(run_ids, runInfoData, context, selectLow=None, anodeThres=None, verbosity=None):
    import strax, straxen
    if selectLow  is None: selectLow = True                                         #True for s1only; False for s1s2
    if anodeThres is None: anodeThres = 100.0
    if verbosity  is None: verbosity = 1
    slowControl = straxen.SCADAInterface(context, use_progress_bar=(1<=verbosity))  #NOTE: take resouces from SC server 
    scPars = {'anode0': 'XE1T.CTPC.Board14.Chan001.VMon', 'anode1': 'XE1T.CTPC2.Board00.Chan000.VMon'}
    scDF = slowControl.get_scada_values(scPars, run_id=list(run_ids), every_nth_value=200, fill_gaps="forwardfill",\
                                        query_type_lab=False)
    scDF = scDF.reset_index()
    scData = scDF.to_dict('series')
    for key in scData: scData[key] = np.array(scData[key])
    scData['time']  = np.array([pandas_datetime.to_numpy() for pandas_datetime in scData['time UTC']])
    scData['anode'] = np.maximum(scData['anode0'], scData['anode1'])
    if verbosity >= 1: print(scDF.head(10))
    
    new_run_ids = []
    for run_id in run_ids:
        runInfoMask = (runInfoData['run_number'] == run_id)
        scMask  = (runInfoData['start'][runInfoMask][0] <= scData['time']) 
        scMask &= (scData['time'] < runInfoData['end'][runInfoMask][0])
        if len(scData['anode'][scMask]) >= 2:
            if ((selectLow == True)  and (anodeThres <= np.mean(scData['anode'][scMask]))) or\
               ((selectLow == False) and (np.mean(scData['anode'][scMask] < anodeThres))):
                new_run_ids.append(run_id)
        else: warnings.warn('select_runs_anode_voltage(): run '+str(run_id)+' removed for having no SC data', Warning)
    return np.array(new_run_ids)
#########################################################################################################################
def main():
    args        = ARGS
    config_dict = CONFIG_DICT
    config_dict['general']['config_file']   = args.config_file
    config_dict['general']['analysis_mode'] = copy.copy(args.analysis_mode)
    if 'load_deadtime' not in config_dict['general']: config_dict['general']['load_deadtime'] = True
    if 'deadtime_type' not in config_dict['general']: config_dict['general']['deadtime_type'] = 'deadtime_tpc'
    if 'generalBinN'   not in config_dict['general']: config_dict['general']['generalBinN']   = 400
    config_dict = trim_dict_entries(config_dict)
    for editStr in args.config_edit:
        keyList = editStr.split(':')
        values = keyList[-1].split('=')[1].split(',')
        output_value = []
        for value_ in values:
            value = copy.copy(value_)
            valueCheck = value.replace('-', '').replace('.', '').replace('e', '').replace('E', '')
            if valueCheck.isdigit() == True:
                if ('.' in value) or ('e' in value.lower()): value = float(value)
                else:                                        value = int(value)
            else:
                value = value.replace(' ', '').replace('\'', '')
                if   value.lower() == 'true':  value = True
                elif value.lower() == 'false': value = False
            output_value.append(value)
        if len(output_value) == 1: output_value = copy.copy(output_value[0])
        keyList[-1] = keyList[-1].split('=')[0]
        config_dict_pointer = config_dict
        for key in keyList[:-1]: config_dict_pointer = config_dict_pointer[key]
        config_dict_pointer[keyList[-1]] = output_value
    if 'selection_str' not in config_dict['general']['data_df']:
        config_dict['general']['data_df']['selection_str'] = config_dict['general']['data_df']['selection']
    output_dir = config_dict['general']['output_dir']
    config_dict = replace_dict_string_string(config_dict, '$output_dir', output_dir)
    config_dict = replace_dict_string_value( config_dict, 'None',        None)
    verbosity     = config_dict['general']['verbosity']
    analysis_mode = config_dict['general']['analysis_mode']
    analysis_main_steps_to_process = []
    for key in config_dict['analysis_main']['process_step']:
        if config_dict['analysis_main']['process_step'][key] == True: analysis_main_steps_to_process.append(key)

    ###
    if verbosity >= 1: 
        print('Reading config_file:', args.config_file, '\nStarting radon analysis...')
        print('  input arguments:', args)  
    run_args = copy.deepcopy(args)
    run_args.runlist = str(run_args.runlist)
    run_args.runlist = run_args.runlist.replace(' ', '').replace('[', '').replace(']', '').replace('\'', '')
    if isinstance(run_args.config_edit, list) == True:
        if len(run_args.config_edit) > 0:
            config_edit_str = ''
            for config_item in run_args.config_edit:
                config_item = config_item.replace(' ', '').replace('\'', '')
                config_edit_str += str(config_item) + '+++'
            run_args.config_edit = config_edit_str[:-3]
    if config_dict['data_building']['process_step'] == True:
        if verbosity >= 1: print('Proceeding with step: data_building')
        data_builder = RadonAnalysisDataBuilder(config_dict=config_dict, run_list=args.runlist,\
                                                load_func=load_context_with_run_info_s1only,\
                                                log_dir=config_dict['general']['submit']['log'], output_dir=output_dir)
        if analysis_mode == 'submit': 
            submit_args                  = copy.deepcopy(config_dict['general']['submit'])
            submit_args['argument_list'] = copy.deepcopy(vars(run_args))
            del submit_args['argument_list']['analysis_mode'], submit_args['argument_list']['runlist'],submit_args['log']
            if 'cutax_path' in submit_args: del submit_args['cutax_path']
            if submit_args['argument_list']['config_edit'] == []: submit_args['argument_list']['config_edit']  = ''
            else:                                                 submit_args['argument_list']['config_edit'] += '+++'
            for key in config_dict['analysis_main']['process_step']:
                submit_args['argument_list']['config_edit'] += 'analysis_main:process_step:'+key+'=false+++'
            submit_args['argument_list']['config_edit'] = submit_args['argument_list']['config_edit'][:-3]
            if submit_args['container'] == 'CONDA_DEFAULT_ENV':
                containerPath = str(os.environ['CONDA_DEFAULT_ENV'])
                containerPath = containerPath.replace('XENONnT_',\
                                                      '/project2/lgrandi/xenonnt/singularity-images/xenonnt-')+'.simg'
                submit_args['container'] = containerPath
            submit_args['jobname'] += (date.today().strftime('%y%m%d') + '_data_builder_s1only')
            submit_args['jobname'] += config_dict['analysis_main']['output_name_tag']
            if verbosity >= 1: print('  submission parameters:', submit_args)
            data_builder.submit(**{key:submit_args[key] for key in submit_args})
        elif analysis_mode == 'run': data_builder.run( Namespace(**{'runlist':data_builder.get_runlist()}) )
    if len(analysis_main_steps_to_process) > 0:
        analysis_main = S1onlyRadonAnalysisMain(config_dict=config_dict, run_list=args.runlist,\
                                                load_func=load_context_with_run_info_s1only,\
                                                log_dir=config_dict['general']['submit']['log'], output_dir=output_dir)
        if analysis_mode == 'submit': 
            submit_args                  = copy.deepcopy(config_dict['general']['submit'])
            submit_args['argument_list'] = copy.deepcopy(vars(run_args))
            del submit_args['argument_list']['analysis_mode'], submit_args['log']
            if 'cutax_path' in submit_args: del submit_args['cutax_path']
            if submit_args['argument_list']['config_edit'] == []: submit_args['argument_list']['config_edit']  = ''
            else:                                                 submit_args['argument_list']['config_edit'] += '+++' 
            submit_args['argument_list']['config_edit'] += 'data_building:process_step=false'
            if submit_args['container'] == 'CONDA_DEFAULT_ENV':
                containerPath = str(os.environ['CONDA_DEFAULT_ENV'])
                containerPath = containerPath.replace('XENONnT_',\
                                                      '/project2/lgrandi/xenonnt/singularity-images/xenonnt-')+'.simg'
                submit_args['container'] = containerPath
            submit_args['jobname'] += (date.today().strftime('%y%m%d') + '_analysis_main_s1only_')
            submit_args['jobname'] += config_dict['analysis_main']['output_name_tag']
            if verbosity >= 1: print('  submission parameters:', submit_args)
            analysis_main.submit(**{key:submit_args[key] for key in submit_args})
        elif analysis_mode == 'run': analysis_main.run(run_args)
    if verbosity >= 1: print('End of code')
#########################################################################################################################
if __name__ == '__main__': main()









