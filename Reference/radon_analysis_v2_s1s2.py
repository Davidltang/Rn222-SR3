#!/usr/bin/env python
# -*- coding: UTF-8 -*- 
from _GlobalFuncs import *
PARSER = ArgumentParser() #analysis_mode and --runlist are required by nton/analysis.py
PARSER.add_argument('analysis_mode', type=str, choices=['submit', 'run'],
                    help='analysis script mode:\n'+\
                         '  submit - submit the main analysis as a job to the cluster\n'+\
                         '  run - run the analysis locally')
PARSER.add_argument('--runlist', nargs='*', action='store', type=str, default=['None'],\
                    help=('list of runs to process; overwrites run_list_file from the configuration files; '+\
                          'load_pickle_data needs to be set false'))
PARSER.add_argument('--config_file', type=str, default='configuration_files/radon_analysis_config_v2_s1s2.json',
                    help='path to analysis configuration file.')
PARSER.add_argument('--config_edit', action='append', default=[],\
                    help='editing configuration file; can be repeated.\n  format: \"class:subclass...:value\"')
ARGS = PARSER.parse_args()
if ARGS.runlist[0] == 'None': ARGS.runlist = None
else:
    if (len(ARGS.runlist) == 1) and (',' in ARGS.runlist[0]): ARGS.runlist = ARGS.runlist[0].split(',')
    ARGS.runlist = [int(item) for item in ARGS.runlist]
ARGS.config_edit = [item.replace('[', '').replace(']', '') for item in ARGS.config_edit]
config_edit_list = []
if (isinstance(ARGS.config_edit, list)) and (len(ARGS.config_edit) > 0):
    for config_edit_item in ARGS.config_edit: config_edit_list.extend(config_edit_item.split('+++'))
    ARGS.config_edit = config_edit_list
CONFIG_DICT = json.load(open(ARGS.config_file, 'r'))
if (CONFIG_DICT['general']['submit']['container'] != 'CONDA_DEFAULT_ENV') and (ARGS.analysis_mode == 'run'):
    if 'cutax_path' in CONFIG_DICT['general']['submit']:
        if os.path.exists(CONFIG_DICT['general']['submit']['cutax_path']):
            sys.path.insert(0, CONFIG_DICT['general']['submit']['cutax_path'])
        else:
            warnings.warn('\n\ncutax_path does not exist: '+str(CONFIG_DICT['general']['submit']['cutax_path'])+'\n',\
                          Warning)
from job_submitter.analysis import Analysis, AnalysisRunBatch
from job_submitter.utils import run_numbers_to_run_names

#########################################################################################################################
class S1S2RadonAnalysisMain(Analysis):
    """ Manage main radon analysis steps. """
    def __init__(self, config_dict, run_list, load_func, **kwargs):
        super().__init__(**kwargs)
        #self.logging.setLevel(logging.DEBUG)
        radon_analysis_class_initialization(self, config_dict, run_list, load_func, **kwargs)
        if (self.run_list is None) and (self.config['analysis_main']['load_from_run_list'] == True):
            self.run_list = read_run_list_file(self.config['analysis_main']['run_list_file'])
        self.step_functions = {'basic_rate_estimation':       self._basic_rate_estimation,
                               'correction_derivation':       self._correction_derivation,
                               'peak_region_rate_estimation': self._peak_region_rate_estimation,
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
        else: print('\n\nNote: general:python_config_file not set\n\n')
    def run(self, args):
        """ Run requested analysis step. """
        for analysis_step_name in self.step_functions:
            if self.config['analysis_main']['process_step'][analysis_step_name] == True: 
                if self.verbosity >= 1: print('Proceeding with step: '+analysis_step_name)
                self.step_functions[analysis_step_name](args)
    def _load_data(self, args, rateDict_path=None):
        """ Manage data loading. """
        updateNew, preSubN = False, 0
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
            _, runInfoData, _, procData = loadPickleRunEventData(self.data_dirs, iter_range=iter_range,\
                run_exclude=run_exclude, run_include=run_include, timestamp_range=timestamp_range,\
                runInfoDF_selection=runInfoDF_selection, run_chunks=existingRunChunks,\
                procDF_selection=procDF_selection, verbosityLevel=(self.verbosity-1)*10)
        else:
            self.load_context_with_run_info(self)
            if updateNew == True:
                if (existingRunChunks is not None):
                    iter_range = [-np.sum(np.array(self.run_info['number']) > int(existingRunChunks[-1][-1]))-1, None]
                    iter_range[0] -= preSubN
                    if self.verbosity >= 1: print('\n\n_load_data(): updateNew => iter_range =', iter_range, '\n\n') 
                else:print('\n\nNote: updateNew does nothing because there are no pickle file to extract runs from\n\n')
            if iter_range[0] is None:
                iter_range[0] = -np.sum( np.array(self.run_info['number']) > int(existingRunChunks[-1][-1]) ) - 1
            runInfoDF, runInfoData, procDF, procData = load_run_event_data(\
                self, iter_range=iter_range, run_exclude=run_exclude,\
                run_include=run_include, timestamp_range=timestamp_range, runInfoDF_selection=runInfoDF_selection,\
                run_chunks=existingRunChunks, verbosityLevel=(self.verbosity-1)*10)
        if self.verbosity >= 1: print('  generating user defined data via gen_procData...')
        self.python_module.gen_procData(procData)
        return runInfoData, procData
    def _basic_rate_estimation(self, args):
        self.config['analysis_main']['basic_rate_estimation']['deadtime_type'] = self.config['general']['deadtime_type']
        self.config['analysis_main']['basic_rate_estimation']['generalBinN']   = self.config['general']['generalBinN']
        runInfoData, procData = self._load_data(args)
        mask_dict = self.python_module.gen_masks(procData, verbosity=self.verbosity)
        basic_rate_estimation(runInfoData, procData, mask_dict, output_dir=self.result_dir,\
                              plot_dir=self.plot_dir, save_plot=self.save_plot, plot_format=self.plot_format,\
                              output_name_tag=self.config['analysis_main']['output_name_tag'],\
                              config_dict=self.config['analysis_main']['basic_rate_estimation'],\
                              plot_extra=self.python_module.plot_extra, verbosity=self.verbosity)
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
    def _peak_region_rate_estimation(self, args):
        """ Estimate the alpha rate through the counting within peak regions. """
        self.config['analysis_main']['peak_region_rate_estimation']['deadtime_type'] = \
            self.config['general']['deadtime_type']
        self.config['analysis_main']['peak_region_rate_estimation']['generalBinN'] =self.config['general']['generalBinN']
        efficiency_dict = {}
        if self.config['analysis_main']['peak_region_rate_estimation']['apply_alpha_efficiency_correction'] == True:
            alpha_eff_dict = self.python_module.gen_efficiency_reference(verbosity=self.verbosity)
            for key in alpha_eff_dict: efficiency_dict[key] = alpha_eff_dict[key]
        rateDict_path = os.path.join(self.result_dir,\
                                     self.config['analysis_main']['output_name_tag']+'peak_region_rate_estimation.pkl')
        runInfoData, procData = self._load_data(args, rateDict_path=rateDict_path)
        s1corr_paths = self.config['analysis_main']['peak_region_rate_estimation']['s1corr_paths']
        corrParPickle = None
        for s1corrPath in s1corr_paths:
            with open(s1corrPath, 'rb') as handle:
                corrParPickle = pickle.load(handle)
        procData['s1_corr1_AFT']  = np.array(procData['s1_area'], dtype='float64')
        procData['s1_corr1_AFT'] *= corrParPickle['S1_areaCorr1Pars_AFT'][-1]
        procData['s1_corr1_AFT'] /= exp2Func(procData['s1_area_fraction_top'], *corrParPickle['S1_areaCorr1Pars_AFT'])
        procData['s1_corr1_DT']   = np.array(procData['s1_area'], dtype='float64')
        procData['s1_corr1_DT']  *= corrParPickle['S1_areaCorr1Pars_DT'][-1]
        procData['s1_corr1_DT']  /= exp2Func(procData['drift_time'], *corrParPickle['S1_areaCorr1Pars_DT'])

        procData['s1_corr1']  = procData['s1_corr1_AFT']
        procData['s1_corrF']  = np.array(procData['s1_corr1'], dtype='float64')
        procData['s1_corrF'] *= corrParPickle['S1_areaCorr2Pars'][0]
        procData['s1_corrF'] /= linearFunc(procData['r2'], *corrParPickle['S1_areaCorr2Pars'])
        mask_dict = self.python_module.gen_masks(procData, verbosity=self.verbosity)
        peak_region_rate_estimation(runInfoData, procData, efficiency_dict, mask_dict, s1corr_paths,\
                                    output_dir=self.result_dir, plot_dir=self.plot_dir, save_plot=self.save_plot,\
                                    plot_format=self.plot_format,\
                                    output_name_tag=self.config['analysis_main']['output_name_tag'],\
                                    config_dict=self.config['analysis_main']['peak_region_rate_estimation'],\
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
def basic_rate_estimation(runInfoData, procData, mask_dict, output_dir=None, plot_dir=None, save_plot=None,\
                          plot_format=None, output_name_tag=None, config_dict=None, plot_extra=None, verbosity=None):
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

    deadtime_type    = None
    generalBinN      = 400
    chunkCheckTags   = False
    combTimeThres    = 6*3600.0
    s1_ampKey        = 'cs1'
    s1_ampBinN       = 800
    s1_ampRange      = [0.0, 1.0E5]
    dateDisplayRange = ['-inf', 'inf']
    if config_dict is not None:
        if 'deadtime_type'     in config_dict: deadtime_type    = config_dict['deadtime_type']
        if 'generalBinN'       in config_dict: generalBinN      = config_dict['generalBinN']
        if 'chunkCheckTags'    in config_dict: chunkCheckTags   = config_dict['chunkCheckTags']
        if 'combTimeThres'     in config_dict: combTimeThres    = config_dict['combTimeThres']
        if 's1_ampKey'         in config_dict: s1_ampKey        = config_dict['s1_ampKey']
        if 's1_ampBinN'        in config_dict: s1_ampBinN       = config_dict['s1_ampBinN']
        if 's1_ampRange'       in config_dict: s1_ampRange      = config_dict['s1_ampRange']
        if 'dateDisplayRange'  in config_dict: dateDisplayRange = config_dict['dateDisplayRange']
    #####################################################################################################################
    mask0 = copy.deepcopy(mask_dict['mask_base'])
    mask1 = copy.deepcopy(mask_dict['mask_final_FV'])
    #####################################################################################################################
    if verbosity >= 1: print('  base data:')
    plot_name = 'basic_rate_estimation_0base' 
    if save_plot == True:
        plt, ax = plotStandard_s1s2_highA(procData, mask0, mask1, s1ampKey=s1_ampKey, s1ampRange=s1_ampRange,\
                                          verbosity=verbosity, binN=generalBinN)
        plot_path = os.path.join(plot_dir, dateStr+output_name_tag+plot_name+'.'+plot_format)
        plt.savefig(plot_path, format=plot_format)
        if verbosity >= 1: print('    saving:', plot_path)
        plt.close('all')
    
    #####
    # Rate evaluation
    #####################################################################################################################
    plot_name = 'basic_rate_estimation_1rate'
    if len(procData['run_number'][mask0]) == 0: raise AssertionError('mask0 cuts out all events')
    time_init = toBeginningOfDay(1.0E-9*procData['time'][mask0][0])
    time_fin  = 1.0E-9*procData['endtime'][mask0][-1]
    if verbosity >= 1: 
        print('  rate data:')
        print('    XeDensity   =', mask_dict['XeDensity'], 'kg/cm^3')
        print('    XeFidMass   =', mask_dict['XeFidMass'], 'kg')    

    fitRangeSetPickle = ['mask_final_FV']

    runChunks, runModes = genRunChunks(np.unique(procData['run_number']), runInfoData, combTimeThres,\
                                       checkTags=chunkCheckTags, verbosity=(verbosity-1))
    timestampRanges, liveTimes, dayArr = [], [], []
    eventDataS1ampArr, eventDataS1ampArr_noBiPo, eventData_chunk = [], [], []
    eventCounts = {peakName:[] for peakName in fitRangeSetPickle}
    eventRates  = {peakName:[] for peakName in fitRangeSetPickle}
    if (deadtime_type is not None) and ((deadtime_type in runInfoData) == False): deadtime_type == None
    for runChunk in tqdm(runChunks):
        timestampRangeChunk, liveTimeChunk = [], ufloat(0, 0)
        mask_run_sel = (procData['run_number'] == -1)
        for run_number in runChunk:      
            timestampRangeChunk.append([runInfoData['start_timestamp'][runInfoData['name']==run_number][0],
                                        runInfoData['end_timestamp'][runInfoData['name']==run_number][0]])
            liveTime  = pd.Timedelta(runInfoData['livetime'][runInfoData['name']==run_number][0]).total_seconds()
            if deadtime_type is not None: liveTime-=runInfoData[deadtime_type][runInfoData['name']==run_number][0]*1.0E-9
            else: warnings.warn('basic_rate_estimation(): deadtime plugin NOT implemented', Warning) 
            liveTimeChunk += ufloat(liveTime, 0.5)          #NOTE: 0.5 gps uncertainty for livetime
            mask_run_sel |= (procData['run_number'] == run_number)
        timestampRanges.append(timestampRangeChunk)
        liveTimes      .append(liveTimeChunk)
        dayArr         .append((timestampRangeChunk[0][0] - time_init)/(3600.0*24.0))
        
        eventDataS1amp = np.array(procData[s1_ampKey][mask1 & mask_run_sel], dtype='float64')
        eventDataS1ampArr.append(eventDataS1amp)
        eventData_chunk.append({key:procData[key][mask1 & mask_run_sel] for key in procData})
        if 'mask_BiPoRemoval' in mask_dict:
            eventDataS1amp_noBiPo = np.array(procData[s1_ampKey][mask1 & mask_run_sel & mask_dict['mask_BiPoRemoval']],\
                                             dtype='float64')
            eventDataS1ampArr_noBiPo.append(eventDataS1amp_noBiPo) 
        countDataSet = {}
        for peakName in fitRangeSetPickle:
            peakData = eventDataS1amp
            countDataSet[peakName] = len(peakData)
            eventCounts[peakName].append(len(peakData))
            activity = 1.0E6*ufloat(countDataSet[peakName], np.sqrt(countDataSet[peakName]))
            activity /= (liveTimeChunk*mask_dict['XeFidMass'])
            eventRates[peakName].append(activity)
 
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
        hist = ax[axIdx].hist(procData[s1_ampKey][mask1], bins=s1_ampBinN, range=s1_ampRange, color='steelblue',\
                              linewidth=1)
        ax[axIdx].set_xlabel(s1_ampKey, fontsize=18)
        ax[axIdx].set_ylabel('counts', fontsize=18)
        ax[axIdx].set_xlim(s1_ampRange)
        ax[axIdx].set_ylim(bottom=0.5)
        ax[axIdx].set_yscale('log')

        for peakName in fitRangeSetPickle:
            hist = ax[axIdx].hist(procData[s1_ampKey][mask0], bins=s1_ampBinN, range=s1_ampRange,\
                                  color='gainsboro', linewidth=1, alpha=0.3)
        ######
        axIdx = 1
        tracePlot, tracePlotText = {}, {}
        for peakName in fitRangeSetPickle:
            tracePlot[peakName] = ax[axIdx].errorbar(dayArr, [item.n for item in eventRates[peakName]],\
                                                             [item.s for item in eventRates[peakName]],\
                                                     color='blue', fmt='s', **ERROR_PLOT_PARS)[0]
            tracePlotText[peakName]  = peakName
        setXaxisDateLabel(ax[axIdx], time_init, [time_init, time_fin])
        ax[axIdx].set_ylabel('activity concentration (uBq/kg_Xe)', fontsize=18)
        ax[axIdx].set_ylim(bottom=3.0E-2)
        ax[axIdx].set_yscale('log')
        ax[axIdx].secondary_yaxis('right')
        if plot_extra is not None: plot_extra(plot_name, ax, mask_dict) 
        plot_path = os.path.join(plot_dir, dateStr+output_name_tag+plot_name+'.'+plot_format)
        plt.savefig(plot_path, format=plot_format)
        if verbosity >= 1: print('    saving:', plot_path)
        plt.close('all')
    
    ###########################################################
    output_path = os.path.join(output_dir, output_name_tag+'basic_rate_estimation.pkl')
    RnRatePickle, inputDictToExtend = None, None    
    inputDictToExtend = {'runChunks': runChunks, 'runModes': runModes, 'timestampRanges': timestampRanges,\
                         'liveTimes': liveTimes, 'eventCounts': eventCounts, 'eventRates': eventRates}
    RnRatePickle = loadPickleRateData(output_path, dayArr, time_init, inputDictToExtend=inputDictToExtend)
    eventData = {}
    if RnRatePickle is not None:
        time_init = copy.deepcopy(RnRatePickle['time_init'])
        time_fin  = copy.deepcopy(RnRatePickle['time_fin'])
        runChunks       = copy.deepcopy(RnRatePickle['runChunks'])
        runModes        = copy.deepcopy(RnRatePickle['runModes'])
        timestampRanges = copy.deepcopy(RnRatePickle['timestampRanges'])
        liveTimes       = copy.deepcopy(RnRatePickle['liveTimes'])
        dayArr          = copy.deepcopy(RnRatePickle['dayArr'])
        eventCounts       = copy.deepcopy(RnRatePickle['eventCounts'])
        eventRates        = copy.deepcopy(RnRatePickle['eventRates'])
        eventData = copy.deepcopy(RnRatePickle['eventData'])
        timeMask = (eventData['run_number'] < runChunks[0][0])
        eventData = {key:eventData[key][timeMask] for key in eventData}
    for chunkIdx, eventDataS1amp in enumerate(eventDataS1ampArr):
        if len(eventData) == 0: eventData = copy.deepcopy(eventData_chunk[chunkIdx])
        else:
            for key in eventData: 
                if key in eventData_chunk[chunkIdx]:
                    eventData[key] = np.append(eventData[key], eventData_chunk[chunkIdx][key]) 
                else:
                    eventData[key] =np.append(eventData[key],np.nan*np.empty(len(eventData_chunk[chunkIdx]['s1_corrF'])))
    ###########################################################
    timeDisplay = genTimeDisplay(dateDisplayRange, time_init, time_fin)
    if verbosity >= 1: print('    [time_init, time_fin] =', [time_init, time_fin], ', timeDisplay =', timeDisplay)
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

    correlationFactor = mask_dict['validS1S2Ratio']/mask_dict['XeFidMass']
    rateDict['eventRates_noCorrErr'] = {}
    for peakName in rateDict['eventRates']:
        rateDict['eventRates_noCorrErr'][peakName] = []
        for eventRate in rateDict['eventRates'][peakName]:
            eventRate_noCorrErr = ufloat(-1, 0)
            if eventRate > 0: 
                eventRate_noCorrErr = removeMultiplicationErr(eventRate, correlationFactor)*correlationFactor.n
            rateDict['eventRates_noCorrErr'][peakName].append(eventRate_noCorrErr)
    rateDict['parameters'] = {}
    for key in mask_dict:
        if 'mask' not in key: rateDict['parameters'][key] = mask_dict[key]
    rateDict['parameters']['XeDensity_kg/cm^3'] = mask_dict['XeDensity']
    rateDict['parameters']['XeFidVol_cm^3']     = mask_dict['XeFidVol']
    rateDict['parameters']['XeFidMass_kg']      = mask_dict['XeFidMass']
    with open(output_path, 'wb') as handle:
        pickle.dump(rateDict, handle, protocol=pickle.HIGHEST_PROTOCOL)
    if verbosity >=1: print('    saving:', output_path)

    rateDictLite = ufloatToErrDict(rateDict)
    del rateDictLite['eventData']
    with open(output_path.replace('.pkl', '_lite.pkl'), 'wb') as handle:
        pickle.dump(rateDictLite, handle, protocol=pickle.HIGHEST_PROTOCOL)
    if verbosity >=1: print('    saving:', output_path.replace('.pkl', '_lite.pkl'))
#########################################################################################################################
def correction_derivation(procData, mask_dict, output_dir=None, plot_dir=None, save_plot=None, plot_format=None,\
                          output_name_tag=None, config_dict=None, plot_extra=None, verbosity=None):
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
    if verbosity >= 1: print('  base data:')
    plot_name = 'correction_derivation_0base' 
    mask0 = copy.deepcopy(mask_dict['mask_base'])
    mask1 = copy.deepcopy(mask_dict['mask_final_FV'])
    if save_plot == True:
        plt, ax = plotStandard_s1s2_highA(procData, mask0, mask1, verbosity=verbosity, binN=config_dict['generalBinN'])
        if plot_extra is not None: plot_extra(plot_name, ax, mask_dict) 
        plot_path = os.path.join(plot_dir, dateStr+output_name_tag+plot_name+'.'+plot_format)
        plt.savefig(plot_path, format=plot_format)
        if verbosity >= 1: print('    saving:', plot_path)
        plt.close('all')
    #####
    # First Correction
    #####################################################################################################################
    plot_name = 'correction_derivation_1base_fit1'

    useAFTcorrection = True                            if config_dict is None else config_dict['useAFTcorrection']
    selAFTupperBound = [30.0E4, -0.15, -0.1, 3.8E4]    if config_dict is None else config_dict['selAFTupperBound']
    selAFTlowerBound = [30.0E4, -0.12, -0.1, 5.8E4]    if config_dict is None else config_dict['selAFTlowerBound']
    selDTupperBound  = [1.15E5, 2.58E6, 5.45E5, 3.5E4] if config_dict is None else config_dict['selDTupperBound'] 
    selDTlowerBound  = [1.15E5, 2.58E6, 5.45E5, 5.7E4] if config_dict is None else config_dict['selDTlowerBound']
    if verbosity >= 1: print('  first correction:\n    useAFTcorrection:', useAFTcorrection)
    ###Corr1 Fit Cuts########################################################################
    ### mask_amp_AFT_sel
    mask_gen = copy.deepcopy(mask_dict['mask_final_FV'])
    mask_gen &= (exp2Func(procData['s1_area_fraction_top'], *selAFTupperBound) <= procData['s1_area'])
    mask_gen &= (procData['s1_area'] < exp2Func(procData['s1_area_fraction_top'], *selAFTlowerBound))
    mask_amp_AFT_sel = copy.deepcopy(mask_gen)
    if verbosity >= 1: print('    mask_amp_AFT_sel,', np.sum(mask_amp_AFT_sel), ', done')
    ### mask_amp_DT_sel
    mask_gen = copy.deepcopy(mask_dict['mask_final_FV'])
    mask_gen &= (exp2Func(procData['drift_time'], *selDTupperBound) <= procData['s1_area'])
    mask_gen &= (procData['s1_area'] < exp2Func(procData['drift_time'], *selDTlowerBound))
    mask_amp_DT_sel = copy.deepcopy(mask_gen)
    if verbosity >= 1: print('    mask_amp_DT_sel,', np.sum(mask_amp_DT_sel), ', done')
    ###Corr1 Fit##############################################################################
    maskAFT = copy.deepcopy(mask_amp_AFT_sel)
    dataX = np.array(procData['s1_area'][maskAFT],              dtype='float64')        
    dataY = np.array(procData['s1_area_fraction_top'][maskAFT], dtype='float64')
    costFunc = iminuit.cost.LeastSquares(dataY, dataX, np.ones(len(dataX)), exp2Func)
    objMinuit = iminuit.Minuit(costFunc, *((np.array(selAFTupperBound) + np.array(selAFTlowerBound))/2.0))
    optResult_S1_areaCorr1AFT = objMinuit.migrad()
    S1_areaCorr1Pars_AFT = [val for val in optResult_S1_areaCorr1AFT.values]
    S1_areaCorr1Errs_AFT = [val for val in optResult_S1_areaCorr1AFT.errors]
    if verbosity >= 1: print('    s1_area_fraction_top fit:\n', optResult_S1_areaCorr1AFT)

    maskDT = copy.deepcopy(mask_amp_DT_sel)
    dataX = np.array(procData['s1_area'][maskDT],    dtype='float64')        
    dataY = np.array(procData['drift_time'][maskDT], dtype='float64')
    costFunc = iminuit.cost.LeastSquares(dataY, dataX, np.ones(len(dataX)), exp2Func)
    objMinuit = iminuit.Minuit(costFunc, *((np.array(selDTupperBound) + np.array(selDTlowerBound))/2.0))
    optResult_S1_areaCorr1DT = objMinuit.migrad()
    S1_areaCorr1Pars_DT = [val for val in optResult_S1_areaCorr1DT.values]
    S1_areaCorr1Errs_DT = [val for val in optResult_S1_areaCorr1DT.errors]
    if verbosity >= 1: print('    drift_time fit:\n', optResult_S1_areaCorr1DT)

    procData['s1_corr1_AFT'] = np.array(procData['s1_area'], dtype='float64')
    procData['s1_corr1_AFT']*= S1_areaCorr1Pars_AFT[-1]
    procData['s1_corr1_AFT']/= exp2Func(procData['s1_area_fraction_top'],*S1_areaCorr1Pars_AFT)
    procData['s1_corr1_DT']  = np.array(procData['s1_area'], dtype='float64')
    procData['s1_corr1_DT'] *= S1_areaCorr1Pars_DT[-1]/exp2Func(procData['drift_time'], *S1_areaCorr1Pars_DT)
    mask0 = copy.deepcopy(mask_dict['mask_final'])
    S1_areaCorr1Pars = None
    if useAFTcorrection == True:
        S1_areaCorr1Pars = S1_areaCorr1Pars_AFT
        procData['s1_corr1'] = procData['s1_corr1_AFT']
        mask1 = copy.deepcopy(maskAFT)
    else:
        S1_areaCorr1Pars = S1_areaCorr1Pars_DT
        procData['s1_corr1'] = procData['s1_corr1_DT']
        mask1 = copy.deepcopy(maskDT)
    if save_plot == True:
        plt, ax = plotStandard_s1s2_highA(procData, mask0, mask1, verbosity=verbosity, binN=config_dict['generalBinN'])
        funcY = np.linspace(*ax[3].get_ylim(), 1000)[:-1]
        ax[3].plot(exp2Func(funcY, *selAFTupperBound),     funcY, color='orange', linestyle='dashed', linewidth=1)
        ax[3].plot(exp2Func(funcY, *selAFTlowerBound),     funcY, color='orange', linestyle='dashed', linewidth=1)
        ax[3].plot(exp2Func(funcY, *S1_areaCorr1Pars_AFT), funcY, color='orange', linestyle='solid',  linewidth=2)
        funcY = np.linspace(*ax[6].get_ylim(), 1000)[:-1]
        ax[6].plot(exp2Func(funcY, *selDTupperBound),     funcY, color='orange', linestyle='dashed', linewidth=1)
        ax[6].plot(exp2Func(funcY, *selDTlowerBound),     funcY, color='orange', linestyle='dashed', linewidth=1)
        ax[6].plot(exp2Func(funcY, *S1_areaCorr1Pars_DT), funcY, color='orange', linestyle='solid',  linewidth=2)
        if plot_extra is not None: plot_extra(plot_name, ax, mask_dict) 
        plot_path = os.path.join(plot_dir, dateStr+output_name_tag+plot_name+'.'+plot_format)
        plt.savefig(plot_path, format=plot_format)
        if verbosity >= 1: print('    saving:', plot_path)
        plt.close('all')
    ###########################################################
    plot_name = 'correction_derivation_2corr1'
    mask0 = copy.deepcopy(mask_dict['mask_base'])
    mask1 = copy.deepcopy(mask_dict['mask_final_FV'])
    if save_plot == True:
        plt, ax = plotStandard_s1s2_highA(procData, mask0, mask1, verbosity=verbosity, binN=config_dict['generalBinN'],\
                                          s1ampKey='s1_corr1')
        if plot_extra is not None: plot_extra(plot_name, ax, mask_dict) 
        plot_path = os.path.join(plot_dir, dateStr+output_name_tag+plot_name+'.'+plot_format)
        plt.savefig(plot_path, format=plot_format)
        if verbosity >= 1: print('    saving:', plot_path)
        plt.close('all')
    #####
    # Second Correction
    #####################################################################################################################
    plot_name = 'correction_derivation_3corr2_fit'
    if verbosity >= 1: print('  second correction:') 
    
    selR2upperBound = [4.0E4, -1.4] if config_dict is None else config_dict['selR2upperBound']
    selR2lowerBound = [5.8E4, -1.4] if config_dict is None else config_dict['selR2lowerBound']
    ###Corr2 Fit Cuts########################################################################
    ### mask_amp_r2_sel
    mask_gen = copy.deepcopy(mask_dict['mask_final_FV'])
    mask_gen &= (linearFunc(procData['r2'], *selR2upperBound) <= procData['s1_corr1'])
    mask_gen &= (procData['s1_corr1'] < linearFunc(procData['r2'], *selR2lowerBound))
    mask_amp_r2_sel = copy.deepcopy(mask_gen)
    if verbosity >= 1: print('    mask_amp_r2_sel,', np.sum(mask_amp_r2_sel), ', done')
    #Corr2 Fit###############################################################################
    maskR2 = copy.deepcopy(mask_amp_r2_sel)
    dataX = np.array(procData['s1_corr1'][maskR2], dtype='float64')        
    dataY = np.array(procData['r2'][maskR2],       dtype='float64')
    costFunc = iminuit.cost.LeastSquares(dataY, dataX, np.ones(len(dataX)), linearFunc)
    minuitObj = iminuit.Minuit(costFunc, *((np.array(selR2upperBound) + np.array(selR2upperBound))/2.0)) 
    optResult_S1_areaCorr2 = minuitObj.migrad()
    S1_areaCorr2Pars = [val for val in optResult_S1_areaCorr2.values]
    S1_areaCorr2Errs = [val for val in optResult_S1_areaCorr2.errors]
    if verbosity >= 1: print('    r2 fit:\n', optResult_S1_areaCorr2)

    procData['s1_corr2']  = np.array(procData['s1_corr1'], dtype='float64')
    procData['s1_corr2'] *= S1_areaCorr2Pars[0]/linearFunc(procData['r2'], *S1_areaCorr2Pars)
    mask0 = copy.deepcopy(mask_dict['mask_final'])
    mask1 = copy.deepcopy(maskR2)
    if save_plot == True:
        plt, ax = plotStandard_s1s2_highA(procData, mask0, mask1, verbosity=verbosity, binN=config_dict['generalBinN'],\
                                          s1ampKey='s1_corr1')
        funcY = np.linspace(*ax[8].get_ylim(), 1000)[:-1]
        ax[8].plot(linearFunc(funcY, *selR2upperBound),  funcY, color='orange', linestyle='dashed', linewidth=1)
        ax[8].plot(linearFunc(funcY, *selR2lowerBound),  funcY, color='orange', linestyle='dashed', linewidth=1)
        ax[8].plot(linearFunc(funcY, *S1_areaCorr2Pars), funcY, color='orange', linestyle='solid',  linewidth=1)
        if plot_extra is not None: plot_extra(plot_name, ax, mask_dict) 
        plot_path = os.path.join(plot_dir, dateStr+output_name_tag+plot_name+'.'+plot_format)
        plt.savefig(plot_path, format=plot_format)
        if verbosity >= 1: print('    saving:', plot_path)
        plt.close('all')
    ###########################################################
    plot_name = 'correction_derivation_4corr2'
    mask0 = copy.deepcopy(mask_dict['mask_base'])
    mask1 = copy.deepcopy(mask_dict['mask_final_FV'])
    if save_plot == True:
        plt, ax = plotStandard_s1s2_highA(procData, mask0, mask1, verbosity=verbosity, binN=config_dict['generalBinN'],\
                                          s1ampKey='s1_corr2')
        if plot_extra is not None: plot_extra(plot_name, ax, mask_dict) 
        plot_path = os.path.join(plot_dir, dateStr+output_name_tag+plot_name+'.'+plot_format)
        plt.savefig(plot_path, format=plot_format)
        if verbosity >= 1: print('    saving:', plot_path) 
        plt.close('all')
    #####################################################################################################################
    procData['s1_corrF'] = procData['s1_corr2']
    #####
    # Multi-Peak Fit
    #####################################################################################################################
    plot_name = 'correction_derivation_5peak_fit'
    if verbosity >= 1: print('  multipeak fit:') 
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
    time_diffBinN          = 200              if config_dict is None else config_dict['time_diffBinN']
    time_diffRange_BiPo212 = [0.0E3, 5.0E3]   if config_dict is None else config_dict['time_diffRange_BiPo212']
    time_diffRange_BiPo214 = [0.0E3, 2.0E6]   if config_dict is None else config_dict['time_diffRange_BiPo214']
    BiPo212_peakRange      = [7.0E4, 10.0E4]  if config_dict is None else config_dict['BiPo212_peakRange']
    BiPo212_coinRange      = [6.0E2, 2.7E3]   if config_dict is None else config_dict['BiPo212_coinRange']
    BiPo214_peakRange      = [6.0E4, 7.0E4]   if config_dict is None else config_dict['BiPo214_peakRange']
    BiPo214_coinRange      = [0.7E3, 400.0E3] if config_dict is None else config_dict['BiPo214_coinRange']
    ###BiPo Cuts#############################################################################
    ###  mask_BiPo212_eRange, mask_BiPo212_sel
    mask_gen = copy.deepcopy(mask_dict['mask_base'])
    mask_gen &= ((BiPo212_peakRange[0] <= procData['s1_corrF']) & (procData['s1_corrF'] < BiPo212_peakRange[1]))
    mask_BiPo212_eRange = copy.deepcopy(mask_gen)
    if verbosity >= 1: print('    mask_BiPo212_eRange,', np.sum(mask_BiPo212_eRange), ', done')
    mask_gen &= ((BiPo212_coinRange[0] <= procData['s1_time_diff']) & (procData['s1_time_diff'] < BiPo212_coinRange[1]))
    mask_BiPo212_sel = copy.deepcopy(mask_gen)
    if verbosity >= 1: print('    mask_BiPo212_sel,', np.sum(mask_BiPo212_sel), ', done')
    ###  mask_BiPo214_eRange, mask_BiPo214_sel
    mask_gen = copy.deepcopy(mask_dict['mask_base'])
    mask_gen &= ((BiPo214_peakRange[0] <= procData['s1_corrF']) & (procData['s1_corrF'] < BiPo214_peakRange[1]))
    mask_BiPo214_eRange = copy.deepcopy(mask_gen)
    if verbosity >= 1: print('    mask_BiPo214_eRange,', np.sum(mask_BiPo214_eRange), ', done')
    mask_gen &= ((BiPo214_coinRange[0] <= procData['s1_time_diff']) & (procData['s1_time_diff'] < BiPo214_coinRange[1]))
    mask_BiPo214_sel = copy.deepcopy(mask_gen)
    if verbosity >= 1: print('    mask_BiPo214_sel,', np.sum(mask_BiPo214_sel), ', done')

    mask0 = mask_dict['mask_final_FV']
    mask0_ranged = mask0 & np.logical_not(mask_BiPo212_eRange) & np.logical_not(mask_BiPo214_eRange)
    mask_BiPo212 = mask0 & mask_BiPo212_sel
    mask_BiPo214 = mask0 & mask_BiPo214_sel
    mask1 = mask0_ranged | mask_BiPo212 | mask_BiPo214
    ###Multi-peak Fit#############################################################################
    fitPars, fitErrs, fitCov = None, None, None
    timeDiff_Po212Pars, timeDiff_Po212Errs, timeDiff_Po212Cov = None, None, None
    timeDiff_Po214Pars, timeDiff_Po214Errs, timeDiff_Po214Cov = None, None, None
    
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
        print('    multipeak fit:\n      initPars =', initPars, '\n', optResult)
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
    fitX     = copy.deepcopy(s1_areaX)
    fitY     = peakFitFunc(fitX, fitPars)
    fitYinit = peakFitFunc(fitX, initPars)
    fitResidual = np.array([np.sqrt(iminuit.cost.poisson_chi2(np.array([1.0*histYval]), np.array([1.0*fitYval])))
                            for histYval, fitYval in zip(s1_areaHist, peakFitFunc(s1_areaX, fitPars))])
    ###BiPo Lifetime Fit###################################################################################
    BiPo212_fitRange = BiPo212_coinRange
    timeDiff_Po212 = np.array(procData['s1_time_diff'][mask_BiPo212], dtype='float64')
    timeDiff_Po212_hist, timeDiff_Po212_X = np.histogram(timeDiff_Po212, bins=time_diffBinN, 
                                                         range=time_diffRange_BiPo212)
    timeDiff_Po212_X = timeDiff_Po212_X[:-1]
    timeDiff_Po212_mask = ((BiPo212_fitRange[0] < timeDiff_Po212_X) & (timeDiff_Po212_X < BiPo212_fitRange[1]))
    timeDiff_Po212_X    = timeDiff_Po212_X[timeDiff_Po212_mask]
    timeDiff_Po212_hist = timeDiff_Po212_hist[timeDiff_Po212_mask]
    timeDiff_Po212_err  = np.array([np.sqrt(count) if count != 0 else 1.0 for count in timeDiff_Po212_hist])
    costFunc = iminuit.cost.LeastSquares(timeDiff_Po212_X, timeDiff_Po212_hist, timeDiff_Po212_err, expFunc)
    objMinuit = iminuit.Minuit(costFunc, max(timeDiff_Po212_hist), 300.0, BiPo212_fitRange[0], 0.0)
    objMinuit.fixed['par2'] = True
    objMinuit.fixed['par3'] = True
    optResult = objMinuit.migrad()
    try:
        timeDiff_Po212Pars = [val for val in optResult.values]
        timeDiff_Po212Errs = [err for err in optResult.errors]
        timeDiff_Po212Cov  = [row for row in optResult.covariance]
    except:
        pass
    if verbosity >= 2: print('    timeDiff_Po212 fit:\n', optResult)

    BiPo214_fitRange = BiPo214_coinRange
    timeDiff_Po214 = np.array(procData['s1_time_diff'][mask_BiPo214], dtype='float64')
    timeDiff_Po214_hist, timeDiff_Po214_X = np.histogram(timeDiff_Po214, bins=time_diffBinN, 
                                                         range=time_diffRange_BiPo214)
    timeDiff_Po214_X = timeDiff_Po214_X[:-1]
    timeDiff_Po214_mask = ((BiPo214_fitRange[0] < timeDiff_Po214_X) & (timeDiff_Po214_X < BiPo214_fitRange[1]))
    timeDiff_Po214_X    = timeDiff_Po214_X[timeDiff_Po214_mask]
    timeDiff_Po214_hist = timeDiff_Po214_hist[timeDiff_Po214_mask]
    timeDiff_Po214_err  = np.array([np.sqrt(count) if count != 0 else 1.0 for count in timeDiff_Po214_hist])
    costFunc = iminuit.cost.LeastSquares(timeDiff_Po214_X, timeDiff_Po214_hist, timeDiff_Po214_err, expFunc)
    objMinuit = iminuit.Minuit(costFunc, max(timeDiff_Po214_hist), 164.0E3, BiPo214_fitRange[0], 0.0)
    objMinuit.fixed['par2'] = True
    objMinuit.fixed['par3'] = True
    optResult = objMinuit.migrad()
    try:
        timeDiff_Po214Pars = [val for val in optResult.values]
        timeDiff_Po214Errs = [err for err in optResult.errors]
        timeDiff_Po214Cov  = [row for row in optResult.covariance]
    except:
        pass
    if verbosity >= 2: print('    timeDiff_Po214 fit:\n', optResult)
    ###########################################################
    if save_plot == True:
        gridSpec = [3, 1]   
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
        hist = ax[axIdx].hist(procData['s1_corrF'][mask_BiPo212], bins=s1_ampBinN, range=s1_ampRange,
                              color='blueviolet', linewidth=1)
        hist = ax[axIdx].hist(procData['s1_corrF'][mask_BiPo214], bins=s1_ampBinN, range=s1_ampRange,
                              color='orange', linewidth=1)
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
                ax[axIdx].plot(fitX[fitMask], fitY[fitMask], color='red', linewidth=1.5, alpha=0.8, zorder=3)
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
                fitMask = ((fitRangeSet[peakName][0] <= s1_areaX) & (s1_areaX < fitRangeSet[peakName][1]))
                axRes.plot(s1_areaX[fitMask], fitResidual[fitMask], color='black')
        ######
        axIdx = 1
        hist = ax[axIdx].hist(procData['s1_time_diff'][mask0 & mask_BiPo212_eRange], bins=time_diffBinN, 
                              range=time_diffRange_BiPo212, color='gainsboro', linewidth=1)
        hist = ax[axIdx].hist(procData['s1_time_diff'][mask_BiPo212], bins=time_diffBinN, range=time_diffRange_BiPo212,
                              color='blueviolet', linewidth=1)
        ax[axIdx].set_xlabel('s1_time - alt_s1_time', fontsize=18)
        ax[axIdx].set_ylabel('BiPo212 counts', fontsize=18)
        ax[axIdx].set_xlim(time_diffRange_BiPo212)
        ax[axIdx].set_ylim(bottom=0.5)
        ax[axIdx].set_yscale('log')
        if timeDiff_Po212Cov is not None:
            ax[axIdx].plot(timeDiff_Po212_X, expFunc(timeDiff_Po212_X, *timeDiff_Po212Pars), 
                           color='red', linestyle='solid', linewidth=2)
            fit_life      = ufloat(timeDiff_Po212Pars[1], timeDiff_Po212Errs[1])
            fit_half_life = 0.693*fit_life
            fit_time_zero = ufloat(timeDiff_Po212Pars[2], timeDiff_Po212Errs[2])
            timeDiff_Po212_efficiency  =  uncertainties.umath.exp(-(BiPo212_fitRange[1]-fit_time_zero)/fit_life)
            timeDiff_Po212_efficiency -=  uncertainties.umath.exp(-(BiPo212_fitRange[0]-fit_time_zero)/fit_life)
            timeDiff_Po212_efficiency /= -uncertainties.umath.exp(-(0                  -fit_time_zero)/fit_life)
            boxText  = 'Po212 fit half life = ' + scientificStr(fit_half_life.n)
            boxText +=                  '$\pm$' + scientificStr(fit_half_life.s, 2) + 'ns\n'
            boxText += 'cut efficiency      = ' + scientificStr(100.0*timeDiff_Po212_efficiency.n)
            boxText +=                  '$\pm$' + scientificStr(100.0*timeDiff_Po212_efficiency.s, 2) + '%'
            ax[axIdx].text(0.6, 0.95, boxText, transform=ax[axIdx].transAxes, fontsize=18,\
                           horizontalalignment='left', verticalalignment='top',\
                           bbox=dict(boxstyle='square', facecolor='white', alpha=0.5))
        ######
        axIdx = 2
        hist = ax[axIdx].hist(procData['s1_time_diff'][mask0 & mask_BiPo214_eRange], bins=time_diffBinN, 
                              range=time_diffRange_BiPo214, color='gainsboro', linewidth=1)
        hist = ax[axIdx].hist(procData['s1_time_diff'][mask_BiPo214], bins=time_diffBinN, range=time_diffRange_BiPo214,
                              color='orange', linewidth=1)
        ax[axIdx].set_xlabel('s1_time - alt_s1_time', fontsize=18)
        ax[axIdx].set_ylabel('BiPo214 counts', fontsize=18)
        ax[axIdx].set_xlim(time_diffRange_BiPo214)
        ax[axIdx].set_ylim(bottom=0.5)
        ax[axIdx].set_yscale('log')
        if timeDiff_Po214Cov is not None:
            ax[axIdx].plot(timeDiff_Po214_X, expFunc(timeDiff_Po214_X, *timeDiff_Po214Pars), 
                           color='red', linestyle='solid', linewidth=2)
            fit_life      = ufloat(timeDiff_Po214Pars[1], timeDiff_Po214Errs[1])
            fit_half_life = 0.693*fit_life
            fit_time_zero = ufloat(timeDiff_Po214Pars[2], timeDiff_Po214Errs[2])
            timeDiff_Po214_efficiency  =  uncertainties.umath.exp(-(BiPo214_fitRange[1]-fit_time_zero)/fit_life)
            timeDiff_Po214_efficiency -=  uncertainties.umath.exp(-(BiPo214_fitRange[0]-fit_time_zero)/fit_life)
            timeDiff_Po214_efficiency /= -uncertainties.umath.exp(-(0                  -fit_time_zero)/fit_life)
            boxText  = 'Po214 fit half life = ' + scientificStr(fit_half_life.n/1000.0) 
            boxText +=                  '$\pm$' + scientificStr(fit_half_life.s/1000.0, 2) + '$\mu$s\n'
            boxText += 'cut efficiency      = ' + scientificStr(100.0*timeDiff_Po214_efficiency.n)
            boxText +=                  '$\pm$' + scientificStr(100.0*timeDiff_Po214_efficiency.s, 2) + '%'
            ax[axIdx].text(0.6, 0.95, boxText, transform=ax[axIdx].transAxes, fontsize=18,\
                           horizontalalignment='left', verticalalignment='top',\
                           bbox=dict(boxstyle='square', facecolor='white', alpha=0.5))
        if plot_extra is not None: plot_extra(plot_name, ax, mask_dict) 
        plot_path = os.path.join(plot_dir, dateStr+output_name_tag+plot_name+'.'+plot_format)
        plt.savefig(plot_path, format=plot_format)
        if verbosity >= 1: print('    saving:', plot_path)
        plt.close('all')
    #####################################################################################################################
    corrParDict = {}
    corrParDict['BiPo212_peakRange'] = BiPo212_peakRange
    corrParDict['BiPo212_coinRange'] = BiPo212_coinRange
    corrParDict['BiPo214_peakRange'] = BiPo214_peakRange
    corrParDict['BiPo214_coinRange'] = BiPo214_coinRange
    corrParDict['S1_areaCorr1Pars_AFT'] = S1_areaCorr1Pars_AFT
    corrParDict['S1_areaCorr1Errs_AFT'] = S1_areaCorr1Errs_AFT
    corrParDict['S1_areaCorr1Pars_DT']  = S1_areaCorr1Pars_DT
    corrParDict['S1_areaCorr1Errs_DT']  = S1_areaCorr1Errs_DT
    corrParDict['S1_areaCorr2Pars']     = S1_areaCorr2Pars
    corrParDict['S1_areaCorr2Errs']     = S1_areaCorr2Errs
    corrParDict['fitRangeSet'], corrParDict['peak_normPar'], corrParDict['peak_normErr'] = {}, {}, {}
    for peakName in fitRangeSet:
        corrParDict['fitRangeSet'][peakName]  = fitRangeSet[peakName]
        corrParDict['peak_normPar'][peakName] = peak_normPar[peakName]
        corrParDict['peak_normErr'][peakName] = peak_normErr[peakName]
    corrParDict['cali_factor'] = cali_factor
    corrParDict['parameters'] = {}
    for key in mask_dict:
        if 'mask' not in key: corrParDict['parameters'][key] = mask_dict[key]
    output_path = os.path.join(output_dir, output_name_tag+'correction_derivation.pkl')
    with open(output_path, 'wb') as handle:
        pickle.dump(corrParDict, handle, protocol=pickle.HIGHEST_PROTOCOL)
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
    activity *= 1.0E6/(XeFidMass*dataTime)
    return activity
def peak_region_rate_estimation(runInfoData, procData, efficiency_dict, mask_dict, s1corr_paths, output_dir=None,\
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
    deadtime_type     = None
    generalBinN       = 400
    countRangeSigN    = 3 
    chunkCheckTags    = False
    combTimeThres     = 6*3600.0
    lightYieldObsPar  = 'cs1'
    chargeYieldObsPar = 'cs2'
    s1_ampBinN        = 800
    s1_ampRange       = [1.0E4, 1.0E5]
    chi2DisplayRange  = [0.0,    'inf']
    dateDisplayRange  = ['-inf', 'inf']
    peakDisplaySubset = ['Rn222']
    if config_dict is not None:
        if 'light_yield_ref_path' in config_dict: lightYieldRefPath = config_dict['light_yield_ref_path']
        if 'deadtime_type'     in config_dict: deadtime_type     = config_dict['deadtime_type']
        if 'generalBinN'       in config_dict: generalBinN       = config_dict['generalBinN']
        if 'countRangeSigN'    in config_dict: countRangeSigN    = config_dict['countRangeSigN']
        if 'chunkCheckTags'    in config_dict: chunkCheckTags    = config_dict['chunkCheckTags']
        if 'combTimeThres'     in config_dict: combTimeThres     = config_dict['combTimeThres']
        if 'lightYieldObsPar'  in config_dict: lightYieldObsPar  = config_dict['lightYieldObsPar']
        if 'chargeYieldObsPar' in config_dict: chargeYieldObsPar = config_dict['chargeYieldObsPar']
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
                 'lightYieldObsPar':  lightYieldObsPar,
                 'chargeYieldObsPar': chargeYieldObsPar,
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
    
    if 'ave_efficiency' not in efficiency_dict:
        efficiency_dict['ave_efficiency'] = ufloat(1.0, 0.0)
    if 'run_number' not in efficiency_dict:
        efficiency_dict['run_number']     = np.array(runInfoData['number'])
        efficiency_dict['run_efficiency'] = [efficiency_dict['ave_efficiency'] \
                                             for run_number in efficiency_dict['run_number']]
    runsNotAvailInAlphaEff = []
    if config_dict is not None:
        efficiency_full = []
        for run_number in runInfoData['number']:
            if run_number in efficiency_dict['run_number']:
                eff_run_sel = np.where((efficiency_dict['run_number'] == int(run_number)))[0]
                if len(eff_run_sel) > 1:
                    raise AssertionError('peak_region_rate_estimation(): run '+str(run_number)+\
                                         ' is not unique in efficiency_dict[\'run_number\']')
                efficiency_full.append(efficiency_dict['run_efficiency'][eff_run_sel[0]])
            elif (run_number < efficiency_dict['run_number'][0]) or \
                 (efficiency_dict['run_number'][-1] < run_number):
                efficiency_full.append(efficiency_dict['ave_efficiency'])
            else:
                runsNotAvailInAlphaEff.append(run_number)
        efficiency_dict['run_number']     = np.array(runInfoData['number'])
        efficiency_dict['run_efficiency'] = efficiency_full
    if len(runsNotAvailInAlphaEff) > 0:
        raise AssertionError('peak_region_rate_estimation(): inbetween runs are missing in efficiency_dict:\n'+\
                             str(runsNotAvailInAlphaEff))

    fitLightYieldCorrRanges   = None        # for scaling to find the best fit ranges, default=fitLightYieldCorrPeaks
    fitLightYieldCorrPeaks    = ['Rn222', 'Po218', 'Rn220', 'Po216', 'Po214']
    lightYieldCorrCountRanges = []          # for count in range not included in multipeak fit
    lightChargeYieldObsPeaks  = ['Rn222']
    fitLightYieldCorrCountThres = 20
    fitInitRangeFromConfig      = False
    fitInitRangeSet             = {}
    if config_dict is not None:
        if 'fitLightYieldCorrRanges' in config_dict: fitLightYieldCorrRanges = config_dict['fitLightYieldCorrRanges']
        fitLightYieldCorrPeaks = config_dict['fitLightYieldCorrPeaks']
        if 'lightYieldCorrCountRanges' in config_dict: lightYieldCorrCountRanges=config_dict['lightYieldCorrCountRanges']
        if 'lightChargeYieldObsPeaks'  in config_dict: lightChargeYieldObsPeaks = config_dict['lightChargeYieldObsPeaks']
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
    mask0 = copy.deepcopy(mask_dict['mask_base'])
    mask1 = copy.deepcopy(mask_dict['mask_final_FV'])
    #mask1 = mask0 & mask_dict['mask_z_cut'] & mask_dict['mask_r2_cut'] & np.logical_not(mask_dict['mask_final_FV'])
    #mask1 = mask1 & (2.0E4 <= procData['s2_area']) & (procData['s2_area'] < 3.5E4)
    #mask1 = mask1 & (3.8E4 <= procData['s1_corrF']) & (procData['s1_corrF'] < 3.95E4)
    #####################################################################################################################
    if verbosity >= 1: print('  base data:')
    plot_name = 'peak_region_rate_estimation_0base' 
    if save_plot == True:
        plt, ax = plotStandard_s1s2_highA(procData, mask0, mask1, verbosity=verbosity, binN=generalBinN)
        funcY = np.linspace(*ax[3].get_ylim(), 1000)[:-1]
        ax[3].plot(exp2Func(funcY, *corrPar['S1_areaCorr1Pars_AFT']), funcY,color='orange',linestyle='solid',linewidth=2)
        funcY = np.linspace(*ax[6].get_ylim(), 1000)[:-1]
        ax[6].plot(exp2Func(funcY, *corrPar['S1_areaCorr1Pars_DT']),  funcY,color='orange',linestyle='solid',linewidth=2)
        if plot_extra is not None: plot_extra(plot_name, ax, mask_dict) 
        plot_path = os.path.join(plot_dir, dateStr+output_name_tag+plot_name+'.'+plot_format)
        plt.savefig(plot_path, format=plot_format)
        if verbosity >= 1: print('    saving:', plot_path)
        plt.close('all')
    ###S1 Corr###############################################################################
    plot_name = 'peak_region_rate_estimation_1s1corr'
    if save_plot == True:
        plt, ax = plotStandard_s1s2_highA(procData, mask0,mask1,verbosity=verbosity,binN=generalBinN,s1ampKey='s1_corrF')
        if plot_extra is not None: plot_extra(plot_name, ax, mask_dict) 
        plot_path = os.path.join(plot_dir, dateStr+output_name_tag+plot_name+'.'+plot_format)
        plt.savefig(plot_path, format=plot_format)
        if verbosity >= 1: print('    saving:', plot_path)
        plt.close('all')
    #####
    # Chunk Multi-Peak Rate with No Light Yeild Correction
    #####################################################################################################################
    plot_name = 'peak_region_rate_estimation_2peak_fit0'
    if verbosity >= 1: print('  multipeak fit, no LY corr:')
    if len(procData['run_number'][mask0]) == 0: 
        raise AssertionError('peak_region_rate_estimation(): mask0 cuts out all events')
    time_init = toBeginningOfDay(1.0E-9*procData['time'][mask0][0])
    time_fin  = 1.0E-9*procData['endtime'][mask0][-1]

    peak_amps, peak_width_percent = {}, {}
    for peakIdx, peakName in enumerate(peak_normParPickle):
        peak_amps[peakName]          = peak_normParPickle[peakName][0]
        peak_width_percent[peakName] = 100.0*peak_normParPickle[peakName][1]/peak_normParPickle[peakName][0]
    if verbosity >= 1: 
        print('    XeDensity   =', mask_dict['XeDensity'], 'kg/cm^3')
        print('    XeFidMass   =', mask_dict['XeFidMass'], 'kg')    

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
        for runChunk_ref_path in reference_runChunk_paths:
            if os.path.isfile(runChunk_ref_path) == False:
                raise AssertionError('peak_region_rate_estimation(): pickle file:\n  '+runChunk_ref_path+\
                                     '\ndoes not exist')
            runChunks_, runModes_ = genRunChunksWithRef(np.unique(procData['run_number']), runChunk_ref_path,\
                                                        verbosity=(verbosity-1))
            runChunks.extend(runChunks_)
            runModes .extend(runModes_)
    else:
        runChunks, runModes = genRunChunks(np.unique(procData['run_number']), runInfoData, combTimeThres,\
                                           checkTags=chunkCheckTags, verbosity=(verbosity-1))
    ###########################################################
    timestampRanges, liveTimes, dayArr = [], [], []
    eventDataS1ampArr, eventDataS1ampArr_noBiPo, eventDataS1ampArr_multipeakFit, eventData_chunk = [], [], [], []
    eventCounts = {peakName:[] for peakName in peak_normParPickle}
    eventRates  = {peakName:[] for peakName in peak_normParPickle}
    eventAlphaEff = []
    if (deadtime_type is not None) and ((deadtime_type in runInfoData) == False): deadtime_type == None
    for runChunk in tqdm(runChunks):
        timestampRangeChunk, liveTimeChunk = [], ufloat(0, 0)
        alphaEff = 0.0
        mask_run_sel = (procData['run_number'] == -1)
        for run_number in runChunk:      
            timestampRangeChunk.append([runInfoData['start_timestamp'][runInfoData['name']==run_number][0],
                                        runInfoData['end_timestamp']  [runInfoData['name']==run_number][0]])
            liveTime = pd.Timedelta(runInfoData['livetime'][runInfoData['name']==run_number][0]).total_seconds()
            if deadtime_type is not None: liveTime-=runInfoData[deadtime_type][runInfoData['name']==run_number][0]*1.0E-9
            else: warnings.warn('peak_region_rate_estimation(): deadtime plugin NOT implemented', Warning) 
            liveTimeChunk += ufloat(liveTime, 0.5)          #NOTE: 0.5 gps uncertainty for livetime
            mask_run_sel |= (procData['run_number'] == run_number)
            eff_run_sel = np.where((efficiency_dict['run_number'] == int(run_number)))[0]
            if len(eff_run_sel) > 1: 
                raise AssertionError('peak_region_rate_estimation(): run '+run_number+\
                                     ' in not unique in the updated efficiency_dict[\'run_number\']')
            alphaEff += liveTime*efficiency_dict['run_efficiency'][eff_run_sel[0]]
        timestampRanges.append(timestampRangeChunk)
        liveTimes      .append(liveTimeChunk)
        dayArr         .append((timestampRangeChunk[0][0] - time_init)/(3600.0*24.0))
        eventAlphaEff  .append(alphaEff/liveTimeChunk) 

        eventDataS1amp  = np.array(procData['s1_corrF'][mask1 & mask_run_sel], dtype='float64')
        eventDataS1ampArr .append(eventDataS1amp)
        eventData_chunk.append({key:procData[key][mask1 & mask_run_sel] for key in procData})
        if 'mask_BiPoRemoval' in mask_dict:
            eventDataS1amp_noBiPo = np.array(procData['s1_corrF'][mask1 & mask_run_sel & mask_dict['mask_BiPoRemoval']],\
                                             dtype='float64')
            eventDataS1ampArr_noBiPo.append(eventDataS1amp_noBiPo) 
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
        activityPars = [countDataSet, countRangeSet, peak_normParPickle, liveTimeChunk, mask_dict['XeFidMass']]
        for peakName in peak_normParPickle: eventRates[peakName].append(get_activity(peakName, *activityPars))
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
            hist = ax[axIdx].hist(procData['s1_corrF'][mask_isotope], bins=s1_ampBinN, range=s1_ampRange,\
                                  color='gainsboro', linewidth=1, alpha=0.3)

        linePlot, linePlotText= {}, {}
        for peakName in peak_normParPickle:
            linePlot[peakName] = ax[axIdx].axvline(x=peak_amps[peakName], ymin=0.0, ymax=ax[axIdx].get_ylim()[1],\
                                                   color=ALPHA_LINE_C[peakName], linestyle='dashed', linewidth=1.5)
            linePlotText[peakName] = peakName + ', w=' + scientificStr(peak_width_percent[peakName]) + '%'
        for peakName in ALPHA_LINE_E:
            if peakName not in peak_normParPickle:
                peak_amp = cali_factor*ALPHA_LINE_E[peakName]
                linePlot[peakName] = ax[axIdx].axvline(x=peak_amp, ymin=0.0, ymax=ax[axIdx].get_ylim()[1],\
                                                       color=ALPHA_LINE_C[peakName], linestyle='dashed', linewidth=1.5)
                linePlotText[peakName] = peakName + ''
        legObj = ax[axIdx].legend([linePlot[peakName] for peakName in ALPHA_LINE_E],\
                                  [linePlotText[peakName] for peakName in ALPHA_LINE_E], loc='upper left', fontsize=16)
        ######
        axIdx = 1
        tracePlot, tracePlotText = {}, {}
        for peakName in peak_normParPickle:
            tracePlot[peakName] = ax[axIdx].errorbar(dayArr, [item.n for item in eventRates[peakName]],\
                                                             [item.s for item in eventRates[peakName]],\
                                                     color=ALPHA_LINE_C[peakName], fmt='s', **ERROR_PLOT_PARS)[0]
            tracePlotText[peakName]  = peakName
        setXaxisDateLabel(ax[axIdx], time_init, [time_init, time_fin])
        ax[axIdx].set_ylabel('activity concentration (uBq/kg_Xe)', fontsize=18)
        ax[axIdx].set_ylim(bottom=3.0E-2)
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
    plot_dir_name   = 'peak_region_rate_estimation_3chunk_fit'
    plot_name_chunk = 'peak_region_rate_estimation_3chunk_multipeak_fit'

    peakFitFunc = multiSigLinkedGausPDF_scipy
    output_path = os.path.join(output_dir, output_name_tag+'peak_region_rate_estimation.pkl')
    ###########################################################
    fitLightYieldCorrPeaks = [fitLightYieldCorrPeak for fitLightYieldCorrPeak in fitLightYieldCorrPeaks
                              if fitLightYieldCorrPeak in fitRangeSetPickle]
    if fitLightYieldCorrRanges is None: fitLightYieldCorrRanges = fitLightYieldCorrPeaks
    lightYieldCorrCountRanges = {peakName:lightYieldCorrCountRanges[peakName] for peakName in lightYieldCorrCountRanges\
                                 if peakName not in fitLightYieldCorrPeaks}
    lightChargeYieldObsPeaks = [peakName for peakName in lightChargeYieldObsPeaks if peakName in fitLightYieldCorrPeaks]
    peakDisplaySubset        = [peakName for peakName in peakDisplaySubset if peakName in fitLightYieldCorrPeaks]
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
            fitRanges_gainCorr[peakName] = copy.deepcopy(RnRatePickle['fitRanges_gainCorr'][peakName])
        optFitPars            = copy.deepcopy(RnRatePickle['optFitPars'])
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
    if len(eventDataS1ampArr) == 0: raise AssertionError('peak_region_rate_estimation(): no run chunks to process')
    for chunkIdx, eventDataS1amp in enumerate(eventDataS1ampArr):
        eventDataS1ampC = np.array(eventDataS1amp, dtype='float64')
        eventDataS1ampFit = copy.deepcopy(eventDataS1ampC)
        if 'mask_BiPoRemoval' in mask_dict:
            eventDataS1ampFit = np.array(eventDataS1ampArr_noBiPo[chunkIdx], dtype='float64')
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
                optCount = copy.deepcopy(countInRange)
                optIdx   = copy.deepcopy(segIdx)
            if verbosity >= 3: print('   ', segIdx, segIdx*scaleStep, countInRange, optCount)
        optIdx -= int(optIdx_sameN/2)
        for peakName in fitLightYieldCorrPeaks:
            fitRangeChunk[peakName] = np.array(fitRangeSetPickle[peakName])*(1.0 + optIdx*scaleStep)
            if peakName not in fitRanges_gainCorr: fitRanges_gainCorr[peakName] = []
            fitRanges_gainCorr[peakName].append(fitRangeChunk[peakName])
        ##################################################################
        if verbosity >= 1: 
            print(chunkIdx, '/', len(eventDataS1ampArr), runChunks[chunkIdx])
            print('  day:       ', dayArr[chunkIdx])
            print('  tot eventN:', len(eventDataS1ampC))
            print('  livetime:  ', liveTimes[chunkIdx])
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
            if verbosity >= 1: print('    multipeak fit:\n', optResult)
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
                fitResidual = np.array([np.sqrt(iminuit.cost.poisson_chi2(np.array([1.0*histYval]),
                                                                          np.array([1.0*fitYval])))
                                        for histYval, fitYval in zip(s1ampHist, peakFitFunc(s1ampX, fitPars))])
            else:
                peak_normPar_chunk = {}
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
            figSize, marginRatio = getSizeMargin(gridSpec, subplotSize=[18.0, 5.0], plotMargin=[0.06, 0.14, 0.07, 0.12])
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
            ax[axIdx].set_title(titleStr, fontsize=16, y=1.01)
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
                    fitMask = ((fitRangeChunk[peakName][0] < s1ampX) & (s1ampX < fitRangeChunk[peakName][1]))
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
    plot_name = 'peak_region_rate_estimation_4peak_rate'
    if verbosity >= 1: print('  rate evalution:')    
    ###########################################################
    countRanges_gainCorr = {}
    eventCounts_gainCorr = {peakName:[] for peakName in fitLightYieldCorrPeaks}
    eventRates_gainCorr  = {peakName:[] for peakName in fitLightYieldCorrPeaks}
    eventRates_gainCorr_effCorr  = {peakName:[] for peakName in fitLightYieldCorrPeaks}
    light_yields, charge_yields = [], []
    mask_isotope = {}
    for peakName in fitLightYieldCorrPeaks: mask_isotope[peakName] = copy.deepcopy(mask1)
    ###########################################################
    if RnRatePickle is not None:
        for peakName in RnRatePickle['countRanges_gainCorr']:
            countRanges_gainCorr[peakName] = copy.deepcopy(RnRatePickle['countRanges_gainCorr'][peakName])
        for peakName in RnRatePickle['eventCounts_gainCorr']:
            eventCounts_gainCorr[peakName]        = copy.deepcopy(RnRatePickle['eventCounts_gainCorr'][peakName]) 
            eventRates_gainCorr[peakName]         = copy.deepcopy(RnRatePickle['eventRates_gainCorr'][peakName])
            eventRates_gainCorr_effCorr[peakName] = copy.deepcopy(RnRatePickle['eventRates_gainCorr_effCorr'][peakName])
        light_yields  = copy.deepcopy(RnRatePickle['light_yields'])
        charge_yields = copy.deepcopy(RnRatePickle['charge_yields'])
    chunkIdxOffset = len(light_yields)
    for chunkIdx, eventDataS1amp in enumerate(eventDataS1ampArr):  
        if verbosity >= 1: print(chunkIdx, runChunks[chunkIdx])
        mask_run_sel = np.array([(run_number in runChunks[chunkIdx]) for run_number in procData['run_number']])
        eventDataS1ampC = np.array(eventDataS1amp, dtype='float64')
        cs1data = procData[lightYieldObsPar] [mask1 & mask_run_sel] 
        cs2data = procData[chargeYieldObsPar][mask1 & mask_run_sel]
        peak_sel_masks = {}
        cs1data_peaks, cs2data_peaks, cs1data_vals, cs2data_vals = {}, {}, [], []
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
                    eventCounts_gainCorr       [peakName] = []
                    eventRates_gainCorr        [peakName] = []
                    eventRates_gainCorr_effCorr[peakName] = []
                if peakName not in mask_isotope: mask_isotope[peakName] = copy.deepcopy(mask1) 
            for peakName in countRangeSet:
                if peakName not in countRanges_gainCorr: countRanges_gainCorr[peakName] = []
                countRanges_gainCorr[peakName].append(np.array(countRangeSet[peakName]))
                peakData = eventDataS1ampC[(countRangeSet[peakName][0] <= eventDataS1ampC)&
                                           (eventDataS1ampC < countRangeSet[peakName][1])]
                eventCounts_gainCorr[peakName].append(len(peakData))
            countDataSet = {peakName:eventCounts_gainCorr[peakName][-1] for peakName in eventCounts_gainCorr}
            activityPars = [countDataSet, countRangeSet, peakMuSigSet, liveTimes[chunkIdx], mask_dict['XeFidMass']]
            for peakName in countRangeSet: 
                eventRates_gainCorr        [peakName].append(get_activity(peakName, *activityPars))
                eventRates_gainCorr_effCorr[peakName].append(get_activity(peakName, *activityPars)/\
                                                             eventAlphaEff[chunkIdx])   # efficiency correction
                peak_sel_masks[peakName] = mask_run_sel & (countRangeSet[peakName][0] <= procData['s1_corrF']) &\
                                                          (procData['s1_corrF'] < countRangeSet[peakName][1]) 
                mask_isotope[peakName] &= (np.logical_not(mask_run_sel) | peak_sel_masks[peakName])
            if verbosity >= 1:
                print('    liveTime =', liveTimes[chunkIdx], 's')
                print('    XeFidMass =', mask_dict['XeFidMass'], 'kg')
                for peakName in peakMuSigSet:
                    print('    peakMuSigSet['+peakName+'] =', peakMuSigSet[peakName])
                for peakName in countRanges_gainCorr:
                    print('    countRanges_gainCorr['+peakName+'] =',        countRanges_gainCorr[peakName][-1])
                for peakName in eventCounts_gainCorr:
                    print('    eventCounts_gainCorr['+peakName+'] =',        eventCounts_gainCorr[peakName][-1])
                for peakName in eventRates_gainCorr:
                    print('    eventRates_gainCorr['+peakName+'] =',         eventRates_gainCorr[peakName][-1])
                for peakName in eventRates_gainCorr_effCorr:
                    print('    eventRates_gainCorr_effCorr['+peakName+'] =', eventRates_gainCorr_effCorr[peakName][-1])
            ###
            for peakName in lightChargeYieldObsPeaks:
                cs1data_peaks[peakName] = procData[lightYieldObsPar] [mask1 & peak_sel_masks[peakName]]
                cs2data_peaks[peakName] = procData[chargeYieldObsPar][mask1 & peak_sel_masks[peakName]] 
                light_yield_val   = np.median(cs1data_peaks[peakName])
                light_yield_err   = np.std(cs1data_peaks[peakName],ddof=1)/np.sqrt(max(1.0,len(cs1data_peaks[peakName])))
                light_yield_err  *= MEDIAN_SE_CONST
                charge_yield_val  = np.median(cs2data_peaks[peakName])
                charge_yield_err  = np.std(cs2data_peaks[peakName],ddof=1)/np.sqrt(max(1.0,len(cs2data_peaks[peakName])))
                charge_yield_err *= MEDIAN_SE_CONST
                cs1data_vals.append(ufloat(light_yield_val,  light_yield_err) /(ALPHA_LINE_E[peakName]*1.0E3)) #PE/keV
                cs2data_vals.append(ufloat(charge_yield_val, charge_yield_err)/(ALPHA_LINE_E[peakName]*1.0E3)) #PE/keV
            light_yields .append(weightedAverage(cs1data_vals))
            charge_yields.append(weightedAverage(cs2data_vals))
            if verbosity >= 1:
                print('    light_yield  =', light_yields[-1],  'PE/keV')
                print('    charge_yield =', charge_yields[-1], 'PE/keV')
        else:
            for peakName in countRanges_gainCorr:
                countRanges_gainCorr[peakName].append(-1)
            for peakName in fitLightYieldCorrPeaks:
                eventCounts_gainCorr[peakName]       .append(0)
                eventRates_gainCorr [peakName]       .append(ufloat(-1.0, 0.0))
                eventRates_gainCorr_effCorr[peakName].append(ufloat(-1.0, 0.0))
            for peakName in lightYieldCorrCountRanges:
                if peakName not in eventCounts_gainCorr:
                    eventCounts_gainCorr       [peakName] = []
                    eventRates_gainCorr        [peakName] = []
                    eventRates_gainCorr_effCorr[peakName] = []
                    mask_isotope               [peakName] = copy.deepcopy(mask1)
                eventCounts_gainCorr       [peakName].append(0)
                eventRates_gainCorr        [peakName].append(ufloat(-1.0, 0.0))
                eventRates_gainCorr_effCorr[peakName].append(ufloat(-1.0, 0.0))
                if peakName not in mask_isotope: mask_isotope[peakName] = copy.deepcopy(mask1)
                mask_isotope[peakName] &= np.logical_not(mask_run_sel)
            light_yields .append(ufloat(-1, 0))
            charge_yields.append(ufloat(-1, 0))
        ###########################################################
        if (save_plot == True) and (len(lightChargeYieldObsPeaks) != 0):
            gridSpec = [len(lightChargeYieldObsPeaks), 1]
            figSize, marginRatio = getSizeMargin(gridSpec, subplotSize=[18.0, 5.0], plotMargin=[0.06, 0.14, 0.03, 0.08])
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
            for axIdx, peakName in enumerate(lightChargeYieldObsPeaks):
                rangeXY = [s1_ampRange, [0.0, 6.0E4]]
                hist = ax[axIdx].hist2d(cs1data, cs2data, bins=(s1_ampBinN, int(s1_ampBinN/4)), cmap=CMAP, alpha=0.2,\
                                        cmin=1.0, range=[rangeXY[0], rangeXY[1]])
                if len(cs1data_peaks) != 0:
                    hist = ax[axIdx].hist2d(cs1data_peaks[peakName], cs2data_peaks[peakName],\
                                            bins=(s1_ampBinN, int(s1_ampBinN/4)),\
                                            cmap=CMAP, cmin=1.0, range=[rangeXY[0], rangeXY[1]])
                divider = make_axes_locatable(ax[axIdx])
                cax = divider.append_axes('right', size='5%', pad=0.1)
                cbar = fig.colorbar(hist[3], ax=ax[axIdx], cax=cax)
                runChunkTitle = [str(runChunk[runIdx:runIdx          + 12])
                                 for runIdx in range(0, len(runChunk), 12)]
                dateInit = timestampToDate(timestampRanges[chunkIdx][0][0]).split('\n')
                dateInit = dateInit[0] + '_' + dateInit[1]
                dateFin  = timestampToDate(timestampRanges[chunkIdx][-1][1]).split('\n')
                dateFin  = dateFin[0]  + '_' + dateFin[1]
                titleStr  = peakName + ': ' + str(eventCounts_gainCorr[peakName][-1]) + ' counts => '
                titleStr += str(eventRates_gainCorr[peakName][-1]) + ' uBq/kg'
                if axIdx == 0:
                    titleStr += '\nLY: ' + str(light_yields[-1])  + 'PE/keV, '
                    titleStr +=   'CY: ' + str(charge_yields[-1]) + 'PE/keV\n'
                    ax[axIdx].set_title(titleStr, fontsize=18, y=0.85)
                else:
                    ax[axIdx].set_title(titleStr, fontsize=18)
                ax[axIdx].set_xlabel(lightYieldObsPar,  fontsize=18)
                ax[axIdx].set_ylabel(chargeYieldObsPar, fontsize=18)
                ax[axIdx].set_xlim(*rangeXY[0])
                ax[axIdx].set_ylim(*rangeXY[1])
            if plot_extra is not None: plot_extra(plot_name_chunk, ax, mask_dict)
            plot_subdir = dateStr + output_name_tag + plot_dir_name
            pathlib.Path(os.path.join(plot_dir, plot_subdir)).mkdir(parents=True, exist_ok=True)
            plot_path = os.path.join(plot_dir, plot_subdir, dateStr+plot_name_chunk+\
                                     str(runChunks[chunkIdx][0])+'___.'+plot_format)
            plt.savefig(plot_path, format=plot_format)
            fuseFigVer(plot_path.replace('___', ''), plot_path)
            if verbosity >= 1: print('    saving:', plot_path.replace('___', ''))
            plt.close('all')
    if RnRatePickle is not None:
        time_init = copy.deepcopy(RnRatePickle['time_init'])
        runChunks       = copy.deepcopy(RnRatePickle['runChunks'])
        runModes        = copy.deepcopy(RnRatePickle['runModes'])
        timestampRanges = copy.deepcopy(RnRatePickle['timestampRanges'])
        liveTimes       = copy.deepcopy(RnRatePickle['liveTimes'])
        dayArr          = copy.deepcopy(RnRatePickle['dayArr'])
        eventCounts = copy.deepcopy(RnRatePickle['eventCounts'])
        eventRates  = copy.deepcopy(RnRatePickle['eventRates'])
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
        gridSpec = [5, 1]
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
        for peakName in eventRates_gainCorr_effCorr:
            tracePlot[peakName] = \
              ax[axIdx].errorbar(np.array(dayArr)[validMask],\
                                 [item.n for item in np.array(eventRates_gainCorr_effCorr[peakName])[validMask]],\
                                 [item.s for item in np.array(eventRates_gainCorr_effCorr[peakName])[validMask]],\
                                 color=ALPHA_LINE_C[peakName], fmt='s', **ERROR_PLOT_PARS)[0]
            tracePlotText[peakName] = peakName
        setXaxisDateLabel(ax[axIdx], time_init, timeDisplay, displayOff=True)
        ax[axIdx].set_ylabel('activity concentration\n(uBq/kg_Xe)', fontsize=18)
        ax[axIdx].set_yscale('log')
        ax[axIdx].secondary_yaxis('right')
        if 'Po218' in tracePlotText: tracePlotText['Po218'] = tracePlotText['Po218'].replace('Po218', 'Po218/Bi212')
        legObj = ax[axIdx].legend([tracePlot[peakName]     for peakName in eventRates_gainCorr_effCorr\
                                   if (peakName in ALPHA_LINE_E) and (peakName not in peakDisplaySubset)],\
                                  [tracePlotText[peakName] for peakName in eventRates_gainCorr_effCorr\
                                   if (peakName in ALPHA_LINE_E) and (peakName not in peakDisplaySubset)],\
                                  loc='upper left', fontsize=16, markerscale=1.5)
        ######
        axIdx = 1
        subTracePlot, subTracePlotText = {}, {}
        for peakName in peakDisplaySubset:
            subTracePlot[peakName] = \
              ax[axIdx].errorbar(np.array(dayArr)[validMask],\
                                 [item.n for item in np.array(eventRates_gainCorr_effCorr[peakName])[validMask]],
                                 [item.s for item in np.array(eventRates_gainCorr_effCorr[peakName])[validMask]],
                                 color=ALPHA_LINE_C[peakName], fmt='s', **ERROR_PLOT_PARS)[0]
            subTracePlotText[peakName] = peakName
        setXaxisDateLabel(ax[axIdx], time_init, timeDisplay, displayOff=True)
        ax[axIdx].set_ylabel('activity concentration\n(uBq/kg_Xe)', fontsize=18)
        ax[axIdx].set_yscale('log')
        ax[axIdx].secondary_yaxis('right')
        legObj =ax[axIdx].legend([subTracePlot[peakName] for peakName in peakDisplaySubset if peakName in ALPHA_LINE_E],\
                                 [subTracePlotText[peakName]\
                                  for peakName in peakDisplaySubset if peakName in ALPHA_LINE_E],\
                                 loc='upper right', fontsize=14, markerscale=1.5)
        ######
        axIdx = 2
        ax[axIdx].errorbar(np.array(dayArr)[validMask],\
                           [item.n for item in np.array(eventAlphaEff)[validMask]],
                           [item.s for item in np.array(eventAlphaEff)[validMask]],
                           color='blue', fmt='s', **ERROR_PLOT_PARS)[0]
        setXaxisDateLabel(ax[axIdx], time_init, timeDisplay)
        ax[axIdx].set_ylabel('efficiency', fontsize=18)
        ######
        axIdx = 3
        ax[axIdx].errorbar(np.array(dayArr)[validMask],\
                           [item.n for item in np.array(light_yields)[validMask]], 
                           [item.s for item in np.array(light_yields)[validMask]], 
                           color='blue', fmt='s', **ERROR_PLOT_PARS)[0]
        setXaxisDateLabel(ax[axIdx], time_init, timeDisplay, displayOff=True)
        ax[axIdx].set_ylabel('light yield (PE/keV)', fontsize=18, color='tab:blue', alpha=0.7)
        axRight = ax[axIdx].twinx()
        
        axRight.errorbar(np.array(dayArr)[validMask],\
                         [item.n for item in np.array(charge_yields)[validMask]],
                         [item.s for item in np.array(charge_yields)[validMask]],
                         color='green', fmt='s', **ERROR_PLOT_PARS)[0]
        axRight.set_ylabel('charge yield (PE/keV)', fontsize=18, color='tab:green')
        ax[axIdx].set_zorder(axRight.get_zorder()+1)
        ax[axIdx].set_frame_on(False)
        ######
        axIdx = 4
        ax[axIdx].scatter(np.array(dayArr)[validMask], np.array(chi2_ndofs)[validMask], color='blue')
        setXaxisDateLabel(ax[axIdx], time_init, timeDisplay)
        ax[axIdx].set_ylabel('multipeak fit chi2/ndof', fontsize=18)
        ax[axIdx].set_ylim(bottom=0.0)
        ax[axIdx].secondary_yaxis('right')
        ax[axIdx].axhline(y=1.0, color='grey', linestyle='dashed', zorder=0)
        ######
        if plot_extra is not None: plot_extra(plot_name, ax, mask_dict, extra_ax=[axRight], time_init=time_init) 
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
    rateDict['eventRates_gainCorr_noCorrErr'] = {}
    rateDict['eventRates_gainCorr_effCorr']  = eventRates_gainCorr_effCorr
    rateDict['eventRates_gainCorr_effCorr_noCorrErr'] = {}    
 
    correlationFactor = 1.0/mask_dict['XeFidMass']
    for peakName in rateDict['eventRates_gainCorr']:
        rateDict['eventRates_gainCorr_noCorrErr'][peakName] = []
        for eventRate in rateDict['eventRates_gainCorr'][peakName]:
            eventRate_noCorrErr = ufloat(-1, 0)
            if eventRate > 0: 
                eventRate_noCorrErr = removeMultiplicationErr(eventRate, correlationFactor)*correlationFactor.n
            rateDict['eventRates_gainCorr_noCorrErr'][peakName].append(eventRate_noCorrErr)
        rateDict['eventRates_gainCorr_effCorr_noCorrErr'][peakName] = []
        for eventRate in rateDict['eventRates_gainCorr_effCorr'][peakName]:
            eventRate_noCorrErr = ufloat(-1, 0)
            if eventRate > 0: 
                eventRate_noCorrErr = removeMultiplicationErr(eventRate, correlationFactor)*correlationFactor.n
            rateDict['eventRates_gainCorr_effCorr_noCorrErr'][peakName].append(eventRate_noCorrErr)
    rateDict['light_yields']  = light_yields
    rateDict['charge_yields'] = charge_yields
    rateDict['correction_derivation_parameters'] = corrParPickle
    rateDict['parameters'] = {}
    for key in inputPars: rateDict['parameters'][key] = inputPars[key]
    for key in mask_dict:
        if ('mask' not in key) and ('Xe' not in key): rateDict['parameters'][key] = mask_dict[key]
    rateDict['parameters']['XeDensity_kg/cm^3'] = mask_dict['XeDensity']
    rateDict['parameters']['XeFidVol_cm^3']     = mask_dict['XeFidVol']
    rateDict['parameters']['XeFidMass_kg']      = mask_dict['XeFidMass']
    with open(output_path, 'wb') as handle:
        pickle.dump(rateDict, handle, protocol=pickle.HIGHEST_PROTOCOL)
    if verbosity >=1: print('    saving:', output_path)
    
    rateDictLite = ufloatToErrDict(rateDict)
    del rateDictLite['eventData']
    with open(output_path.replace('.pkl', '_lite.pkl'), 'wb') as handle:
        pickle.dump(rateDictLite, handle, protocol=pickle.HIGHEST_PROTOCOL)
    if verbosity >=1: print('    saving:', output_path.replace('.pkl', '_lite.pkl'))
#########################################################################################################################
class RadonAnalysisDataBuilder(AnalysisRunBatch):
    """ Manage building of the data necessary for analysis. """
    def __init__(self, config_dict, run_list, load_func, **kwargs):
        super().__init__(**kwargs)
        radon_analysis_class_initialization(self, config_dict, run_list, load_func, **kwargs)
        self.runs_per_job = self.config['data_building']['runs_per_job']
        if (self.config['data_building']['save_run_pickle'] == True) or\
           (self.config['data_building']['save_event_pickle'] == True):
            self.output_dir = self.config['data_building']['pickle_dir_name']
            pathlib.Path(self.output_dir).mkdir(parents=True, exist_ok=True)
        if self.analysis_mode == 'run':
            self.load_context_with_run_info(self)
            if self.config['general']['load_deadtime'] == True: 
                load_deadtime_to_run_info(self, self.config['general']['deadtime_type'])
            if self.config['data_building']['save_run_pickle'] == True:
                output_name = self.output_dir + '/__' + self.targets[0] + '_runInfo.pkl'
                self.run_info.to_pickle(output_name)
                if self.verbosity >= 1: print('  saving:', output_name)
    def get_runlist(self):
        if self.verbosity >= 1: print('  loading run list:')
        run_list = None
        if self.run_list is not None: run_list = self.run_list
        elif self.config['data_building']['load_from_run_list'] == True:
            run_list = read_run_list_file(self.config['data_building']['run_list_file'])
            if self.verbosity >= 1: print('Generating run list from', self.config['data_building']['run_list_file'])
        else:
            if self.context is None: self.load_context_with_run_info(self)
            run_info_dict = self.run_info.to_dict('series')
            run_list = run_info_dict['number'].tolist()
        run_list = list(set(run_list))
        run_list.sort()
        if self.verbosity >= 1: print(run_list)
        return run_list
    def submit_single(self, runlist, **kwargs):
        runlist_str = ' '.join(runlist)
        job_process_str = f'python {self.file} run --runlist {runlist_str}'
        if '000' not in kwargs['jobname']: 
            job_process_str += f' --config_edit \"data_building:save_run_pickle=false\"'
            job_process_str += f' --config_edit \"general:load_deadtime=false\"'
        self.logging.debug(f'Call \"_submit\" with job: {job_process_str} and kwargs: {kwargs}.')
        self._submit(job_process_str, **kwargs)
    def run_single(self, run):
        data_df_args = copy.deepcopy(self.config['general']['data_df'])
        del data_df_args['targets']
        del data_df_args['selection_str']
        if self.verbosity >= 1: 
            print('Processing run '+str(run))
            print('  targets:         ', self.targets)
            print('  selection_str:   ', self.selection_str)
            print('  other get_df args', data_df_args)
        availTarget = {}
        run_id = '0' + str(run)
        if self.config['data_building']['save_event_pickle'] == True:
            output_name = self.output_dir + '/' + self.targets[0] + '_' + run_id + '.pkl'
            if os.path.exists(output_name) is False:
                if self.context.is_stored(run_id, target=self.run_info_available_target) == True:
                    data_df = None
                    if 'peaks' in self.targets:
                        self.context.make(run_id, 'peaks', save='peaks')
                        self.context.storage = self.context.storage[:-1]
                    try:
                        data_df = self.context.get_df(run_id, targets=self.targets, selection_str=self.selection_str,\
                                                      **{key:data_df_args[key] for key in data_df_args})
                    except Exception as errMsg:
                        ### NOTE: checking which cuts do not exist
                        for target in self.targets:
                            try:
                                _ = self.context.get_df(run_id, targets=target)
                                if self.verbosity >= 1: print('  true:', target)
                                availTarget[target] = True
                            except:
                                if self.verbosity >= 1: print('\n  false:', target)
                                availTarget[target] = False
                        print('  run_id: '+run_id+', availTarget: ', availTarget)
                        ### cut_mv_veto specification
                        only_cut_mv_veto_false = True
                        for target in self.targets:
                            if target != 'cut_mv_veto':
                                only_cut_mv_veto_false &= availTarget[target]
                        if only_cut_mv_veto_false == True:
                            try:
                                data_df = self.context.get_df(run_id, targets=[target for target in self.targets \
                                                                               if (target != 'cut_mv_veto')],\
                                                              selection_str=self.selection_str,\
                                                              **{key:data_df_args[key] for key in data_df_args})
                                data_df['cut_mv_veto'] = True
                            except Exception as errMsg:
                                print('\nDATA LOADING FAILS:\n', errMsg, '\n')
                        else:
                            print('\nDATA LOADING FAILS:\n', errMsg, '\n')
                    if data_df is not None:
                        data_df.to_pickle(output_name)
                        del data_df
                        if self.verbosity >= 1: print('  saving:', output_name)
                elif self.verbosity >= 1: print('  run '+run_id+': data not found for target='+str(self.targets))
            elif self.verbosity >= 1: print('  run '+run_id+' pickle file already exist')
            if self.verbosity >= 1: print('')
#########################################################################################################################
# helper functions
def radon_analysis_class_initialization(self, config_dict, run_list, load_func, **kwargs):
    """ Execute initialization code common to both radon analysis classes. """
    self.config = copy.deepcopy(config_dict)
    self.verbosity     = self.config['general']['verbosity']
    self.save_plot     = self.config['general']['save_plot']
    self.plot_format   = self.config['general']['plot_format']
    self.analysis_mode = self.config['general']['analysis_mode']
    self.run_list = copy.deepcopy(run_list)
    self.load_context_with_run_info = load_func
    self.file = self.file.replace('_s1s2.',
        str(self.load_context_with_run_info).split('load_context_with_run_info')[1].split(' ')[0]+'.')
    self.context = None
    if self.verbosity >= 1: print('Environment:', os.environ['CONDA_DEFAULT_ENV'])
def load_context_with_run_info_s1s2(self):
    import strax, straxen, cutax
    if self.verbosity >= 1: straxen.print_versions(('strax', 'straxen', 'cutax', 'xedocs', 'rframe'))
    #####################################################################################################################
    @strax.takes_config(strax.Option('s1_area_hard_high_thres', default=2.0E4,\
                        help='To reveal only high energy events'),)
    class cut_s1_area_hard_high_thres(strax.CutPlugin):
        depends_on = ('event_info',)
        provides = 'cut_s1_area_hard_high_thres'
        cut_name = 'cut_s1_area_hard_high_thres'
        cut_description = 'To reveal only high energy events'
        __version__ = '0.0.1'
        def cut_by(self, events):
            greater_than_thres = events['s1_area'] > self.config['s1_area_hard_high_thres']
            return greater_than_thres
    self.register_lists = [cutax.cut_lists.basic, cutax.cut_lists.science]                          #plugin registers
    try:
        self.register_lists.append(cutax.plugins.deadtime)
    except:
        pass
    self.register_singles    = [cut_s1_area_hard_high_thres]
    self.targets_fixed       = ["cut_s1_area_hard_high_thres"]
    self.selection_str_fixed = ["cut_s1_area_hard_high_thres"]
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
        raise AssertionError('load_context_with_run_info_s1s2(): no runs have been found in self.run_info, '+\
                             'check runInfoDF_selection function from:\n  '+self.config['general']['python_config_file'])
    self.run_info.reset_index(drop=True, inplace=True)
    if self.verbosity >= 1: 
        print('  available modes:', set(self.run_info['mode']), '\n  leftover tags:', set(self.run_info['tags']))

    self.targets = copy.deepcopy(self.config['general']['data_df']['targets'])
    self.targets.extend([target for target in self.targets_fixed if target not in self.targets])
    self.selection_str = copy.deepcopy(self.config['general']['data_df']['selection_str'])
    self.selection_str.extend([selection for selection in self.selection_str_fixed\
                               if selection not in self.selection_str])
def load_deadtime_to_run_info(self, deadtimeType):
    if self.verbosity >= 2: print('  loading deadtime_info:')
    self.run_info = self.run_info.sort_values('number', ascending=False)
    self.run_info.reset_index(drop=True, inplace=True)

    run_info_dict = self.run_info.to_dict('series')
    run_name_list = run_info_dict['name'].tolist()
    if len(run_name_list) == 0:
        raise AssertionError('load_deadtime_to_run_info(): no runs have been found in self.run_info')
    if self.verbosity >= 2: 
        run_name_list_display = copy.deepcopy(run_name_list)
        run_name_list_display.sort()
        print('  available run names from context:\n   ', run_name_list_display)
    
    deadtime_info_list, deadtime_NaN_list, deadtimeType_NaN_list = [], [], []
    for run_id in run_name_list:
        deadDF = None
        try:
            deadDF = self.context.get_df(run_id, 'deadtime_tpc', progress_bar=(self.verbosity >= 2))
            deadDF.rename(columns={'deadtime':      'deadtime_tpc'},           inplace=True)
            deadDF.rename(columns={'lifetime_loss': 'lifetime_loss_tpc_only'}, inplace=True)
        except:
            warnings.warn('load_deadtime_to_run_info(): deadtime_tpc not available for run_id: '+run_id, Warning)
            deadtime_NaN_list.append(run_id)
        if deadDF is not None:  
            if deadtimeType not in deadDF.columns:
                deadDF_ = None
                try:
                    deadDF_ = self.context.get_df(run_id, 'deadtime_tpc_mv', progress_bar=(self.verbosity >= 2))
                except:
                    warnings.warn('load_deadtime_to_run_info(): deadtime_tpc_mv not available for run_id: '+run_id,
                                  Warning)
                if deadDF_ is not None:
                    deadDF_.rename(columns={'deadtime': 'deadtime_tpc_mv'}, inplace=True)
                    deadDF_.drop('lifetime_loss', axis=1, inplace=True)
                    deadDF = pd.concat([deadDF, deadDF_.loc[:, ~deadDF_.columns.isin(deadDF.columns)]], axis=1)
            if deadtimeType not in deadDF.columns:
                deadDF_ = None
                try:
                    deadDF_ = self.context.get_df(run_id, 'deadtime', progress_bar=(self.verbosity >= 2))       
                except:
                    warnings.warn('load_deadtime_to_run_info(): full deadtime not available for run_id: '+run_id,Warning)
                if deadDF_ is not None:
                    deadDF = pd.concat([deadDF, deadDF_.loc[:, ~deadDF_.columns.isin(deadDF.columns)]], axis=1)
                else:
                    if 'deadtime_tpc_mv' in deadDF: 
                        deadDF['deadtime']      = deadDF['deadtime_tpc_mv']
                        deadDF['lifetime_loss'] = deadDF['lifetime_loss_tpc_mv']
                    else:
                        deadDF['deadtime']      = deadDF['deadtime_tpc']
                        deadDF['lifetime_loss'] = deadDF['lifetime_loss_tpc_only']
            if deadtimeType not in deadDF.columns:
                deadtimeType_NaN_list.append(run_id) 
                raise AssertionError('load_deadtime_to_run_info(): deadtimeType '+deadtimeType+' not found '+\
                                     'for run_id: '+str(run_id))

            deadDF['run_id'] = [run_id]
            deadtime_info_list.append(deadDF)
    if len(deadtime_NaN_list) != 0:
        raise AssertionError('load_deadtime_to_run_info(): deadtime not available for run_id: '+run_id)
    if len(deadtimeType_NaN_list) != 0:
        raise AssertionError('load_deadtime_to_run_info(): deadtimeType '+deadtimeType+' not found '+\
                             'for run_ids:\n'+str(deadtimeType_NaN_list))

    deadtime_info = pd.concat(deadtime_info_list)
    if deadtimeType not in deadtime_info.columns:
        raise AssertionError('load_deadtime_to_run_info(): deadtimeType not in any run_info loaded: '+str(deadtimeType)+\
                             '\n  amoung the list: '+str(deadtime_info.columns))
    #if 'lifetime_loss' in deadtimeType:
    #    deadtime_info[deadtimeType] = deadtime_info[deadtimeType].fillna(deadtime_info['lifetime_loss'])
    #elif 'deadtime' in deadtimeType:
    #    deadtime_info[deadtimeType] = deadtime_info[deadtimeType].fillna(deadtime_info['deadtime'])
    deadtime_info = deadtime_info.sort_values('run_id', ascending=False)
    deadtime_info.reset_index(drop=True, inplace=True)
    self.run_info = pd.concat([self.run_info, deadtime_info], axis=1)
    if self.verbosity >= 2: print(self.run_info.head(20), '\n', self.run_info.tail(20))
def load_run_event_data(self, iter_range=None, timestamp_range=None, run_exclude=None, run_include=None,\
                        runInfoDF_selection=None, run_chunks=None, verbosityLevel=None):
    self.run_info = self.run_info.sort_values('number', ascending=True)
    self.run_info.reset_index(drop=True, inplace=True)
    start_list, end_list = None, None 
    if str(self.run_info['start'].dtypes) == 'datetime64[ns]':
        start_list = self.run_info['start'].values.astype(np.int64)/1.0E9
        end_list   = self.run_info['end']  .values.astype(np.int64)/1.0E9
    elif str(self.run_info['start'].dtypes) == 'object':
        start_list = [pdTimestamp.timestamp() for pdTimestamp in self.run_info['start'].values]
        end_list   = [pdTimestamp.timestamp() for pdTimestamp in self.run_info['end'].values]
    self.run_info['start_timestamp'] = start_list
    self.run_info['end_timestamp']   = end_list
    run_info_temp = copy.deepcopy(self.run_info)

    if self.run_list is not None:
        run_list_names = run_numbers_to_run_names(self.run_list)
        run_info_temp = run_info_temp[run_info_temp.name.apply(lambda name: name in run_list_names)]
    run_info_temp = run_info_temp[timestamp_range[0] < run_info_temp['start_timestamp']]
    run_info_temp = run_info_temp[run_info_temp['end_timestamp'] < timestamp_range[1]]            
    if hasattr(self.python_module, 'runInfoDF_selection'):
        run_info_temp = self.python_module.runInfoDF_selection(run_info_temp)
    else:
        print('  Note: runInfoDF_selection() does not exist in general:python_config_file')
    run_info_sel = copy.deepcopy(run_info_temp['name'].tolist())
    run_info_sel = [run_number_str for run_number_str in run_info_sel if run_number_str not in run_exclude]
    if run_include is not None:
        run_info_sel = [run_number_str for run_number_str in run_info_sel if run_number_str in run_include]
    run_info_sel = run_info_sel[iter_range[0]:iter_range[1]]
    if len(run_info_sel) == 0:
        raise AssertionError('load_run_event_data(): no runs selected, check '+\
                             'iter_range(=\"updateNew\"?)/timestamp_range/run_exclude/run_include/'+\
                             'runInfoDF_selection from:\n  '+\
                             self.config['general']['config_file']+'\n  '+self.config['general']['python_config_file'])
    elif verbosityLevel >= 10: print('  loading the following runs:\n   ', run_info_sel) 
    first_run_number_str = run_info_sel[0] 

    extra_runs_to_include = []
    if (run_chunks is not None) and (len(run_chunks) > 0):
        chunk_mask = np.array([(int(first_run_number_str) <= int(runChunk[-1])) for runChunk in run_chunks])
        if verbosityLevel >= 20:
            print('  handling existing run chunks:\n    first_run_number_str =', first_run_number_str)
            for chunk_idx in range(len(run_chunks)): 
                print('    reprocessing:', chunk_mask[chunk_idx], run_chunks[chunk_idx])
        if np.sum(chunk_mask) > 0: extra_runs_to_include = copy.deepcopy(run_chunks[chunk_mask][0])
        else:                      extra_runs_to_include = copy.deepcopy(run_chunks[-1])
    run_info_sel = np.append(run_info_sel, extra_runs_to_include)
    run_info_sel = np.sort(np.unique(run_info_sel))
           
    self.run_info = self.run_info[self.run_info.name.apply(lambda name: name in run_info_sel)]
    if self.config['general']['load_deadtime'] == True: load_deadtime_to_run_info(self)     #loading deadtime here
    runInfoDF = self.run_info
    runInfoData = runInfoDF.to_dict('series')
    for key in runInfoData: runInfoData[key] = np.array(runInfoData[key]) 

    data_df_args = copy.deepcopy(self.config['general']['data_df'])
    del data_df_args['targets']
    del data_df_args['selection_str']
    if verbosityLevel >= 1:
        print('Processing runs: ')
        print('  targets:         ', self.targets)
        print('  selection_str:   ', self.selection_str)
        print('  other get_df args', data_df_args) 
    
    if 'peaks' in self.targets:
        for run_id in run_info_sel:
            self.context.make(run_id, 'peaks', save='peaks')
            self.context.storage = self.context.storage[:-1]
    procDF = self.context.get_df(run_info_sel, targets=self.targets, selection_str=self.selection_str,\
                                 **{key:data_df_args[key] for key in data_df_args})
    if len(run_info_sel) == 1: procDF['run_id'] = [run_info_sel[0] for eventIter in range(procDF.shape[0])]
    procData = procDF.to_dict('series')
    for key in procData: procData[key] = np.array(procData[key])
    try:
        procData['run_number'] = np.array([run_id.decode("utf-8") for run_id in procData['run_id']])
    except:
        procData['run_number'] = np.array(procData['run_id']).astype(str)
    procDF['run_number'] = procData['run_number']

    if verbosityLevel >= 0:
        print('  number of runs:  ', len(runInfoData['name']))
        print('  number of events:', len(procData[list(procData.keys())[0]]))
        if verbosityLevel >= 1: print('  loading run numbers:', runInfoData['number'])
        display(runInfoDF.head(verbosityLevel))
        display(runInfoDF.tail(verbosityLevel))
        display(procDF.head(verbosityLevel))
        display(procDF.tail(verbosityLevel))
    return runInfoDF, runInfoData, procDF, procData
def trim_dict_entries(input_dict, rm_key_contain='//'):
    output_dict = copy.deepcopy(input_dict)
    for key, item in copy.deepcopy(input_dict).items():
        if isinstance(output_dict[key], dict) == True: 
            output_dict[key] = trim_dict_entries(output_dict[key], rm_key_contain=rm_key_contain)
        elif rm_key_contain in key: del output_dict[key]
    return output_dict
def replace_dict_string_string(input_dict, replace_from, replace_to):
    output_dict = copy.deepcopy(input_dict)
    for key, item in copy.deepcopy(input_dict).items():
        if isinstance(output_dict[key], dict) == True:
            output_dict[key] = replace_dict_string_string(output_dict[key], replace_from, replace_to)
        elif isinstance(output_dict[key], list) == True:
            for listIdx, listItem in enumerate(output_dict[key]):
                if isinstance(listItem, str) == True:
                    output_dict[key][listIdx] = listItem.replace(replace_from, replace_to)
        elif isinstance(output_dict[key], str) == True:
            output_dict[key] = output_dict[key].replace(replace_from, replace_to)
    return output_dict
def replace_dict_string_value(input_dict, replace_from, replace_to):
    output_dict = copy.deepcopy(input_dict)
    for key, item in copy.deepcopy(input_dict).items():
        if isinstance(output_dict[key], dict) == True:
            output_dict[key] = replace_dict_string_value(output_dict[key], replace_from, replace_to)
        elif isinstance(output_dict[key], list) == True:
            for listIdx, listItem in enumerate(output_dict[key]):
                if isinstance(listItem, str) == True:
                    if output_dict[key][listIdx] == replace_from: output_dict[key][listIdx] = replace_to
        elif isinstance(output_dict[key], str) == True:
            if output_dict[key] == replace_from: output_dict[key] = replace_to
    return output_dict
def genTimeDisplay(dateDisplayRange, time_init, time_fin):
    timeDisplay = [None, None]
    if isinstance(dateDisplayRange[0], str) == False:
        if 10000 <= dateDisplayRange[0]: dateDisplayRange[0] = str(dateDisplayRange[0])
    if isinstance(dateDisplayRange[1], str) == False:
        if 10000 <= dateDisplayRange[1]: dateDisplayRange[1] = str(dateDisplayRange[1])
 
    if isinstance(dateDisplayRange[0], str) == False:
        timeDisplay = [None, YYMMDDtoTimestamp(dateDisplayRange[1])]
        if np.isinf(timeDisplay[1]) == True: timeDisplay[1] = copy.deepcopy(toBeginningOfDay(time_fin)+24*3600.0)
        timeDisplay[0] = timeDisplay[1] + dateDisplayRange[0]*24*3600.0
        if np.isinf(timeDisplay[0]) == True: timeDisplay[0] = copy.deepcopy(toBeginningOfDay(time_init))
    elif isinstance(dateDisplayRange[1], str) == False:
        timeDisplay = [YYMMDDtoTimestamp(dateDisplayRange[0]), None]
        if np.isinf(timeDisplay[0]) == True: timeDisplay[0] = copy.deepcopy(toBeginningOfDay(time_init))
        timeDisplay[1] = timeDisplay[0] + dateDisplayRange[1]*24*3600.0
        if np.isinf(timeDisplay[1]) == True: timeDisplay[1] = copy.deepcopy(toBeginningOfDay(time_fin)+24*3600.0)
    else:
        timeDisplay = [YYMMDDtoTimestamp(dateDisplayRange[0]), YYMMDDtoTimestamp(dateDisplayRange[1])]
        if np.isinf(timeDisplay[0]) == True: timeDisplay[0] = copy.deepcopy(toBeginningOfDay(time_init))
        if np.isinf(timeDisplay[1]) == True: timeDisplay[1] = copy.deepcopy(toBeginningOfDay(time_fin)+24*3600.0)
    if timeDisplay[0] >= timeDisplay[1]:
        raise AssertionError('genTimeDisplay(): dateDisplayRange[1] < dateDisplayRange[0]')
    return timeDisplay
#########################################################################################################################
def main():
    args        = ARGS
    config_dict = CONFIG_DICT
    config_dict['general']['config_file']   = copy.deepcopy(args.config_file)
    config_dict['general']['analysis_mode'] = copy.deepcopy(args.analysis_mode)
    if 'load_deadtime' not in config_dict['general']: config_dict['general']['load_deadtime'] = True
    if 'deadtime_type' not in config_dict['general']: config_dict['general']['deadtime_type'] = 'lifetime_loss_tpc_only'
    if 'generalBinN'   not in config_dict['general']: config_dict['general']['generalBinN']   = 400
    config_dict = trim_dict_entries(config_dict)
    for editStr in args.config_edit:
        keyList = editStr.split(':')
        values = keyList[-1].split('=')[1].split(',')
        output_value = []
        for value_ in values:
            value = copy.deepcopy(value_)
            valueCheck = value.replace('-', '').replace('.', '').replace('e', '').replace('E', '')
            if valueCheck.isdigit() == True:
                if ('.' in value) or ('e' in value.lower()): value = float(value)
                else:                                        value = int(value)
            else:
                value = value.replace(' ', '').replace('\'', '')
                if   value.lower() == 'true':  value = True
                elif value.lower() == 'false': value = False
            output_value.append(value)
        if len(output_value) == 1: output_value = copy.deepcopy(output_value[0])
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
    if verbosity >= 1: 
        print('Reading config_file:', args.config_file, '\nStarting radon analysis...')
        print('  input arguments:', args)  
    
    ###
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
                                                load_func=load_context_with_run_info_s1s2,\
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
            submit_args['jobname'] += (date.today().strftime('%y%m%d') + '_data_builder_s1s2')
            submit_args['jobname'] += config_dict['analysis_main']['output_name_tag']
            if verbosity >= 1: print('  submission parameters:', submit_args)
            data_builder.submit(**{key:submit_args[key] for key in submit_args})
        elif analysis_mode == 'run': data_builder.run( Namespace(**{'runlist':data_builder.get_runlist()}) )
    if len(analysis_main_steps_to_process) > 0:
        analysis_main = S1S2RadonAnalysisMain(config_dict=config_dict, run_list=args.runlist,\
                                              load_func=load_context_with_run_info_s1s2,\
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
            submit_args['jobname'] += (date.today().strftime('%y%m%d') + '_analysis_main_s1s2_')
            submit_args['jobname'] += config_dict['analysis_main']['output_name_tag']
            if verbosity >= 1: print('  submission parameters:', submit_args)
            analysis_main.submit(**{key:submit_args[key] for key in submit_args})
        elif analysis_mode == 'run': analysis_main.run(run_args)
    if verbosity >= 1: print('End of code')
#########################################################################################################################
if __name__ == '__main__': main()









