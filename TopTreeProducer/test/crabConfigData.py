from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.transferLogs = False
config.General.requestName   = 'DataMu'
config.General.requestName   = 'DataEle'
config.General.requestName   = 'DataDoubleMu'
#config.General.requestName   = 'DataDoubleEG'
#config.General.requestName   = 'DataMuEG'

config.section_("JobType")
config.JobType.pluginName  = 'Analysis'
config.JobType.psetName    = 'skimToTreeSUSYDatatfs_JEC.py'
config.JobType.outputFiles = ['Tree.root']
config.JobType.inputFiles  = ['PY8_RunIISpring15DR74_bx50_MC.db',
#'Winter14_V5_DATA_L2L3Residual_AK5PFchs.txt','Winter14_V5_DATA_L3Absolute_AK5PFchs.txt','Winter14_V5_DATA_L2Relative_AK5PFchs.txt','Winter14_V5_DATA_L1FastJet_AK5PFchs.txt',
'Summer13_V5_DATA_UncertaintySources_AK5PFchs.txt']
#'PHYS14_V4_MC_L3Absolute_AK4PFchs.txt','PHYS14_V4_MC_L2Relative_AK4PFchs.txt','PHYS14_V4_MC_L1FastJet_AK4PFchs.txt']

config.section_("Data")

config.Data.inputDataset    = '/SingleMuon/Run2015B-PromptReco-v1/MINIAOD'
config.Data.inputDataset    = '/SingleElectron/Run2015B-PromptReco-v1/MINIAOD'
config.Data.inputDataset    = '/DoubleMuon/Run2015B-PromptReco-v1/MINIAOD'
#config.Data.inputDataset    = '/DoubleEG/Run2015B-PromptReco-v1/MINIAOD'
#config.Data.inputDataset    = '/MuonEG/Run2015B-PromptReco-v1/MINIAOD'

config.Data.lumiMask = '/afs/cern.ch/user/j/jfernan2/work2/public/miniAOD/CMSSW_7_4_2_patch1/src/DR74X.50ns/json_DCSONLY_Run2015B.txt'
config.Data.splitting       = 'LumiBased'
config.Data.unitsPerJob     = 100
#config.Data.publishDataName = 'NoSkim'
config.Data.publishDataName = 'Skim_2L_Pt17_8_Puppy'
config.Data.publication     = False
config.Data.outLFNDirBase = '/store/user/jfernan/DR74X/'
config.Data.ignoreLocality  = True

config.section_("Site")
config.Site.storageSite = 'T2_ES_IFCA'

config.section_("User")
