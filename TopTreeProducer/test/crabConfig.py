from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.transferLogs = False
config.General.requestName  = 'TTbar_Powheg'

config.section_("JobType")
config.JobType.pluginName  = 'Analysis'
config.JobType.psetName    = 'skimToTreeSUSYMCtfs_JEC.py'
config.JobType.outputFiles = ['Tree.root','Histos.root']
config.JobType.inputFiles  = ['PY8_RunIISpring15DR74_bx50_MC.db',
#'Winter14_V5_DATA_L2L3Residual_AK5PFchs.txt','Winter14_V5_DATA_L3Absolute_AK5PFchs.txt','Winter14_V5_DATA_L2Relative_AK5PFchs.txt','Winter14_V5_DATA_L1FastJet_AK5PFchs.txt',
'Summer13_V5_DATA_UncertaintySources_AK5PFchs.txt']
#'PHYS14_V4_MC_L3Absolute_AK4PFchs.txt','PHYS14_V4_MC_L2Relative_AK4PFchs.txt','PHYS14_V4_MC_L1FastJet_AK4PFchs.txt']

config.section_("Data")
config.Data.inputDataset    = '/TT_TuneCUETP8M1_13TeV-powheg-pythia8/RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v4/MINIAODSIM'

#config.Data.splitting       = 'LumiBased'
#config.Data.unitsPerJob     = 50
#config.Data.publishDataName = 'NoSkim'
config.Data.splitting       = 'FileBased'
config.Data.unitsPerJob     = 1
config.Data.publishDataName = 'Skim_2L_Pt17_8'
config.Data.publication     = False
config.Data.outLFNDirBase = '/store/user/jfernan/DR74X/'
config.Data.ignoreLocality  = True

config.section_("Site")
config.Site.storageSite = 'T2_ES_IFCA'

config.section_("User")
