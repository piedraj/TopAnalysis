from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.transferLogs = True
config.General.requestName  = 'TTJets'

config.section_("JobType")
config.JobType.pluginName  = 'Analysis'
config.JobType.psetName    = 'skimToTreeSUSYMCtfs.py'
config.JobType.outputFiles = ['Tree_13TeV.root']
config.JobType.inputFiles  = ['Winter14_V5_DATA_L2L3Residual_AK5PFchs.txt','Winter14_V5_DATA_L3Absolute_AK5PFchs.txt','Winter14_V5_DATA_L2Relative_AK5PFchs.txt','Winter14_V5_DATA_L1FastJet_AK5PFchs.txt','Summer13_V5_DATA_UncertaintySources_AK5PFchs.txt','PHYS14_V2_MC_L3Absolute_AK4PFchs.txt','PHYS14_V2_MC_L2Relative_AK4PFchs.txt','PHYS14_V2_MC_L1FastJet_AK4PFchs.txt']

config.section_("Data")
config.Data.inputDataset    = '/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM'
config.Data.splitting       = 'FileBased'
config.Data.unitsPerJob     = 1
config.Data.publishDataName = 'NoSkim'
config.Data.ignoreLocality  = False

config.section_("Site")
config.Site.storageSite = 'T2_ES_IFCA'

config.section_("User")
