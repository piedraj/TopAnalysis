from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
#config.General.requestName = 'DYJetsToLL_PU20bx25_PHYS14'
#config.General.requestName = 'WJetsToLNu_PU20bx25_PHYS14'
#config.General.requestName = 'GluGluToHToWW_PU20bx25_PHYS14'
#config.General.requestName = 'QCD_Pt-20toInf_PU20bx25_PHYS14'
#config.General.requestName = 'TTDMDMJets_M1000GeV_PU20bx25_PHYS14'
config.General.requestName = 'TTDMDMJets_M600GeV_PU20bx25_PHYS14'
#config.General.requestName = 'TTDMDMJets_M200GeV_PU20bx25_PHYS14'
#config.General.requestName = 'TTDMDMJets_M100GeV_PU20bx25_PHYS14'
#config.General.requestName = 'TTDMDMJets_M50GeV_PU20bx25_PHYS14'
#config.General.requestName = 'TTDMDMJets_M10GeV_PU20bx25_PHYS14'
#config.General.requestName = 'TTDMDMJets_M1GeV_PU20bx25_PHYS14'

config.section_("JobType")
config.JobType.pluginName  = 'Analysis'
config.JobType.psetName    = 'skimToTreeSUSYMCtfs.py'
config.JobType.outputFiles = ['TreeSUSYtfs.root']
config.JobType.inputFiles  = ['Winter14_V5_DATA_UncertaintySources_AK5PFchs.txt','PHYS14_V2_MC_L3Absolute_AK4PFchs.txt','PHYS14_V2_MC_L2Relative_AK4PFchs.txt','PHYS14_V2_MC_L1FastJet_AK4PFchs.txt']

config.section_("Data")
#config.Data.inputDataset    = '/DYJetsToLL_M-50_13TeV-madgraph-pythia8/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM'
#config.Data.inputDataset    = '/WJetsToLNu_13TeV-madgraph-pythia8-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM'
#config.Data.inputDataset    = '/GluGluToHToWWTo2LAndTau2Nu_M-125_13TeV-powheg-pythia6/Phys14DR-PU20bx25_tsg_PHYS14_25_V1-v1/MINIAODSIM'
#config.Data.inputDataset    = '/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/Phys14DR-PU20bx25_PHYS14_25_V1-v3/MINIAODSIM'
#config.Data.inputDataset    = '/TTDMDMJets_M1000GeV_Tune4C_13TeV-madgraph-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM'
config.Data.inputDataset    = '/TTDMDMJets_M600GeV_Tune4C_13TeV-madgraph-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM'
#config.Data.inputDataset    = '/TTDMDMJets_M200GeV_Tune4C_13TeV-madgraph-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM'
#config.Data.inputDataset    = '/TTDMDMJets_M100GeV_Tune4C_13TeV-madgraph-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM'
#config.Data.inputDataset    = '/TTDMDMJets_M50GeV_Tune4C_13TeV-madgraph-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM'
#config.Data.inputDataset    = '/TTDMDMJets_M10GeV_Tune4C_13TeV-madgraph-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM'
#config.Data.inputDataset    = '/TTDMDMJets_M1GeV_Tune4C_13TeV-madgraph-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM'
config.Data.splitting       = 'FileBased'
config.Data.unitsPerJob     = 1
config.Data.publishDataName = 'Skim_2L_Pt17_8'
#config.Data.publishDataName = 'NoSkim'
config.Data.ignoreLocality  = True

config.section_("Site")
config.Site.storageSite = 'T2_ES_IFCA'

config.section_("User")

