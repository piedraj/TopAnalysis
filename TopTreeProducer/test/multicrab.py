import os
from WMCore.Configuration import Configuration
config = Configuration()

config.section_('General')
config.General.transferLogs = True
config.General.workArea     = 'crab_7May_PHYS14'

config.section_('JobType')
config.JobType.pluginName  = 'Analysis'
config.JobType.psetName    = 'skimToTreeSUSYMCtfs.py'
config.JobType.outputFiles = ['Tree_13TeV.root']
config.JobType.inputFiles  = ['Winter14_V5_DATA_L2L3Residual_AK5PFchs.txt','Winter14_V5_DATA_L3Absolute_AK5PFchs.txt','Winter14_V5_DATA_L2Relative_AK5PFchs.txt','Winter14_V5_DATA_L1FastJet_AK5PFchs.txt','Summer13_V5_DATA_UncertaintySources_AK5PFchs.txt','PHYS14_V2_MC_L3Absolute_AK4PFchs.txt','PHYS14_V2_MC_L2Relative_AK4PFchs.txt','PHYS14_V2_MC_L1FastJet_AK4PFchs.txt']

config.section_('Data')    
config.Data.inputDBS        = 'global'
config.Data.splitting       = 'FileBased'
config.Data.unitsPerJob     = 1
config.Data.publishDataName = 'Skim_2L_Pt17_8'
config.Data.ignoreLocality  = True

config.section_('Site')
config.Site.storageSite = 'T2_ES_IFCA'


if __name__ == '__main__':

    from CRABAPI.RawCommand import crabCommand

    def submit(config):
        res = crabCommand('submit', config = config)

    ####################### Public samples to be analysed ######################
                   
    config.General.requestName = 'TTDMDMJets_M1GeV'
    config.Data.inputDataset   = '/TTDMDMJets_M1GeV_Tune4C_13TeV-madgraph-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM'
    submit(config)

    config.General.requestName = 'TTDMDMJets_M10GeV'
    config.Data.inputDataset   = '/TTDMDMJets_M10GeV_Tune4C_13TeV-madgraph-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM'
    submit(config)

    config.General.requestName = 'TTDMDMJets_M50GeV'
    config.Data.inputDataset   = '/TTDMDMJets_M50GeV_Tune4C_13TeV-madgraph-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM'
    submit(config)

    config.General.requestName = 'TTDMDMJets_M100GeV'
    config.Data.inputDataset   = '/TTDMDMJets_M100GeV_Tune4C_13TeV-madgraph-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM'
    submit(config)

    config.General.requestName = 'TTDMDMJets_M200GeV'
    config.Data.inputDataset   = '/TTDMDMJets_M200GeV_Tune4C_13TeV-madgraph-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM'
    submit(config)

    config.General.requestName = 'TTDMDMJets_M600GeV'
    config.Data.inputDataset   = '/TTDMDMJets_M600GeV_Tune4C_13TeV-madgraph-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM'
    submit(config)

    config.General.requestName = 'TTDMDMJets_M1000GeV'
    config.Data.inputDataset   = '/TTDMDMJets_M1000GeV_Tune4C_13TeV-madgraph-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM'
    submit(config)
