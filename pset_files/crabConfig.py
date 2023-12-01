import CRABClient
from CRABClient.UserUtilities import config

config = config()

config.General.requestName = 'Cosmics_2023C'
config.General.workArea = 'crab'
config.General.transferOutputs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'run3_data.py'
config.JobType.numCores = 8

config.Data.inputDataset = '/Cosmics/Run2023C-PromptReco-v1/AOD'
config.Data.inputDBS = 'global'
# config.Data.useParent = True
# config.Data.partialDataset = True
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 10
# config.Data.publication = True
config.Data.outputDatasetTag = 'Cosmics_2023C'
config.Site.storageSite = 'T3_CH_CERNBOX'
