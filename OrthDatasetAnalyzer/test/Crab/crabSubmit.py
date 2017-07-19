from CRABClient.UserUtilities import config, getUsernameFromSiteDB
import sys, os
import datetime

now = datetime.datetime.now()
date = now.strftime('%Y%m%d')

submitVersion = 'OrthDataset' ### modified
mainOutputDir = '/store/user/%s/%s' % (getUsernameFromSiteDB(),submitVersion) ### modified
PrimaryDataset = ['JetHT','MET','SingleElectron','SinglePhoton']

DataSample = [
                 #('Bv1','/Run2016B-03Feb2017_ver1-v1/MINIAOD',''),
                 ('RunBv2','/Run2016B-03Feb2017_ver2-v2/MINIAOD','80X_dataRun2_2016SeptRepro_v7'),
                 ('RunC','/Run2016C-03Feb2017-v1/MINIAOD','80X_dataRun2_2016SeptRepro_v7'),
                 ('RunD','/Run2016D-03Feb2017-v1/MINIAOD','80X_dataRun2_2016SeptRepro_v7'),
                 ('RunE','/Run2016E-03Feb2017-v1/MINIAOD','80X_dataRun2_2016SeptRepro_v7'),
                 ('RunF','/Run2016F-03Feb2017-v1/MINIAOD','80X_dataRun2_2016SeptRepro_v7'),
                 ('RunG','/Run2016G-03Feb2017-v1/MINIAOD','80X_dataRun2_2016SeptRepro_v7'),
                 ('RunHv2','/Run2016H-03Feb2017_ver2-v1/MINIAOD','80X_dataRun2_Prompt_v16'),
                 ('RunHv3','/Run2016H-03Feb2017_ver3-v1/MINIAOD','80X_dataRun2_Prompt_v16'),
              ]


MCSample = [
                 #('DY800_1400','/ZToMuMu_NNPDF30_13TeV-powheg_M_800_1400/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'),
                 #('DY1400_2300','/ZToMuMu_NNPDF30_13TeV-powheg_M_1400_2300/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
            ]



if 'Data' in sys.argv:
  Sample = DataSample
if 'MC' in sys.argv:
  Sample = MCSample

if __name__ == '__main__':  # and 'submit' in sys.argv:
  print "Date : ", now

  crab_cfg = '''
from WMCore.Configuration import Configuration
config = Configuration()

config.section_('General')
config.General.transferOutputs = True
config.General.transferLogs = False
config.General.workArea = 'crab_OrthDataset_%(PDataset)s_%(date)s'
config.General.requestName = 'OrthDataset_%(FullName)s_%(date)s'

config.section_('JobType')
config.JobType.psetName = 'run_crab.py'
config.JobType.pluginName = 'Analysis'
#config.JobType.disableAutomaticOutputCollection = True
#config.JobType.outputFiles = ['OrthDatasetDimuonEvents_M900.txt','OrthDatasetHistos.root']
config.JobType.numCores = 8
config.JobType.maxMemoryMB = 2500
config.JobType.maxJobRuntimeMin = 2500

config.section_('Data')
config.Data.outLFNDirBase = '%(mainOutputDir)s/'
config.Data.inputDataset = '%(datasetPath)s'
config.Data.publication = False
job_control
config.section_('User')

config.section_('Site')
config.Site.blacklist = ['T2_IT_Legnaro']
config.Site.storageSite = 'T2_KR_KNU'
#config.Site.storageSite = 'T3_KR_KISTI'
'''

  if 'Data' in sys.argv:
    crab_cfg = crab_cfg.replace('job_control','''
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 3000
config.Data.totalUnits = -1
config.Data.inputDBS = 'https://cmsweb.cern.ch/dbs/prod/global/DBSReader'
#config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON_MuonPhys.txt'
config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
''')

  if 'MC' in sys.argv:
    crab_cfg = crab_cfg.replace('job_control','''
config.Data.splitting = 'EventAwareLumiBased'
#config.Data.splitting = 'FileBased'
config.Data.totalUnits = -1
config.Data.unitsPerJob  = 20000
config.Data.inputDBS = 'https://cmsweb.cern.ch/dbs/prod/global/DBSReader'
''')

  if 'submit' in sys.argv:
    for PDataset in PrimaryDataset:
      print "\n\n", PDataset
      for name, dataset, GT in Sample:
        FullName = PDataset + '_' + name
        datasetPath = '/' + PDataset + dataset
        print "\n", datasetPath

        f = open('run.py', 'r')
        fData = f.read()
        f.close()
        newData = fData.replace('GlobalTagReplace',GT)
        newData = newData.replace('PDReplace',PDataset)
        newData = newData.replace('PeriodReplayce',name)
        fc = open('run_crab.py','w')
        fc.write(newData)
        fc.close()

        open('crabConfig.py', 'wt').write(crab_cfg % locals())
        #os.system('crab submit -c crabConfig.py')


