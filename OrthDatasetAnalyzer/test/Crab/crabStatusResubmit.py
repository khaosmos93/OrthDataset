import sys, os
if __name__ == '__main__':

  dirs = [

                (' JetHT_Bv2 ','crab_OrthDataset_JetHT_20170714/crab_OrthDataset_JetHT_RunBv2_20170714'),
                (' JetHT_C ','crab_OrthDataset_JetHT_20170714/crab_OrthDataset_JetHT_RunC_20170714'),
                (' JetHT_D ','crab_OrthDataset_JetHT_20170714/crab_OrthDataset_JetHT_RunD_20170714'),
                (' JetHT_E ','crab_OrthDataset_JetHT_20170714/crab_OrthDataset_JetHT_RunE_20170714'),
                (' JetHT_F ','crab_OrthDataset_JetHT_20170714/crab_OrthDataset_JetHT_RunF_20170714'),
                (' JetHT_G ','crab_OrthDataset_JetHT_20170714/crab_OrthDataset_JetHT_RunG_20170714'),
                (' JetHT_Hv2','crab_OrthDataset_JetHT_20170714/crab_OrthDataset_JetHT_RunHv2_20170714'),
                (' JetHT_Hv3','crab_OrthDataset_JetHT_20170714/crab_OrthDataset_JetHT_RunHv3_20170714'),

                (' MET_Bv2 ','crab_OrthDataset_MET_20170714/crab_OrthDataset_MET_RunBv2_20170714'),
                (' MET_C ','crab_OrthDataset_MET_20170714/crab_OrthDataset_MET_RunC_20170714'),
                (' MET_D ','crab_OrthDataset_MET_20170714/crab_OrthDataset_MET_RunD_20170714'),
                (' MET_E ','crab_OrthDataset_MET_20170714/crab_OrthDataset_MET_RunE_20170714'),
                (' MET_F ','crab_OrthDataset_MET_20170714/crab_OrthDataset_MET_RunF_20170714'),
                (' MET_G ','crab_OrthDataset_MET_20170714/crab_OrthDataset_MET_RunG_20170714'),
                (' MET_Hv2','crab_OrthDataset_MET_20170714/crab_OrthDataset_MET_RunHv2_20170714'),
                (' MET_Hv3','crab_OrthDataset_MET_20170714/crab_OrthDataset_MET_RunHv3_20170714'),

                (' SingleElectron_Bv2 ','crab_OrthDataset_SingleElectron_20170714/crab_OrthDataset_SingleElectron_RunBv2_20170714'),
                (' SingleElectron_C ','crab_OrthDataset_SingleElectron_20170714/crab_OrthDataset_SingleElectron_RunC_20170714'),
                (' SingleElectron_D ','crab_OrthDataset_SingleElectron_20170714/crab_OrthDataset_SingleElectron_RunD_20170714'),
                (' SingleElectron_E ','crab_OrthDataset_SingleElectron_20170714/crab_OrthDataset_SingleElectron_RunE_20170714'),
                (' SingleElectron_F ','crab_OrthDataset_SingleElectron_20170714/crab_OrthDataset_SingleElectron_RunF_20170714'),
                (' SingleElectron_G ','crab_OrthDataset_SingleElectron_20170714/crab_OrthDataset_SingleElectron_RunG_20170714'),
                (' SingleElectron_Hv2','crab_OrthDataset_SingleElectron_20170714/crab_OrthDataset_SingleElectron_RunHv2_20170714'),
                (' SingleElectron_Hv3','crab_OrthDataset_SingleElectron_20170714/crab_OrthDataset_SingleElectron_RunHv3_20170714'),

                (' SinglePhoton_Bv2 ','crab_OrthDataset_SinglePhoton_20170714/crab_OrthDataset_SinglePhoton_RunBv2_20170714'),
                (' SinglePhoton_C ','crab_OrthDataset_SinglePhoton_20170714/crab_OrthDataset_SinglePhoton_RunC_20170714'),
                (' SinglePhoton_D ','crab_OrthDataset_SinglePhoton_20170714/crab_OrthDataset_SinglePhoton_RunD_20170714'),
                (' SinglePhoton_E ','crab_OrthDataset_SinglePhoton_20170714/crab_OrthDataset_SinglePhoton_RunE_20170714'),
                (' SinglePhoton_F ','crab_OrthDataset_SinglePhoton_20170714/crab_OrthDataset_SinglePhoton_RunF_20170714'),
                (' SinglePhoton_G ','crab_OrthDataset_SinglePhoton_20170714/crab_OrthDataset_SinglePhoton_RunG_20170714'),
                (' SinglePhoton_Hv2','crab_OrthDataset_SinglePhoton_20170714/crab_OrthDataset_SinglePhoton_RunHv2_20170714'),
                (' SinglePhoton_Hv3','crab_OrthDataset_SinglePhoton_20170714/crab_OrthDataset_SinglePhoton_RunHv3_20170714'),

          ]


  for name, di in dirs:

    name = '|| -- '+name+' -- ||'
    if 'status' in sys.argv:
      command = 'crab status -d ' + di #+ ' --long'
    elif 'resubmit' in sys.argv:
      command = 'crab resubmit -d ' + di
    elif 'kill' in sys.argv:
      command = 'crab kill -d ' + di

    print '\n',name
    os.system(command)
