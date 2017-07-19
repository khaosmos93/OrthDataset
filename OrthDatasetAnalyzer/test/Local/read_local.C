#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TTree.h"
#include "TBranch.h"

#include <iostream>
#include <algorithm>
#include <vector>
#include <string>
#include <map>

#include "DataFormats/Math/interface/deltaR.h"
//#include "MuonHLTNtuples/Analyzers/src/MuonTree.h"
#include "TLorentzVector.h"

double muonmass = 0.10565837;

const int nZPE = 293;
unsigned long long ZPrimeEvents[nZPE] = { 80025299, 102044699,1807633414,1751440345, 168766305, 707174799, 348062871,1502197004, 737696054,2158744629, 163837656,  46575796, 557743636, 129462036, 731395790,2036221243,  29437063, 461506595,2273760003, 724035648, 683941987, 641035873,2570923722, 465718408, 822904207,1491309805,1317604762, 224391879,2141455575, 172333402,4228668678,1043175650, 646321165, 233511973, 781278070, 164328984, 132273293, 398097480,2040124327,1291880061, 148951530, 886135059,  87415813,2561273750, 456814847,1089110527,3430307394,1133923162, 627297837, 668638994,2587883653,2794725084, 946606979, 839648380, 564315872, 424129031,1148427838,1398125706,  56502151, 545944592, 117074036, 207268415, 108863580, 882936293,  19644572, 526964876, 905350875,2235164626, 559385294,1563214828,1651139530, 593031982, 998852388,3225009196,3497296748, 245134853, 262790936,1177952475,1876459698, 833087925, 449811753,1697622247,1408454291,1398441723,1036218642,1703315077, 794052368,1933026733, 790119897, 133061470,1027774646, 103971827, 111248874, 995449716,1501620028,1339179910,1496315538, 293525896,1473721426, 288081692,   2240836, 362224985, 197808469,1612249125,2723865415,2941202420, 106152127, 206296298, 100641251,2332075909, 509311590,1249383398,1210213994,1090725731,1212340877,1901777118,3714674215, 400698612,  92906618,  95510738, 752260265,1196038009,2599229827, 330541383,  87694338, 313121685,  95584845,2162101314,3521968843,3392997950,2643079256,  59656496, 522576774, 324281834, 848563973, 222492431,1215395982, 256484428,2173236349,  35163519, 434635045, 198965344, 148915341,1777648427, 184352359, 475663239,2290869084, 492869259, 131269573,3053197353, 591789001,2211340369, 379265037, 281518062, 501871737,2099116360, 358838863, 351378620,1471282253, 315461884, 209797504,  39332041, 725487265, 950975177,1718269077,2072013665,2004948916,1075020228, 690551089,3529794117, 175965898,1539125792,  69333665, 952236920, 820405788, 108057469,  91857941, 327844969,1072719684, 983845763, 945261250, 130795909,2041970533,1658235670, 841713264, 483321692,2253330273, 746978377, 736995811, 232277083,  91005084,1141189789, 452707854, 328472356, 347808395, 424347292,  69308815,   6064394,  37448457, 653231214,2352609977,2064298422, 567286061,1842698514, 402153206,1354128049, 578279025,3616505346, 964177485,1768269831,1251706278, 103183615,2102021816,1487573631, 739088741,1728533708, 957779359,1922115605, 173190891,  59248890,2007612053,1224689298, 884488198,1203216999, 235704275,1444175555,1296908948,1576123886,  36191395,1990946855,1821429837, 858081924, 119831719, 729492867,2044961549, 402967394,1491796940, 210004024,  49382263,1783486685,1400221353,1481690384,2002900343,2262180829,2832916658, 429190704, 116935228, 318522998, 106461317,2650154968,3027061354, 284789350,4191779446,2616312606, 421956582,3932050225, 143994623, 253705534,2775188492,2334291201,  99430735, 456192114, 319185452,  93532853,  75926460,1070669073,2559866537,1494357742,2833007656,1239422009,1488881347,2778202581,1988145943, 106254202,1911653832,  59358666, 668271901, 485660009,1969855121,2436893471,1708799576, 695032770, 537500969,1441952618,1009753919, 310415205, 550203932, 327531822, 330535165, 978044075, 422723396, 647974202, 324709216 };


void read_local() {

  TFile* inputfile = TFile::Open("OrthDatasetTree_SingleMuon_MaybeG.root", "READ");
  std::cout << "input file: " << inputfile -> GetName() << std::endl;

  TTree *tree = (TTree*) inputfile -> Get("OrthDataset/EvTree");
  if (!tree) {
    std::cout << " *** tree not found *** " << std::endl;
    return;
  }

  Int_t runNum;
  Int_t lumiBlock;
  unsigned long long evtNum;
  Int_t nMuon;
  Int_t isMu50;
  Int_t isTkMu50;
  Int_t isIsoMu;
  Int_t isIsoTkMu;
  Int_t isMuTrg;
  std::vector<std::string> *ptrgPaths = 0;
  Double_t Dimuonmass;
  Double_t LeadingMuPt;
  Double_t LeadingMuEta;
  Double_t LeadingMuPhi;
  Double_t SubLeadingMuPt;
  Double_t SubLeadingMuEta;
  Double_t SubLeadingMuPhi;

  tree->SetBranchStatus("*", 1);

  tree->SetBranchAddress("runNum", &runNum);
  tree->SetBranchAddress("lumiBlock", &lumiBlock);
  tree->SetBranchAddress("evtNum", &evtNum);
  tree->SetBranchAddress("nMuon", &nMuon);
  tree->SetBranchAddress("isMu50", &isMu50);
  tree->SetBranchAddress("isTkMu50", &isTkMu50);
  tree->SetBranchAddress("isIsoMu", &isIsoMu);
  tree->SetBranchAddress("isIsoTkMu", &isIsoTkMu);
  tree->SetBranchAddress("isMuTrg", &isMuTrg);
  tree->SetBranchAddress("trgPaths",&ptrgPaths);
  tree->SetBranchAddress("Dimuonmass", &Dimuonmass);
  tree->SetBranchAddress("LeadingMuPt",&LeadingMuPt);
  tree->SetBranchAddress("LeadingMuEta",&LeadingMuEta);
  tree->SetBranchAddress("LeadingMuPhi",&LeadingMuPhi);
  tree->SetBranchAddress("SubLeadingMuPt",&SubLeadingMuPt);
  tree->SetBranchAddress("SubLeadingMuEta",&SubLeadingMuEta);
  tree->SetBranchAddress("SubLeadingMuPhi",&SubLeadingMuPhi);


  //TBranch*  evBranch = tree->GetBranch("");
  //evBranch -> SetAddress(&ev);


  int nentries = tree->GetEntriesFast();
  std::cout << "Number of entries = " << nentries << std::endl;

  for (Int_t eventNo=0; eventNo < nentries; eventNo++)
  {
    Int_t IgetEvent   = tree   -> GetEvent(eventNo);

    bool isZPrimeTrggerFire = false;
    bool isIsoTriggerFire = false;
    bool isMuTriggerFire = false;
    if ( isMu50==1 || isTkMu50==1 ) isZPrimeTrggerFire = true;
    if ( isIsoMu==1 || isIsoTkMu==1 ) isIsoTriggerFire = true;
    if ( isMuTrg==1 ) isMuTriggerFire = true;

    bool isInZPrimeList = false;
    for (Int_t ZPE = 0; ZPE < nZPE; ++ZPE) {
      if ( evtNum == ZPrimeEvents[ZPE]) isInZPrimeList = true;
    }

    std::cout << std::endl;
    std::cout << "[Event " << eventNo << "] " << runNum << ":" << lumiBlock << ":" << evtNum << std::endl;
    std::cout << "\tInv. Mass= " << Dimuonmass << std::endl;
    std::cout << "\t\tL.pT= " << LeadingMuPt << "\tL.eta= " << LeadingMuEta << "\tL.phi= " << LeadingMuPhi << std::endl;
    std::cout << "\t\tS.pT= " << SubLeadingMuPt << "\tS.eta= " << SubLeadingMuEta << "\tS.phi= " << SubLeadingMuPhi << std::endl;
    std::cout << "\tisMuTriggerFire : " << isMuTriggerFire << "\tisZPrimeTrggerFire : " << isZPrimeTrggerFire << "\tisIsoTriggerFire : " << isIsoTriggerFire << "\tisInZPrimeList : " << isInZPrimeList << std::endl;

    std::cout << "||Fired trigger paths||\n";

    for (int i=0; i<ptrgPaths->size(); ++i) {
    //for (std::vector<std::string>::const_iterator it = ptrgPaths.begin(); it != ptrgPaths.end(); ++it ) {
      std::string path = (*ptrgPaths)[i];
      if ( (path.find("Mu") !=std::string::npos) || (path.find("mu")) ) std::cout << path << ", ";
    }
    std::cout << std::endl;

  }

}
