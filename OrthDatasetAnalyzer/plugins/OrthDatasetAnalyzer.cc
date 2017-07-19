// -*- C++ -*-
//
// Package:    OrthDataset/OrthDatasetAnalyzer
// Class:      OrthDatasetAnalyzer
//
/**\class OrthDatasetAnalyzer OrthDatasetAnalyzer.cc OrthDataset/OrthDatasetAnalyzer/plugins/OrthDatasetAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  OH Minseok
//         Created:  Wed, 12 Jul 2017 16:10:19 GMT
//
//


#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include <iostream>
#include <fstream>
#include <cmath>

//////////////////////////
// -- Track & Vertex -- //
//////////////////////////
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "DataFormats/PatCandidates/interface/Vertexing.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexTools/interface/InvariantMassFromVertex.h"
#include "PhysicsTools/RecoUtils/interface/CandCommonVertexFitter.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/VertexPrimitives/interface/CachingVertex.h"

/////////////////////
// -- For Muons -- //
/////////////////////
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/L1GlobalMuonTrigger/interface/L1MuGMTCand.h"
#include "DataFormats/MuonReco/interface/MuonTimeExtra.h"
#include "DataFormats/MuonReco/interface/MuonTimeExtraMap.h"
#include "DataFormats/MuonReco/interface/MuonCosmicCompatibility.h"
#include "DataFormats/MuonReco/interface/MuonCocktails.h"

///////////////////
// -- For MET -- //
///////////////////
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"

////////////////////
// -- Triggers -- //
////////////////////
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
//#include "DataFormats/PatCandidates/interface/TriggerPrimitive.h"

////////////////
// -- Else -- //
////////////////
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/GeometrySurface/interface/GloballyPositioned.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositFwd.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositVetos.h"
#include "PhysicsTools/IsolationAlgos/interface/IsoDepositVetoFactory.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositDirection.h"
#include "DataFormats/PatCandidates/interface/Isolation.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include <TLorentzVector.h>
#include <TTree.h>

// MiniAOD
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

//
// class declaration
//

using namespace std;
using namespace reco;
using namespace edm;
using namespace pat;
using namespace isodeposit;


// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class OrthDatasetAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit OrthDatasetAnalyzer(const edm::ParameterSet&);
      ~OrthDatasetAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      TTree *EvTree;

      bool verbose_;
      //std::string outputFile_;
      double minMass_;

      int nEvt;
      int nEvt_minMass;

      // -- Tree -- //
      int runNum;
      unsigned long long evtNum;
      int lumiBlock;
      int nMuon;

      int isMu50;
      int isTkMu50;
      int isIsoMu;
      int isIsoTkMu;
      int isMuTrg;
      std::vector<std::string> trgPaths;
      std::vector<std::string> *ptrgPaths;

      double Dimuonmass;

      double LeadingMuPt;
      double LeadingMuEta;
      double LeadingMuPhi;

      double SubLeadingMuPt;
      double SubLeadingMuEta;
      double SubLeadingMuPhi;

      // ----------member data ---------------------------
      edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
      edm::EDGetTokenT<pat::MuonCollection> patMuonToken_;
      edm::EDGetTokenT<reco::VertexCollection> recoVertexToken_;
};

OrthDatasetAnalyzer::OrthDatasetAnalyzer(const edm::ParameterSet& iConfig)
  : //outputFile_(iConfig.getParameter<std::string>("OutputFile")),
    minMass_(iConfig.getParameter<double>("MinMass"))
{

  nEvt = 0;
  nEvt_minMass =0;

  triggerBits_      = consumes<edm::TriggerResults>        (edm::InputTag("TriggerResults","","HLT"));
  patMuonToken_     = consumes<pat::MuonCollection>        (edm::InputTag("slimmedMuons","","PAT"));
  recoVertexToken_  = consumes<reco::VertexCollection>     (edm::InputTag("offlineSlimmedPrimaryVertices","","PAT"));   //MiniAOD

  verbose_ = iConfig.getParameter<bool>("Verbose");

  //usesResource("TFileService");
}

OrthDatasetAnalyzer::~OrthDatasetAnalyzer(){}


void OrthDatasetAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  if (verbose_) std::cout << std::endl << "========== Analysing EVENT: " << iEvent.id() << " ====" << std::endl;

  //std::ofstream fout(outputFile_, std::ofstream::app);

  // Initialize
  nMuon = -1;

  isMu50 = 0;
  isTkMu50 = 0;
  isIsoMu = 0;
  isIsoTkMu = 0;
  isMuTrg = 0;

  Dimuonmass = -1;

  LeadingMuPt = -1;
  LeadingMuEta = -1;
  LeadingMuPhi = -1;

  SubLeadingMuPt = -1;
  SubLeadingMuEta = -1;
  SubLeadingMuPhi = -1;

  edm::Handle <pat::MuonCollection> patMuons;
  iEvent.getByToken(patMuonToken_,patMuons);
  edm::Handle <reco::VertexCollection> recoVertex;
  iEvent.getByToken(recoVertexToken_,recoVertex);
  edm::Handle <edm::TriggerResults> triggerBits;
  iEvent.getByToken(triggerBits_,triggerBits);

  // -- Vertex Information -- //
  reco::VertexCollection goodVertices;
  for (reco::VertexCollection::const_iterator it = recoVertex->begin(); it != recoVertex->end(); ++it) {
    if (it->ndof() > 4 &&
        fabs(it->z()) <= 24 &&
        fabs(it->position().rho()) <= 2) {
          goodVertices.push_back(*it);
    }
  }
  //if (verbose_) std::cout << "goodVertices.size(): " << goodVertices.size() << std::endl;

  if (goodVertices.size()==0){
          if (verbose_) std::cout << "WARNING!: no goodVertices - skipping event" << std::endl;
          return;
  }

  // -- Single Muon Selection -- //
  //const pat::Muon* SelMuon;
  pat::MuonCollection SelMuons;
  for (pat::MuonCollection::const_iterator it_mu = patMuons->begin(); it_mu != patMuons->end(); ++it_mu) {
    if (  it_mu->isGlobalMuon() and
          it_mu->isTrackerMuon() and
          it_mu->globalTrack().isNonnull() and
          it_mu->globalTrack()->hitPattern().numberOfValidMuonHits() > 0 and
          it_mu->globalTrack()->hitPattern().numberOfValidPixelHits() > 0 and
          it_mu->globalTrack()->hitPattern().trackerLayersWithMeasurement() > 5 and
          it_mu->muonBestTrack().isNonnull() and
          (it_mu->muonBestTrack()->ptError()/it_mu->muonBestTrack()->pt()) < 0.3 and
          (fabs(it_mu->muonBestTrack()->dxy(goodVertices.at(0).position())) < 0.2 ) and
          it_mu->pt()>53 and
          it_mu->innerTrack().isNonnull() and
          (it_mu->isolationR03().sumPt/it_mu->innerTrack()->pt())<0.1 and
          //it_mu->numberOfMatchedStations()>1
          (it_mu->numberOfMatchedStations()>1 or
          (it_mu->numberOfMatchedStations()==1 and !(it_mu->stationMask()==1 or it_mu->stationMask()==16)) or
          ((it_mu->numberOfMatchedStations()==1 and (it_mu->stationMask()==1 or it_mu->stationMask()==16)) and it_mu->numberOfMatchedRPCLayers()>2))
      ) {
      if (verbose_) std::cout << "\tMuon passes selection" << std::endl;
      SelMuons.push_back(*it_mu);
    }
  }

  if (SelMuons.size() < 2) {
    if (verbose_) std::cout << "WARNING!: Less than 2 muons - skipping event" << std::endl;
    return;
  }

  // -- Dimoun Selection -- //
  std::vector<pat::MuonCollection> DimuonPairs;
  for (pat::MuonCollection::const_iterator it1 = SelMuons.begin(); it1 != SelMuons.end(); ++it1) {
    for (pat::MuonCollection::const_iterator it2 = it1; it2 != SelMuons.end(); ++it2) {

      pat::MuonCollection TempPair;
      TempPair.clear();
      if (it1 == it2) continue;
      double cos_angle = (it1->px()*it2->px() + it1->py()*it2->py() + it1->pz()*it2->pz())/(it1->p()*it2->p());
      int Multi_charge = it1->charge()*it2->charge();
      if (cos_angle > -0.9998 and Multi_charge == -1) {
        TempPair.push_back(*it1);
        TempPair.push_back(*it2);
        DimuonPairs.push_back(TempPair);
      }

    }
  }
  if (verbose_) std::cout << "\tnumber of dimuon pairs passing selection: " << DimuonPairs.size() << std::endl;
  if (DimuonPairs.size() < 1) {
    if (verbose_) std::cout << "WARNING!: Less than 1 muon pair selected - skipping event" << std::endl;
    return;
  }

  //Highest Mass pair
  double TempMass = -999;
  double HighestMass = -999;
  unsigned int thePair = -999;
  for (unsigned int iPair=0; iPair < DimuonPairs.size(); ++iPair) {
    if (DimuonPairs[iPair].size() != 2) {
      if (verbose_) std::cout << "\t\tDimuon not in pair collection ?!?!" << std::endl;
    }
    auto&& mu0=DimuonPairs[iPair].at(0);
    auto&& mu1=DimuonPairs[iPair].at(1);
    TempMass = sqrt( pow(mu0.energy()+mu1.energy(),2) - pow(mu0.px()+mu1.px(),2) - pow(mu0.py()+mu1.py(),2) - pow(mu0.pz()+mu1.pz(),2) );
    if (TempMass > HighestMass) {
      HighestMass = TempMass;
      thePair=iPair;
    }
  }
  if (verbose_) std::cout << "\t\tChoosing Highest mass pair  highestMass= " << HighestMass << std::endl;

  pat::MuonCollection HighestMassPair;
  HighestMassPair.clear();
  HighestMassPair = DimuonPairs[thePair];
  auto&& HMmu0=HighestMassPair.at(0);
  auto&& HMmu1=HighestMassPair.at(1);
  auto&& Lmu=HMmu0;
  auto&& Smu=HMmu1;
  //auto&& TempMu;
  if (HMmu0.pt() < HMmu1.pt()) {
    Lmu=HMmu1;
    Smu=HMmu0;
  }
  if (verbose_) std::cout << "\t\tChoosing Leading and sub-leading muon  L.pT= " << Lmu.pt() << " S.pT= " << Smu.pt() << std::endl;



  nEvt++;

  double dimuonMass_=sqrt( pow(Lmu.energy()+Smu.energy(),2) - pow(Lmu.px()+Smu.px(),2) - pow(Lmu.py()+Smu.py(),2) - pow(Lmu.pz()+Smu.pz(),2) );

  if (verbose_) std::cout << "\t\tWriting Event ID" << std::endl;
  // -- Print High Mass event's ID , Fill only high mass events-- //
  if (dimuonMass_ > minMass_) {

    nEvt_minMass++;
    runNum = iEvent.id().run();
    lumiBlock = iEvent.id().luminosityBlock();
    evtNum = iEvent.id().event();
    //if (verbose_) std::cout << evtNum << std::endl;
    nMuon = patMuons->size();

    Dimuonmass = dimuonMass_;

    LeadingMuPt = Lmu.pt();
    LeadingMuEta = Lmu.eta();
    LeadingMuPhi = Lmu.phi();

    SubLeadingMuPt = Smu.pt();
    SubLeadingMuEta = Smu.eta();
    SubLeadingMuPhi = Smu.phi();


    ptrgPaths = &trgPaths;
    ptrgPaths->clear();
    const edm::TriggerNames &triggerNames_ = iEvent.triggerNames(*triggerBits);
    for (unsigned int itrig=0; itrig < triggerBits->size(); ++itrig) {
      //if (verbose_) std::cout << "Trigger path : " << triggerNames_.triggerName(itrig) << " : " << triggerBits->accept(itrig) << std::endl;
      if (triggerBits->accept(itrig)) {

        std::string pathName = triggerNames_.triggerName(itrig);
        ptrgPaths->push_back(pathName);
        if (verbose_) std::cout << "Trigger path : " << pathName << std::endl;

        if (pathName.find("HLT_Mu50") !=std::string::npos) isMu50 = 1;
        //else isMu50 = 0;
        if (pathName.find("HLT_TkMu50") !=std::string::npos) isTkMu50 = 1;
        //else isTkMu50 = 0;
        if (pathName.find("HLT_IsoMu") !=std::string::npos) isIsoMu = 1;
        //else isIsoMu = 0;
        if (pathName.find("HLT_IsoTkMu") !=std::string::npos) isIsoTkMu = 1;
        //else isIsoTkMu = 0;
        if ( (pathName.find("Mu") !=std::string::npos) || (pathName.find("mu") !=std::string::npos) ) isMuTrg = 1;
        //else isMuTrg = 0;
      }
    }

    EvTree->Fill();    //Tree
  }

  if (verbose_) std::cout << "\t\t\tFinish this event!!" << std::endl;

}


void OrthDatasetAnalyzer::beginJob() {

  //Tree
  edm::Service<TFileService> fs;
  EvTree = fs->make<TTree>("EvTree","EvTree");
  EvTree->Branch("nEvents_minMass",&nEvt_minMass,"nEvents_minMass/I");
  EvTree->Branch("runNum",&runNum,"runNum/I");
  EvTree->Branch("lumiBlock",&lumiBlock,"lumiBlock/I");
  EvTree->Branch("evtNum",&evtNum,"evtNum/l");
  EvTree->Branch("nMuon",&nMuon,"nMuon/I");
  EvTree->Branch("isMu50",&isMu50,"isMu50/I");
  EvTree->Branch("isTkMu50",&isTkMu50,"isTkMu50/I");
  EvTree->Branch("isIsoMu",&isIsoMu,"isIsoMu/I");
  EvTree->Branch("isIsoTkMu",&isIsoTkMu,"isIsoTkMu/I");
  EvTree->Branch("isMuTrg",&isMuTrg,"isMuTrg/I");
  EvTree->Branch("trgPaths","std::vector<std::string>",&ptrgPaths);
  EvTree->Branch("Dimuonmass",&Dimuonmass,"Dimuonmass/D");
  EvTree->Branch("LeadingMuPt",&LeadingMuPt,"LeadingMuPt/D");
  EvTree->Branch("LeadingMuEta",&LeadingMuEta,"LeadingMuEta/D");
  EvTree->Branch("LeadingMuPhi",&LeadingMuPhi,"LeadingMuPhi/D");
  EvTree->Branch("SubLeadingMuPt",&SubLeadingMuPt,"SubLeadingMuPt/D");
  EvTree->Branch("SubLeadingMuEta",&SubLeadingMuEta,"SubLeadingMuEta/D");
  EvTree->Branch("SubLeadingMuPhi",&SubLeadingMuPhi,"SubLeadingMuPhi/D");

}

void OrthDatasetAnalyzer::endJob(){
  std::cout <<"++++++++++++++++++++++++++++++++++++++" << std::endl;
  std::cout <<"analyzed " << nEvt << " dimuon events: " << std::endl;
  std::cout <<"analyzed " << nEvt_minMass << " dimuon events with mass above " << minMass_ << std::endl;
  std::cout <<"++++++++++++++++++++++++++++++++++++++" << std::endl;
}

void OrthDatasetAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(OrthDatasetAnalyzer);
