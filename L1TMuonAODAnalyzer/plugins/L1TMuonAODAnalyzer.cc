// -*- C++ -*-
//
// Package:    L1TMuonAODAnalyzer/L1TMuonAODAnalyzer
// Class:      L1TMuonAODAnalyzer
//
/**\class L1TMuonAODAnalyzer L1TMuonAODAnalyzer.cc L1TMuonAODAnalyzer/L1TMuonAODAnalyzer/plugins/L1TMuonAODAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Efe Yigitbasi
//         Created:  Sat, 10 Sep 2022 11:08:53 GMT
//
//

#include "L1TMuonAODAnalyzer.h"


//
// constructors and destructor
//
L1TMuonAODAnalyzer::L1TMuonAODAnalyzer(const edm::ParameterSet& iConfig)
    :
    muonToken_(consumes< std::vector< reco::Muon> >(iConfig.getParameter<edm::InputTag>("Muons"))),
    standAloneMuonToken_(consumes< std::vector< reco::Track> >(iConfig.getParameter<edm::InputTag>("standAloneMuons"))),
    cosmicMuonToken_(consumes< std::vector< reco::Track> >(iConfig.getParameter<edm::InputTag>("cosmicMuons"))),
    l1MuonToken_(consumes<l1t::MuonBxCollection>(edm::InputTag("gmtStage2Digis","Muon"))),
    l1BMTFRegionalMuonCandToken_(consumes<BXVector<l1t::RegionalMuonCand>>(edm::InputTag("gmtStage2Digis","BMTF"))),
    verticesToken_(consumes<std::vector<Vertex> > (iConfig.getParameter<edm::InputTag>("Vertices"))),

    MuonPtCut_(iConfig.getParameter<double>("MuonPtCut")),
    SaveTree_(iConfig.getParameter<bool>("SaveTree")),
    IsMC_(iConfig.getParameter<bool>("IsMC")),
    Debug_(iConfig.getParameter<bool>("Debug")),

    muPropagatorSetup1st_(iConfig.getParameter<edm::ParameterSet>("muProp1st"), consumesCollector()),
    muPropagatorSetup2nd_(iConfig.getParameter<edm::ParameterSet>("muProp2nd"), consumesCollector()),
    results_(nullptr)

{
  //now do what ever initialization is needed
  // usesResource("TFileService"); // shared resources

  edm::Service<TFileService> fs;
  outputTree = fs->make<TTree>("tree","tree");

}

L1TMuonAODAnalyzer::~L1TMuonAODAnalyzer() {}

//
// member functions
//

// ------------ method called for each event  ------------
void L1TMuonAODAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  InitandClearStuff();

  muPropagator1st_ = muPropagatorSetup1st_.init(iSetup);
  muPropagator2nd_ = muPropagatorSetup2nd_.init(iSetup);


  _runNb = iEvent.id().run();
  _eventNb = iEvent.id().event();
  _lumiBlock = iEvent.luminosityBlock();
  _bx=iEvent.bunchCrossing();

  // L1 muons
    edm::Handle<l1t::MuonBxCollection> l1muoncoll;
    iEvent.getByToken(l1MuonToken_ , l1muoncoll);
    for(int i = l1muoncoll->getFirstBX() ; i<= l1muoncoll->getLastBX() ;i++){
      for( l1t::MuonBxCollection::const_iterator l1muonit= l1muoncoll->begin(i); l1muonit != l1muoncoll->end(i) ; ++l1muonit){
        if(l1muonit->pt() < 0) continue;
        l1mu_qual.push_back( l1muonit->hwQual() );
        l1mu_charge.push_back( l1muonit->charge() );
        l1mu_pt.push_back( l1muonit->pt() );
        l1mu_pt_dxy.push_back( l1muonit->ptUnconstrained() );
        l1mu_dxy.push_back( l1muonit->hwDXY() );
        l1mu_eta.push_back( l1muonit->eta() );
        l1mu_etaAtVtx.push_back( l1muonit->etaAtVtx() );
        l1mu_phi.push_back( l1muonit->phi() );
        l1mu_phiAtVtx.push_back( l1muonit->phiAtVtx() );
        l1mu_tfIdx.push_back(l1muonit->tfMuonIndex());

        l1mu_bx.push_back( i);
        l1mu_size++;

      }
    }

  // BMTF RegionalMuonCand
  edm::Handle<BXVector<l1t::RegionalMuonCand>> l1RegionalMuoncoll;
  iEvent.getByToken(l1BMTFRegionalMuonCandToken_ , l1RegionalMuoncoll);
  for(int i = l1RegionalMuoncoll->getFirstBX() ; i<= l1RegionalMuoncoll->getLastBX() ;i++){
    for( BXVector<l1t::RegionalMuonCand>::const_iterator l1muonit= l1RegionalMuoncoll->begin(i); l1muonit != l1RegionalMuoncoll->end(i) ; ++l1muonit){
      BMTFMu_processor.push_back(l1muonit->processor());
      BMTFMu_hwPt.push_back(l1muonit->hwPt());
      BMTFMu_hwQual.push_back(l1muonit->hwQual());
      BMTFMu_hwSign.push_back(l1muonit->hwSign());
      BMTFMu_hwSignValid.push_back(l1muonit->hwSignValid());
      BMTFMu_hwEta.push_back(l1muonit->hwEta());
      BMTFMu_hwPhi.push_back(l1muonit->hwPhi());
    }
  }


  //Vertices
  edm::Handle<std::vector<Vertex> > theVertices;
  iEvent.getByToken(verticesToken_,theVertices) ;
  _n_PV = theVertices->size();
  Vertex::Point PV(0,0,0);
  if(_n_PV){ PV = theVertices->begin()->position();}

  //Muons
  edm::Handle< std::vector<reco::Muon> > theMuons;
  iEvent.getByToken(muonToken_,theMuons);
  for( std::vector<reco::Muon>::const_iterator muon = (*theMuons).begin(); muon != (*theMuons).end(); muon++ ) {
    if((&*muon)->pt() <0) continue; //Loose cut  on uncorrected pt 

    double ptmuoncorr= (&*muon)->pt();

    // store all reco muons for now
    muon_size++;
    muon_eta.push_back((&*muon)->eta());
    muon_phi.push_back((&*muon)->phi());
    muon_pt.push_back((&*muon)->pt());
    muon_ptCorr.push_back( ptmuoncorr );
    muon_charge.push_back((&*muon)->charge());
    muon_PassTightID.push_back(  (&*muon)->passed(reco::Muon::CutBasedIdMediumPrompt )&& (&*muon)->passed(reco::Muon::PFIsoTight ) );
    muon_PassLooseID.push_back(  (&*muon)->passed(reco::Muon::CutBasedIdLoose )&& (&*muon)->passed(reco::Muon::PFIsoLoose ) );
    muon_isSAMuon.push_back( (&*muon)->isStandAloneMuon());

    if( !((&*muon)->innerTrack()).isNull()){
      muon_dxy.push_back( (&*muon)->innerTrack()->dxy(PV));
    }
    else{
      muon_dxy.push_back(-999.);
    }


    // extrapolation of muon track coordinates
    TrajectoryStateOnSurface stateAtMuSt1 = muPropagator1st_.extrapolate(*muon);
    if (stateAtMuSt1.isValid()) {
        muon_etaAtSt1.push_back(stateAtMuSt1.globalPosition().eta());
        muon_phiAtSt1.push_back(stateAtMuSt1.globalPosition().phi());
    } else {
        muon_etaAtSt1.push_back(-999);
        muon_phiAtSt1.push_back(-999);
    }

    TrajectoryStateOnSurface stateAtMuSt2 = muPropagator2nd_.extrapolate(*muon);
    if (stateAtMuSt2.isValid()) {
        muon_etaAtSt2.push_back(stateAtMuSt2.globalPosition().eta());
        muon_phiAtSt2.push_back(stateAtMuSt2.globalPosition().phi());
    } else {
        muon_etaAtSt2.push_back(-999);
        muon_phiAtSt2.push_back(-999);
    }
  }

  // standAlone Muons
  edm::Handle< std::vector<reco::Track> > theStandAloneMuons;
  iEvent.getByToken(standAloneMuonToken_,theStandAloneMuons);
  for( std::vector<reco::Track>::const_iterator muon = (*theStandAloneMuons).begin(); muon != (*theStandAloneMuons).end(); muon++ ) {
    if((&*muon)->pt() <0) continue; //Loose cut  on uncorrected pt 

    double ptmuoncorr= (&*muon)->pt();

    standAloneMuon_size++;
    standAloneMuon_eta.push_back((&*muon)->eta());
    standAloneMuon_phi.push_back((&*muon)->phi());
    standAloneMuon_pt.push_back((&*muon)->pt());
    standAloneMuon_ptCorr.push_back( ptmuoncorr );
    standAloneMuon_charge.push_back((&*muon)->charge());
    standAloneMuon_dxy.push_back( (&*muon)->dxy());

    // extrapolation of muon track coordinates
    TrajectoryStateOnSurface stateAtMuSt1 = muPropagator1st_.extrapolate(*muon);
    if (stateAtMuSt1.isValid()) {
        standAloneMuon_etaAtSt1.push_back(stateAtMuSt1.globalPosition().eta());
        standAloneMuon_phiAtSt1.push_back(stateAtMuSt1.globalPosition().phi());
    } else {
        standAloneMuon_etaAtSt1.push_back(-999);
        standAloneMuon_phiAtSt1.push_back(-999);
    }

    TrajectoryStateOnSurface stateAtMuSt2 = muPropagator2nd_.extrapolate(*muon);
    if (stateAtMuSt2.isValid()) {
        standAloneMuon_etaAtSt2.push_back(stateAtMuSt2.globalPosition().eta());
        standAloneMuon_phiAtSt2.push_back(stateAtMuSt2.globalPosition().phi());
    } else {
        standAloneMuon_etaAtSt2.push_back(-999);
        standAloneMuon_phiAtSt2.push_back(-999);
    }
  }

  // cosmicMuons
  edm::Handle< std::vector<reco::Track> > theCosmicMuons;
  iEvent.getByToken(cosmicMuonToken_,theCosmicMuons);
  for( std::vector<reco::Track>::const_iterator muon = (*theCosmicMuons).begin(); muon != (*theCosmicMuons).end(); muon++ ) {
    if((&*muon)->pt() <0) continue; //Loose cut  on uncorrected pt 

    double ptmuoncorr= (&*muon)->pt();

    // store all reco muons for now
    cosmicMuon_size++;
    cosmicMuon_eta.push_back((&*muon)->eta());
    cosmicMuon_phi.push_back((&*muon)->phi());
    cosmicMuon_pt.push_back((&*muon)->pt());
    cosmicMuon_ptCorr.push_back( ptmuoncorr );
    cosmicMuon_charge.push_back((&*muon)->charge());
    cosmicMuon_dxy.push_back( (&*muon)->dxy());

    // extrapolation of muon track coordinates
    TrajectoryStateOnSurface stateAtMuSt1 = muPropagator1st_.extrapolate(*muon);
    if (stateAtMuSt1.isValid()) {
        cosmicMuon_etaAtSt1.push_back(stateAtMuSt1.globalPosition().eta());
        cosmicMuon_phiAtSt1.push_back(stateAtMuSt1.globalPosition().phi());
    } else {
        cosmicMuon_etaAtSt1.push_back(-999);
        cosmicMuon_phiAtSt1.push_back(-999);
    }

    TrajectoryStateOnSurface stateAtMuSt2 = muPropagator2nd_.extrapolate(*muon);
    if (stateAtMuSt2.isValid()) {
        cosmicMuon_etaAtSt2.push_back(stateAtMuSt2.globalPosition().eta());
        cosmicMuon_phiAtSt2.push_back(stateAtMuSt2.globalPosition().phi());
    } else {
        cosmicMuon_etaAtSt2.push_back(-999);
        cosmicMuon_phiAtSt2.push_back(-999);
    }
  }

  if(SaveTree_)outputTree->Fill();

}

// ------------ method called once each job just before starting event loop  ------------
void L1TMuonAODAnalyzer::beginJob() {
  // please remove this method if not needed

  outputTree->Branch("_eventNb",   &_eventNb,   "_eventNb/l");
  outputTree->Branch("_runNb",     &_runNb,     "_runNb/l");
  outputTree->Branch("_lumiBlock", &_lumiBlock, "_lumiBlock/l");
  outputTree->Branch("_bx", &_bx, "_bx/l");


  outputTree->Branch("muon_eta",&muon_eta);
  outputTree->Branch("muon_etaAtSt1",&muon_etaAtSt1);
  outputTree->Branch("muon_etaAtSt2",&muon_etaAtSt2);
  outputTree->Branch("muon_phi",&muon_phi);
  outputTree->Branch("muon_phiAtSt1",&muon_phiAtSt1);
  outputTree->Branch("muon_phiAtSt2",&muon_phiAtSt2);
  outputTree->Branch("muon_pt",&muon_pt);
  outputTree->Branch("muon_ptCorr",&muon_ptCorr);
  outputTree->Branch("muon_charge",&muon_charge);

  outputTree->Branch("muon_dz",&muon_dz);
  outputTree->Branch("muon_dzError",&muon_dzError);
  outputTree->Branch("muon_dxy",&muon_dxy);
  outputTree->Branch("muon_dxyError",&muon_dxyError);
  outputTree->Branch("muon_3dIP",&muon_3dIP);
  outputTree->Branch("muon_3dIPError",&muon_3dIPError);

  outputTree->Branch("muon_PassTightID",&muon_PassTightID);
  outputTree->Branch("muon_PassLooseID",&muon_PassLooseID);
  outputTree->Branch("muon_isSAMuon",&muon_isSAMuon);

  outputTree->Branch("muon_size", &muon_size, "muon_size/I");

  outputTree->Branch("standAloneMuon_eta",&standAloneMuon_eta);
  outputTree->Branch("standAloneMuon_etaAtSt1",&standAloneMuon_etaAtSt1);
  outputTree->Branch("standAloneMuon_etaAtSt2",&standAloneMuon_etaAtSt2);
  outputTree->Branch("standAloneMuon_phi",&standAloneMuon_phi);
  outputTree->Branch("standAloneMuon_phiAtSt1",&standAloneMuon_phiAtSt1);
  outputTree->Branch("standAloneMuon_phiAtSt2",&standAloneMuon_phiAtSt2);
  outputTree->Branch("standAloneMuon_pt",&standAloneMuon_pt);
  outputTree->Branch("standAloneMuon_ptCorr",&standAloneMuon_ptCorr);
  outputTree->Branch("standAloneMuon_charge",&standAloneMuon_charge);
  outputTree->Branch("standAloneMuon_dz",&standAloneMuon_dz);
  outputTree->Branch("standAloneMuon_dzError",&standAloneMuon_dzError);
  outputTree->Branch("standAloneMuon_dxy",&standAloneMuon_dxy);
  outputTree->Branch("standAloneMuon_dxyError",&standAloneMuon_dxyError);
  outputTree->Branch("standAloneMuon_3dIP",&standAloneMuon_3dIP);
  outputTree->Branch("standAloneMuon_3dIPError",&standAloneMuon_3dIPError);
  outputTree->Branch("standAloneMuon_PassTightID",&standAloneMuon_PassTightID);
  outputTree->Branch("standAloneMuon_PassLooseID",&standAloneMuon_PassLooseID);
  outputTree->Branch("standAloneMuon_size", &standAloneMuon_size, "standAloneMuon_size/I");

  outputTree->Branch("cosmicMuon_eta",&cosmicMuon_eta);
  outputTree->Branch("cosmicMuon_etaAtSt1",&cosmicMuon_etaAtSt1);
  outputTree->Branch("cosmicMuon_etaAtSt2",&cosmicMuon_etaAtSt2);
  outputTree->Branch("cosmicMuon_phi",&cosmicMuon_phi);
  outputTree->Branch("cosmicMuon_phiAtSt1",&cosmicMuon_phiAtSt1);
  outputTree->Branch("cosmicMuon_phiAtSt2",&cosmicMuon_phiAtSt2);
  outputTree->Branch("cosmicMuon_pt",&cosmicMuon_pt);
  outputTree->Branch("cosmicMuon_ptCorr",&cosmicMuon_ptCorr);
  outputTree->Branch("cosmicMuon_charge",&cosmicMuon_charge);
  outputTree->Branch("cosmicMuon_dz",&cosmicMuon_dz);
  outputTree->Branch("cosmicMuon_dzError",&cosmicMuon_dzError);
  outputTree->Branch("cosmicMuon_dxy",&cosmicMuon_dxy);
  outputTree->Branch("cosmicMuon_dxyError",&cosmicMuon_dxyError);
  outputTree->Branch("cosmicMuon_3dIP",&cosmicMuon_3dIP);
  outputTree->Branch("cosmicMuon_3dIPError",&cosmicMuon_3dIPError);
  outputTree->Branch("cosmicMuon_PassTightID",&cosmicMuon_PassTightID);
  outputTree->Branch("cosmicMuon_PassLooseID",&cosmicMuon_PassLooseID);
  outputTree->Branch("cosmicMuon_size", &cosmicMuon_size, "cosmicMuon_size/I");

  outputTree->Branch("l1mu_qual",&l1mu_qual);
  outputTree->Branch("l1mu_charge",&l1mu_charge);
  outputTree->Branch("l1mu_pt",&l1mu_pt);
  outputTree->Branch("l1mu_pt_dxy",&l1mu_pt_dxy);
  outputTree->Branch("l1mu_dxy",&l1mu_dxy);
  outputTree->Branch("l1mu_eta",&l1mu_eta);
  outputTree->Branch("l1mu_etaAtVtx",&l1mu_etaAtVtx);
  outputTree->Branch("l1mu_phi",&l1mu_phi);
  outputTree->Branch("l1mu_phiAtVtx",&l1mu_phiAtVtx);
  outputTree->Branch("l1mu_tfIdx",&l1mu_tfIdx);
  outputTree->Branch("l1mu_bx",&l1mu_bx);
  outputTree->Branch("l1mu_size", &l1mu_size, "l1mu_size/I");

  outputTree->Branch("BMTFMu_processor",&BMTFMu_processor);
  outputTree->Branch("BMTFMu_hwPt",&BMTFMu_hwPt);
  outputTree->Branch("BMTFMu_hwQual",&BMTFMu_hwQual);
  outputTree->Branch("BMTFMu_hwSign",&BMTFMu_hwSign);
  outputTree->Branch("BMTFMu_hwSignValid",&BMTFMu_hwSignValid);
  outputTree->Branch("BMTFMu_hwEta",&BMTFMu_hwEta);
  outputTree->Branch("BMTFMu_hwPhi",&BMTFMu_hwPhi);

}

void L1TMuonAODAnalyzer::endJob() {}

void L1TMuonAODAnalyzer::InitandClearStuff() {

  muon_eta.clear();
  muon_etaAtSt1.clear();
  muon_etaAtSt2.clear();
  muon_phi.clear();
  muon_phiAtSt1.clear();
  muon_phiAtSt2.clear();
  muon_pt.clear();
  muon_ptCorr.clear();
  muon_charge.clear();

  muon_dz.clear();
  muon_dzError.clear();
  muon_dxy.clear();
  muon_dxyError.clear();
  muon_3dIP.clear();
  muon_3dIPError.clear();

  muon_PassTightID.clear();
  muon_PassLooseID.clear();
  muon_isSAMuon.clear() ;

  muon_size = 0;

  standAloneMuon_eta.clear();
  standAloneMuon_etaAtSt1.clear();
  standAloneMuon_etaAtSt2.clear();
  standAloneMuon_phi.clear();
  standAloneMuon_phiAtSt1.clear();
  standAloneMuon_phiAtSt2.clear();
  standAloneMuon_pt.clear();
  standAloneMuon_ptCorr.clear();
  standAloneMuon_charge.clear();
  standAloneMuon_dz.clear();
  standAloneMuon_dzError.clear();
  standAloneMuon_dxy.clear();
  standAloneMuon_dxyError.clear();
  standAloneMuon_3dIP.clear();
  standAloneMuon_3dIPError.clear();
  standAloneMuon_PassTightID.clear();
  standAloneMuon_PassLooseID.clear();
  standAloneMuon_isSAMuon.clear() ;
  standAloneMuon_size = 0;

  cosmicMuon_eta.clear();
  cosmicMuon_etaAtSt1.clear();
  cosmicMuon_etaAtSt2.clear();
  cosmicMuon_phi.clear();
  cosmicMuon_phiAtSt1.clear();
  cosmicMuon_phiAtSt2.clear();
  cosmicMuon_pt.clear();
  cosmicMuon_ptCorr.clear();
  cosmicMuon_charge.clear();
  cosmicMuon_dz.clear();
  cosmicMuon_dzError.clear();
  cosmicMuon_dxy.clear();
  cosmicMuon_dxyError.clear();
  cosmicMuon_3dIP.clear();
  cosmicMuon_3dIPError.clear();
  cosmicMuon_PassTightID.clear();
  cosmicMuon_PassLooseID.clear();
  cosmicMuon_isSAMuon.clear() ;
  cosmicMuon_size = 0;

  l1mu_qual.clear();
  l1mu_charge.clear();
  l1mu_pt.clear();
  l1mu_pt_dxy.clear();
  l1mu_dxy.clear();
  l1mu_eta.clear();
  l1mu_etaAtVtx.clear();
  l1mu_phi.clear();
  l1mu_phiAtVtx.clear();
  l1mu_tfIdx.clear();
  l1mu_bx.clear();
  l1mu_size = 0;

  BMTFMu_processor.clear();
  BMTFMu_hwPt.clear();
  BMTFMu_hwQual.clear();
  BMTFMu_hwSign.clear();
  BMTFMu_hwSignValid.clear();
  BMTFMu_hwEta.clear();
  BMTFMu_hwPhi.clear();

}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void L1TMuonAODAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(L1TMuonAODAnalyzer);
