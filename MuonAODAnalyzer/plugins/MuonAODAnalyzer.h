// -*- C++ -*-
//
// Package:    MuonAODAnalyzer/MuonAODAnalyzer
// Class:      MuonAODAnalyzer
//
/**\class MuonAODAnalyzer MuonAODAnalyzer.cc MuonAODAnalyzer/MuonAODAnalyzer/plugins/MuonAODAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Efe Yigitbasi
//         Created:  Sat, 10 Sep 2022 11:08:53 GMT
//
//

// system include files
#include <memory>
#include <iostream>
#include <fstream>
#include <string>
#include "TLorentzVector.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TTree.h"
#include "TMath.h"
#include <fmt/printf.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CondFormats/DataRecord/interface/L1TUtmTriggerMenuRcd.h"
#include "CondFormats/L1TObjects/interface/L1TUtmTriggerMenu.h"


// Data formats
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/L1Trigger/interface/Muon.h"
#include "DataFormats/L1Trigger/interface/BXVector.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/L1TGlobal/interface/GlobalAlgBlk.h"
#include "DataFormats/L1TGlobal/interface/GlobalExtBlk.h"
#include "DataFormats/L1TMuon/interface/RegionalMuonCand.h"

// Information about stations
#include "DataFormats/MuonReco/interface/Muon.h"

// muon track extrapolation
#include "MuonAnalysis/MuonAssociators/interface/PropagateToMuonSetup.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"


//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.
using namespace edm;
// using namespace std;
using namespace reco;


class MuonAODAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
  public:
    explicit MuonAODAnalyzer(const edm::ParameterSet&);
    ~MuonAODAnalyzer() override;

  private:
    void beginJob() override;
    void analyze(const edm::Event&, const edm::EventSetup&) override;
    void endJob() override;
    // void beginRun(const edm::Run&, const edm::EventSetup&);
    // void endRun(const edm::Run&, const edm::EventSetup&);
    virtual void InitandClearStuff();

    // void fillTree();
    // void makeTree();


    // ----------member data ---------------------------
    edm::EDGetTokenT<std::vector< reco::Muon> > muonToken_;
    edm::EDGetTokenT<l1t::MuonBxCollection>l1MuonToken_;
    edm::EDGetTokenT<BXVector<l1t::RegionalMuonCand>> l1BMTFRegionalMuonCandToken_;
    edm::EDGetTokenT<std::vector<Vertex> > verticesToken_;
    edm::EDGetTokenT<edm::TriggerResults> trgresultsToken_;
    edm::EDGetTokenT<GlobalExtBlkBxCollection> UnprefirableEventToken_;
    edm::EDGetTokenT<BXVector<GlobalAlgBlk>> l1GtToken_;

    edm::EDGetTokenT<edm::TriggerResults> triggerResultsToken_;
    edm::EDGetTokenT<trigger::TriggerEvent> triggerSummaryToken_;


    Float_t MuonPtCut_;
    Bool_t SaveTree_, IsMC_, Debug_;

    const PropagateToMuonSetup muPropagatorSetup1st_;
    const PropagateToMuonSetup muPropagatorSetup2nd_;

    PropagateToMuon muPropagator1st_;
    PropagateToMuon muPropagator2nd_;

    TTree* outputTree;

    unsigned long _eventNb;
    unsigned long _runNb;
    unsigned long _lumiBlock;
    unsigned long _bx;

    //Nb of primary vertices
    int _n_PV;
    // Float_t _LV_x,_LV_y,_LV_z;
    // Float_t _LV_errx,_LV_erry,_LV_errz;
    // Float_t _PUV1_x,_PUV1_y,_PUV1_z;
    int trueNVtx;

    GlobalAlgBlk const *results_;
    unsigned long long cache_id_;

    //MINIAOD original MET filters decisions
    bool Flag_goodVertices;
    bool Flag_globalTightHalo2016Filter;
    bool Flag_globalSuperTightHalo2016Filter;
    bool Flag_BadPFMuonFilter;
    bool Flag_BadPFMuonDzFilter;

    //Muons
    vector<Float_t>  muon_mass;
    vector<Float_t>  muon_eta;
    vector<Float_t>  muon_etaAtSt1;
    vector<Float_t>  muon_etaAtSt2;
    vector<Float_t>  muon_phi;
    vector<Float_t>  muon_phiAtSt1;
    vector<Float_t>  muon_phiAtSt2;
    vector<Float_t>  muon_pt;
    vector<Float_t>  muon_ptCorr;
    vector <int>     muon_charge;
    vector<Float_t>  muon_dz;
    vector<Float_t>  muon_dzError;
    vector<Float_t>  muon_dxy;
    vector<Float_t>  muon_dxyError;
    vector<Float_t>  muon_3dIP;
    vector<Float_t>  muon_3dIPError;
    vector<Bool_t>  muon_PassTightID;
    vector<Bool_t>  muon_PassLooseID;
    vector<Bool_t> muon_isSAMuon;
    vector<Bool_t> muon_isGlobalMuon;
    vector<Bool_t> muon_isTrackerMuon;
    vector<Bool_t> muon_isPFMuon;

    vector<Float_t> muon_vx;
    vector<Float_t> muon_vy;
    vector<Float_t> muon_vz;
    vector<Float_t> muon_px;
    vector<Float_t> muon_py;
    vector<Float_t> muon_pz;

    vector<int> muon_nChambers;
    vector<int> muon_nChambersCSCorDT;
    vector<int> muon_nMatches;
    vector<int> muon_nMatchedStations;
    vector<unsigned int> muon_expectedNumberOfMatchedStations;
    vector<unsigned int> muon_stationMask;
    vector<int> muon_nMatchedRPCLayers;
    vector<unsigned int> muon_RPClayerMask;
    int muon_size;

    //L1 muon
    vector <int> l1mu_qual;
    vector <int> l1mu_charge;
    vector <Float_t> l1mu_pt;
    vector <Float_t> l1mu_pt_dxy;
    vector <int> l1mu_dxy;
    vector <Float_t> l1mu_eta;
    vector <Float_t> l1mu_etaAtVtx;
    vector <Float_t> l1mu_phi;
    vector <Float_t> l1mu_phiAtVtx;
    vector <int> l1mu_tfIdx;
    vector <int> l1mu_bx;
    int l1mu_size;
    
    // BMTF RegionalMuonCand muons
    vector <Float_t> BMTFMu_processor;
    vector <Float_t> BMTFMu_hwPt;
    vector <Float_t> BMTFMu_hwQual;
    vector <Float_t> BMTFMu_hwSign;
    vector <Float_t> BMTFMu_hwSignValid;
    vector <Float_t> BMTFMu_hwEta;
    vector <Float_t> BMTFMu_hwPhi;

    //Triggers
    bool HLT_IsoMu27;
    bool HLT_IsoMu24;
  
    bool Flag_IsUnprefirable;
    bool passL1_Final_bxmin1;
    bool passL1_Final_bxmin2;

    vector<float> TrigObj_eta;
    vector<float> TrigObj_phi;
    vector<int>   TrigObj_id;
    vector<int>   TrigObj_filterBits;


    //
    // constants, enums and typedefs
    //

    const int  N_METFilters=18;
    enum METFilterIndex{
      idx_Flag_goodVertices,
      idx_Flag_globalTightHalo2016Filter,
      idx_Flag_globalSuperTightHalo2016Filter,
      idx_Flag_HBHENoiseFilter,
      idx_Flag_HBHENoiseIsoFilter,
      idx_Flag_EcalDeadCellTriggerPrimitiveFilter,
      idx_Flag_BadPFMuonFilter,
      idx_Flag_BadPFMuonDzFilter,
      idx_Flag_hfNoisyHitsFilter,
      idx_Flag_BadChargedCandidateFilter,
      idx_Flag_eeBadScFilter,
      idx_Flag_ecalBadCalibFilter,
      idx_Flag_ecalLaserCorrFilter,
      idx_Flag_EcalDeadCellBoundaryEnergyFilter,
      idx_PassecalBadCalibFilter_Update,
      idx_PassecalLaserCorrFilter_Update,
      idx_PassEcalDeadCellBoundaryEnergyFilter_Update,
      idx_PassBadChargedCandidateFilter_Update
    };


    //
    // static data member definitions
    //

};