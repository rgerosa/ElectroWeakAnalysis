// ########################################                                                                                                  
// ### -*- C++ -*-                                                                                                                                             
// ### Class:      ElectronIdSelector                                                                                                                        
// ### class ElectronIdSelector ElectronIdSelector.cc ElectroWeakAnalysis/VPlusJets/plugins/ElectronIdSelector.cc                                                       
// ### Original Author:  Raffaele Angelo Gerosa                                                                                                                                
// ### Created:  Tue Oct 08 11:57:22 CEST 2013                                                                                                                               
// ######################################## 

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "EgammaAnalysis/ElectronTools/interface/ElectronEffectiveArea.h"

#include <memory>
#include <vector>
#include <sstream>
#include <cmath>


template<typename T>
class ElectronIdSelector : public edm::EDProducer{

 public:
   // construction/destruction
   ElectronIdSelector(const edm::ParameterSet& iConfig);
   virtual ~ElectronIdSelector();
  
   // member functions
   void produce(edm::Event& iEvent,const edm::EventSetup& iSetup);
   void beginJob();
   void endJob();

private:  
   // member data
   edm::InputTag  src_;
   std::string    moduleLabel_;
   std::string    idLabel_;  
   bool           useMVAbasedID_;
   bool           useDetectorIsolation_;
   bool           applyTightID_;
   bool           applyLooseID_;
   bool           applyQcdID_;
   bool           applyHiggsID_;

   unsigned int nTot_;
   unsigned int nPassed_;
};


//______________________________________________________________________________
template<typename T>
ElectronIdSelector<T>::ElectronIdSelector(const edm::ParameterSet& iConfig){

  if(iConfig.existsAs<edm::InputTag>("src"))
    src_ = iConfig.getParameter<edm::InputTag>("src") ;
  else { std::cerr<<" Missing in put collection --> break the producer "<<std::endl; return ; }

  if(iConfig.existsAs<edm::InputTag>("@module_label"))
    moduleLabel_ = iConfig.getParameter<std::string>("@module_label");

  if(iConfig.existsAs<std::string>("idLabel"))
    idLabel_ = iConfig.getParameter<std::string>("idLabel") ;
  else idLabel_ = "loose" ;

  if(iConfig.existsAs<bool>("useMVAbasedID"))
    useMVAbasedID_ = iConfig.getParameter<bool>("useMVAbasedID") ;
  else useMVAbasedID_ = true ;

  if(iConfig.existsAs<bool>("useDetectorIsolation"))
    useDetectorIsolation_ = iConfig.getParameter<bool>("useDetectorIsolation") ;
  else useDetectorIsolation_ = false ;

  nTot_ = 0; nPassed_ = 0;

  produces<std::vector<T> >();

  /// ------- Decode the ID criteria --------
  applyTightID_  = false;
  applyLooseID_  = false;
  applyQcdID_    = false;
  applyHiggsID_  = false;

  if( (idLabel_.compare("tight")==0) || 
      (idLabel_.compare("Tight")==0) || 
      (idLabel_.compare("TIGHT")==0) ||
      (idLabel_.compare("WP70")==0) ||
      (idLabel_.compare("wp70")==0) )  
    applyTightID_ = true;

  else if( (idLabel_.compare("loose")==0) || 
      (idLabel_.compare("Loose")==0) || 
      (idLabel_.compare("LOOSE")==0) ||
      (idLabel_.compare("WP90")==0) ||
      (idLabel_.compare("wp90")==0) )  applyLooseID_ = true;

  else if( (idLabel_.compare("qcd")==0) || 
      (idLabel_.compare("Qcd")==0) || 
      (idLabel_.compare("QCD")==0) ||
      (idLabel_.compare("QCDid")==0) ||
      (idLabel_.compare("QcdId")==0) )  applyQcdID_ = true;

}


//______________________________________________________________________________
template<typename T>
ElectronIdSelector<T>::~ElectronIdSelector(){}



////////////////////////////////////////////////////////////////////////////////
// implementation of member functions
////////////////////////////////////////////////////////////////////////////////

template<typename T>
void ElectronIdSelector<T>::produce(edm::Event& iEvent,const edm::EventSetup& iSetup){


  /////// Pileup density "rho" in the event from fastJet pileup calculation /////
  edm::Handle<double> rho;
  const edm::InputTag eventrho("kt6PFJetsPFlow", "rho");
  iEvent.getByLabel(eventrho,rho);
  double fastJetRho = *rho;


  std::auto_ptr<std::vector<T> > passingElectrons(new std::vector<T >);

  edm::Handle<edm::View<T> > electrons;
  iEvent.getByLabel(src_,electrons);
  
  bool* isPassing = new bool[electrons->size()];
  for(unsigned int iElec=0;iElec<electrons->size();iElec++) { 

    isPassing[iElec]=false;

    const T& ele = electrons->at(iElec);

    // -------- Make sure that the electron is within acceptance ------
    bool isEB = ele.isEB();
    bool isEE = ele.isEE();
    bool inAcceptance = (isEB || isEE) && (ele.ecalDrivenSeed()==1);
    float pt  = ele.p4().Pt();


    // -------- Compute isolation ------
    const double PI = 4.0*atan(1.0);

    float detector_isolation = (ele.dr03TkSumPt() + std::max(0.,ele.dr03EcalRecHitSumEt()-1.0) + ele.dr03HcalTowerSumEt() - PI*0.3*0.3*fastJetRho) / pt;

    float EffArea = ElectronEffectiveArea::GetElectronEffectiveArea( ElectronEffectiveArea::kEleGammaAndNeutralHadronIso03 , ele.eta() , ElectronEffectiveArea::kEleEAData2012);

    float pfIso03EA = (ele.pfIsolationVariables().chargedHadronIso + max(0.,ele.pfIsolationVariables().neutralHadronIso+ele.pfIsolationVariables().photonIso-EffArea*fastJetRho))/pt;
 
    float isolation = 0.;
    if(useDetectorIsolation_) isolation = detector_isolation;
    else isolation = pfIso03EA;

    // -------- Compute ID ------
    Double_t sihih   = ele.sigmaIetaIeta();
    Double_t Dphi    = fabs(ele.deltaPhiSuperClusterTrackAtVtx());
    Double_t Deta    = fabs(ele.deltaEtaSuperClusterTrackAtVtx());
    Double_t abseta  = fabs(ele.eta());
    Int_t innerHits  = ele.gsfTrack()->trackerExpectedHitsInner().numberOfHits();
     
    Double_t dist = ele.convDist(); // default value is -9999 if conversion partner not found
    Double_t dcot = ele.convDcot(); // default value is -9999 if conversion partner not found
    Bool_t isConv = fabs(dist) < 0.02 && fabs(dcot) < 0.02;

    Double_t mvaTrigV0 = -2;
    Double_t mvaNonTrigV0 = -2;

    const pat::Electron* patele = dynamic_cast<const pat::Electron *>( &ele );
    if(patele->isElectronIDAvailable("mvaTrigV0")) {
           mvaTrigV0    =  patele->electronID("mvaTrigV0");
           mvaNonTrigV0 =  patele->electronID("mvaNonTrigV0");
           isConv = ! ( patele->passConversionVeto() );
    }

    bool isTight = false;  /////// <--- equivalent to WP70
    bool isLoose = false;  /////// <--- equivalent to WP90
    bool isQCD   = false;  /////// <--- invert the isolation cut

    // -------- if cut-based ID -----------------
    if( !useMVAbasedID_) {

       isTight = (pt>30.) && inAcceptance && (innerHits==0) && 
	 (isolation<0.05) && (!isConv) && 
	 ((isEB && sihih<0.01 && Dphi<0.03 && Deta<0.004) || 
	 (isEE && sihih<0.03 && Dphi<0.02 && Deta<0.005));       

       isLoose = (pt>20.) && inAcceptance && (innerHits<=1) && 
	(isolation<0.1) && (!isConv) && 
	((isEB && sihih<0.01 && Dphi<0.8 && Deta<0.007) || 
	 (isEE && sihih<0.03 && Dphi<0.07 && Deta<0.005)); 

       isQCD = (pt>20.) && inAcceptance && (innerHits<=1) && 
	 (isolation>0.1) && (!isConv) && 
	 ((isEB && sihih<0.01 && Dphi<0.8 && Deta<0.007) || 
	  (isEE && sihih<0.03 && Dphi<0.07 && Deta<0.005)); 

     } else {

      //-------- if MVA-based ID -----------------
     // WP 80 pT>20 GeV	0.913 / 0.105 	0.964 / 0.178 	0.899 / 0.150

      isTight = (pt>30.) && inAcceptance && (!isConv) && ( 
        ( abseta <= 0.8 && mvaTrigV0 > 0.913&& pfIso03EA < 0.105) ||
        ( abseta > 0.8 && abseta <= 1.479 && mvaTrigV0 > 0.964 && pfIso03EA < 0.178) ||
        ( abseta > 1.479 && mvaTrigV0 > 0.899 && pfIso03EA < 0.150) );

    // WP 90 	0.877 / 0.426 	0.811 / 0.481 	0.707 / 0.390 

      isLoose = (pt>20.) && inAcceptance && (!isConv) && ( 
        ( abseta <= 0.8 && mvaNonTrigV0 > 0.877 && pfIso03EA < 0.426 ) ||
        ( abseta > 0.8 && abseta <= 1.479 && mvaNonTrigV0 > 0.811 && pfIso03EA <0.481 ) ||
        ( abseta > 1.479 && mvaNonTrigV0 > 0.707 && pfIso03EA < 0.390) );

      isQCD = (pt>20.) && inAcceptance && (!isConv) && ( 
        ( abseta <= 0.8 && mvaTrigV0 > -1. && pfIso03EA > 0.177 ) ||
        ( abseta > 0.8 && abseta <= 1.479 && mvaTrigV0 > -1. && pfIso03EA >0.180 ) ||
        ( abseta > 1.479 && mvaTrigV0 > -1. && pfIso03EA > 0.244) );
     }

     //D/ ------- Finally apply selection --------
     if(applyTightID_ && isTight) isPassing[iElec]= true;
     if(applyLooseID_ && isLoose) isPassing[iElec]= true;
     if(applyQcdID_ && isQCD) isPassing[iElec]= true;
  
 } 

 unsigned int counter=0;
 typename edm::View<T>::const_iterator tIt, endcands = electrons->end();
 for (tIt = electrons->begin(); tIt != endcands; ++tIt, ++counter) {
 if(isPassing[counter]) passingElectrons->push_back( *tIt );  
 }

 nTot_  +=electrons->size();
 nPassed_+=passingElectrons->size();

 delete [] isPassing;  
 iEvent.put(passingElectrons);
}

// ____________________________________________________________________________

template<typename T>
void ElectronIdSelector<T>::beginJob(){}


//______________________________________________________________________________
template<typename T>
void ElectronIdSelector<T>::endJob()
{
  std::stringstream ss;
  ss<<"nTot="<<nTot_<<" nPassed="<<nPassed_
    <<" effPassed="<<100.*(nPassed_/(double)nTot_)<<"%\n";

}


////////////////////////////////////////////////////////////////////////////////
// plugin definition
////////////////////////////////////////////////////////////////////////////////
typedef ElectronIdSelector<reco::GsfElectron>     GsfElectronIdSelector;
typedef ElectronIdSelector<pat::Electron>         PATElectronIdSelector;

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(GsfElectronIdSelector);
DEFINE_FWK_MODULE(PATElectronIdSelector);
