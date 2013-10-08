// ########################################                                                                                                                                             
// ### -*- C++ -*-                                                                                                                                                                     
// ### Class:    JetCleaner                                                                                                                                                       
// ### class JetCleaner JetCleaner.cc ElectroWeakAnalysis/VPlusJets/plugins/JetCleaner.cc                                                                               
// ### Original Author:  Raffaele Angelo Gerosa                                                                                                                                      
// ### Created:  Tue Oct 08 11:57:22 CEST 2013                                                                                                                         
// ######################################## 

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
 
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/JetReco/interface/JetCollection.h"
#include "DataFormats/JetReco/interface/CaloJet.h"  
#include "DataFormats/JetReco/interface/CaloJetCollection.h"  
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/JPTJet.h"
#include "DataFormats/JetReco/interface/JPTJetCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/Math/interface/deltaR.h"

#include <memory>
#include <vector>
#include <sstream>


template<typename T>
class JetCleaner : public edm::EDProducer{
 
 public:

   // construction/destruction
   explicit JetCleaner(const edm::ParameterSet& iConfig);
   virtual ~JetCleaner();
   
   // member functions
   void produce(edm::Event& iEvent,const edm::EventSetup& iSetup);
   void endJob();

 private:  
   // member data
   edm::InputTag              srcJets_;
   std::vector<edm::InputTag> srcObjects_;
   std::string  moduleLabel_;
   double                     deltaRMin_;

   int idLevel_;
   double etaMax_;
   double etaMin_;
   double ptMin_;

   unsigned int nJetsTot_;
   unsigned int nJetsClean_;

   bool debug_ ;
};



template<typename T>
JetCleaner<T>::JetCleaner(const edm::ParameterSet& iConfig){

  if(iConfig.existsAs<edm::InputTag>("srcJets"))
    srcJets_ = iConfig.getParameter<edm::InputTag>("srcJets");
  else srcJets_ = edm::InputTag("selectedPatJetsPFlow");
 
  if(iConfig.existsAs<std::vector<edm::InputTag> >("srcObjects") )
     srcObjects_ = iConfig.getParameter<std::vector<edm::InputTag> > ("srcObjects");
  else{ 
       srcObjects_.push_back(edm::InputTag("looseElectrons"));
       srcObjects_.push_back(edm::InputTag("looseMuons"));
  }

  if(iConfig.existsAs<double>("deltaRMin") )
    deltaRMin_ =  iConfig.getParameter<double> ("deltaRMin");
  else 
    deltaRMin_ = 0.3 ;

  if(iConfig.existsAs<std::string>("@module_label") )
    moduleLabel_ = iConfig.getParameter<std::string> ("@module_label");
  else
    moduleLabel_ = "";

  if(iConfig.existsAs<int>("idLevel"))
    idLevel_ = iConfig.getParameter<int> ("idLevel");   
  else idLevel_ = 0;

  if(iConfig.existsAs<double>("etaMax"))
    etaMax_ = iConfig.getParameter<double> ("etaMax");   
  else etaMax_ = 10.;

  if(iConfig.existsAs<double>("etaMin"))
    etaMin_ = iConfig.getParameter<double> ("etaMin");   
  else etaMin_ = 0.;

  if(iConfig.existsAs<double>("ptMin"))
    ptMin_ = iConfig.getParameter<double> ("ptMin");
  else ptMin_ = 0.;

  if(iConfig.existsAs<bool>("debug"))
    debug_ = iConfig.getParameter<bool> ("debug");
  else debug_ = false;


  nJetsTot_ = 0 ;
  nJetsClean_ = 0 ;

  produces<std::vector<T> > ();

}


//______________________________________________________________________________
template<typename T>
JetCleaner<T>::~JetCleaner(){}


//______________________________________________________________________________
template<typename T>
void JetCleaner<T>::produce(edm::Event& iEvent,const edm::EventSetup& iSetup){

  std::auto_ptr<std::vector<T> > cleanJets(new std::vector<T>);
  
  // take the jet collection
  edm::Handle<reco::JetView> jets;
  iEvent.getByLabel(srcJets_,jets);

  bool* isClean = new bool[jets->size()];

  for (unsigned int iJet=0;iJet<jets->size();iJet++) isClean[iJet] = true;
  
  for (unsigned int iSrc=0;iSrc<srcObjects_.size();iSrc++) {
 
    // get the reco candidate to be used for cleaning
    edm::Handle<reco::CandidateView> objects;
    iEvent.getByLabel(srcObjects_[iSrc],objects);
    
    for (unsigned int iJet=0;iJet<jets->size();iJet++) {
      const reco::Jet& jet = jets->at(iJet);
      for (unsigned int iObj=0;iObj<objects->size();iObj++) {
	const reco::Candidate& obj = objects->at(iObj);
	double deltaR = reco::deltaR(jet,obj);
	if (deltaR<deltaRMin_)  isClean[iJet] = false;
      }
    }
  }
  
  for (unsigned int iJet=0;iJet<jets->size();iJet++)
    if (isClean[iJet]) {
      
      bool passedId=false;
      bool ThisIsClean=false;

      const std::type_info & type = typeid((*jets)[iJet]); 
      if( type == typeid(reco::CaloJet) ) {
	passedId = true;
      }
      //calculate the PF jetID
      else if ( type == typeid(reco::PFJet) ) {
	const reco::PFJet pfjet = static_cast<const reco::PFJet &>((*jets)[iJet]);
	ThisIsClean=true;
	//apply following only if |eta|<2.4: CHF>0, CEMF<0.99, chargedMultiplicity>0   
	if(( pfjet.chargedHadronEnergy()/ pfjet.energy())<= 0.0  
	   && fabs(pfjet.eta())<2.4) ThisIsClean=false; 
	if( (pfjet.chargedEmEnergy()/pfjet.energy())>= 0.99 
	    && fabs(pfjet.eta())<2.4 ) ThisIsClean=false;
	if( pfjet.chargedMultiplicity()<=0 && fabs(pfjet.eta())<2.4 ) 
	  ThisIsClean=false;	
	// always require #Constituents > 1
	if( pfjet.nConstituents() <=1 ) ThisIsClean=false;
	if(ThisIsClean && 
	   (pfjet.neutralHadronEnergy()/pfjet.energy())< 0.99 
	   && (pfjet.neutralEmEnergy()/pfjet.energy())<0.99) 
	  passedId=true;	
      }
      // in case of GenJet apply no jet ID
      else if ( type == typeid(reco::GenJet) ) passedId=true;

      else if ( type == typeid(pat::Jet) ){
	const pat::Jet patjet = static_cast<const pat::Jet &>((*jets)[iJet]);
	ThisIsClean=true;
	//apply following only if |eta|<2.4: CHF>0, CEMF<0.99, chargedMultiplicity>0   
	if(( patjet.chargedHadronEnergy()/ patjet.energy())<= 0.0  
	   && fabs(patjet.eta())<2.4) ThisIsClean=false; 
	if( (patjet.chargedEmEnergy()/patjet.energy())>= 0.99 
	    && fabs(patjet.eta())<2.4 ) ThisIsClean=false;
	if( patjet.chargedMultiplicity()<=0 && fabs(patjet.eta())<2.4 ) 
	  ThisIsClean=false;	
	// always require #Constituents > 1
	if( patjet.nConstituents() <=1 ) ThisIsClean=false;
	if(ThisIsClean && 
	   (patjet.neutralHadronEnergy()/patjet.energy())< 0.99 
	   && (patjet.neutralEmEnergy()/patjet.energy())<0.99) 
	  passedId=true;	

      }
      
      bool isPassing = false;
      if(idLevel_==0) isPassing = true; // only cleaning in dR
      if(idLevel_==1) isPassing = passedId; // also ID loose

      const T& goodJet = static_cast<const T&>((*jets)[iJet]);
      double pt = goodJet.pt();
      double eta = goodJet.eta();
      if(isPassing && fabs(eta)<etaMax_ && fabs(eta)>=etaMin_ && pt>ptMin_) cleanJets->push_back( goodJet );
    }

  nJetsTot_  +=jets->size();
  nJetsClean_+=cleanJets->size();

  delete [] isClean;  
  iEvent.put(cleanJets);
}




//______________________________________________________________________________
template<typename T>
void JetCleaner<T>::endJob(){

  if(debug_){
   std::stringstream ss;
   std::cout<<"##################"<<std::endl;
   std::cout<<"### JetCleaner ###"<<std::endl;
   std::cout<<"##################"<<std::endl;
   ss<<"nJetsTot="<<nJetsTot_<<" nJetsClean="<<nJetsClean_<<" fJetsClean="<<100.*(nJetsClean_/(double)nJetsTot_)<<"%\n";
   std::cout<<ss<<std::endl;
   std::cout<<"##################"<<std::endl;
  }

}


////////////////////////////////////////////////////////////////////////////////
// plugin definition
////////////////////////////////////////////////////////////////////////////////

typedef JetCleaner<reco::CaloJet> CaloJetCleaner;
typedef JetCleaner<reco::PFJet>   PFJetCleaner;
typedef JetCleaner<reco::JPTJet>  JPTJetCleaner;
typedef JetCleaner<reco::GenJet>  GenJetCleaner;
typedef JetCleaner<pat::Jet>      PFPATJetCleaner;

DEFINE_FWK_MODULE(CaloJetCleaner);
DEFINE_FWK_MODULE(PFJetCleaner);
DEFINE_FWK_MODULE(JPTJetCleaner);
DEFINE_FWK_MODULE(GenJetCleaner);
DEFINE_FWK_MODULE(PFPATJetCleaner);
