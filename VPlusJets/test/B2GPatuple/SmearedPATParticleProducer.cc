#ifndef PhysicsTools_PatUtils_SmearedPATParticleProducerT_h
#define PhysicsTools_PatUtils_SmearedPATParticleProducerT_h
/** \class SmearedParticleProducer
*
* smear energy of electrons/muons/tau-jets by +/- 1 standard deviation, 
* in order to estimate resulting uncertainty on MET
*
* NOTE: energy resolution need to be specified in python config
*
* \author Raffaele Angelo Gerosa, INFN MIB
*
* \version $Revision: 1.1 $
*
* $Id: SmearedParticleProducer.h,v 1.2 2013/09/30 15:49:31 Gerosa Exp $
*
*/                                                                                                             

////////////////////////////////////////////////////////////////////////////////                                                                                             
// Includes                                                                                 
//////////////////////////////////////////////////////////////////////////////

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Common/interface/View.h"

#include <memory>
#include <vector>
#include <sstream>
#include <cmath>

#include "TRandom3.h"

////////////////////////////////////////////////////////////////////////////////                                                                                               
// class definition                                                                                                                                                              
////////////////////////////////////////////////////////////////////////////////                                                                                                        

template<typename T>
class SmearedPATParticleProducer : public edm::EDProducer  {

  typedef std::vector<T> ParticleCollection;
 
  public:

   // construction/destruction
   SmearedPATParticleProducer(const edm::ParameterSet& iConfig);
   virtual ~SmearedPATParticleProducer();

   // member functions
   void produce(edm::Event& iEvent,const edm::EventSetup& iSetup);
   void endJob();

  private:

   // member data                                                                                                                                                                         
   edm::InputTag  src_;
   std::string    moduleLabel_;
   double         shiftBy_ ;
   double         smearBy_ ;

   TRandom3       *rand_;

   struct binningEntryType {
    
     binningEntryType(double Resolution, double ResolutionUncertainty):
                                             binSelection_(0),
	                                     binResolution_(Resolution),
	                                     binResolutionUncertainty_(ResolutionUncertainty){}

       binningEntryType(const edm::ParameterSet& cfg): 
                binSelection_(new StringCutObjectSelector<T>(cfg.getParameter<std::string>("binSelection"))),
		binResolution_(cfg.getParameter<double>("binResolution")),
		binResolutionUncertainty_(cfg.getParameter<double>("binResolutionUncertainty")){}

       ~binningEntryType() {delete binSelection_;}

       StringCutObjectSelector<T>* binSelection_;
       double binResolution_;
       double binResolutionUncertainty_;

   };

   std::vector<binningEntryType*> binning_;

};
#endif

////////////////////////////////////////////////////////////////////////////////                                                                                                     
// construction/destruction                                                                                                                                                           
///////////////////////////////////////////////////////////////////////////////

template<typename T>
SmearedPATParticleProducer<T>::SmearedPATParticleProducer(const edm::ParameterSet& iConfig){

  moduleLabel_ = iConfig.getParameter<std::string>("@module_label");
  src_         = iConfig.getParameter<edm::InputTag>("src");

  if(iConfig.exists("shiftBy")) shiftBy_ = iConfig.getParameter<double>("shiftBy");
  else shiftBy_ = 0. ;

  if(iConfig.exists("smearBy")) shiftBy_ = iConfig.getParameter<double>("smearBy");
  else smearBy_ = 1. ;

  if(iConfig.exists("binning")){
    std::vector<edm::ParameterSet> cfgBinning =  iConfig.getParameter<std::vector<edm::ParameterSet> >("binning");
    for(std::vector<edm::ParameterSet>::const_iterator iBinning = cfgBinning.begin(); iBinning!=cfgBinning.end() ; iBinning++)
      binning_.push_back(new binningEntryType(*iBinning));   
  }
  else if(iConfig.exists("resolution") && iConfig.exists("resolutionUnc")){             
    binning_.push_back(new binningEntryType(iConfig.getParameter<double>("resolution"),iConfig.getParameter<double>("resolutionUnc")));
  }
  else{ 
        std::cerr << "Invalid binning definition : you must specify or a VPset: binDefinition, Resolution, ResolutionUncertainty or"
         " .\n Resolution and ResolutionUncertainty ascms.double --> module disabled " << ".\n";
  }
   
  produces<std::vector<T> >();

  rand_ = new TRandom3();
  rand_->SetSeed(); 
  
}
 
template<typename T>
SmearedPATParticleProducer<T>::~SmearedPATParticleProducer(){

  for ( typename std::vector<binningEntryType*>::const_iterator it = binning_.begin();it != binning_.end(); ++it ) delete (*it);
  if(!rand_) delete rand_;
 
}

////////////////////////////////////////////////////////////////////////////////                                                                                                
// implementation of member functions                                                                                                                                             
//////////////////////////////////////////////////////////////////////////////// 


template<typename T>
void SmearedPATParticleProducer<T>::produce(edm::Event& iEvent,const edm::EventSetup& iSetup){

  edm::Handle<edm::View<T> > originalParticles;
  iEvent.getByLabel(src_, originalParticles);

  std::auto_ptr<std::vector<T> > smearedParticles(new std::vector<T>);

  if(binning_.empty()) return ;

  for ( typename edm::View<T>::const_iterator originalParticle = originalParticles->begin();
	originalParticle != originalParticles->end(); ++originalParticle ) {

    double Resolution = 0.;
    double ResolutionUncertainty = 0.;

    for ( typename std::vector<binningEntryType*>::iterator binningEntry = binning_.begin();
          binningEntry != binning_.end(); ++binningEntry ) {
 
     if ( (!(*binningEntry)->binSelection_) || (*(*binningEntry)->binSelection_)(*originalParticle) ) {
          Resolution = (*binningEntry)->binResolution_;
          ResolutionUncertainty = (*binningEntry)->binResolutionUncertainty_;
          break;
     }
   }

   reco::Candidate::LorentzVector smearedParticleP4 = originalParticle->p4();
  
   double smearFactor = 1. ;
   if(smearBy_ > 0. )  smearFactor = smearBy_*Resolution ;
   double smearFactorErr = smearBy_*ResolutionUncertainty ;
   if(shiftBy_ !=0.)  smearFactor = smearFactor+smearFactorErr*shiftBy_ ;
   
   smearedParticleP4 = smearedParticleP4*(1+rand_->Gaus(0.,smearFactor));

   T smearedParticle(*originalParticle); 
     
   smearedParticle.setP4(smearedParticleP4);
   
   smearedParticles->push_back(smearedParticle);
  }
 
  iEvent.put(smearedParticles);

}


template<typename T>
void SmearedPATParticleProducer<T>::endJob(){}


////////////////////////////////////////////////////////////////////////////////                                                                                                       
// plugin definition                                                                                                                                                                 
////////////////////////////////////////////////////////////////////////////////                                                                                                         
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"

typedef SmearedPATParticleProducer<pat::Electron> SmearedPATElectronProducer;
typedef SmearedPATParticleProducer<pat::Photon>   SmearedPATPhotonProducer;
typedef SmearedPATParticleProducer<pat::Muon>     SmearedPATMuonProducer;
typedef SmearedPATParticleProducer<pat::Tau>      SmearedPATTauProducer;

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(SmearedPATElectronProducer);
DEFINE_FWK_MODULE(SmearedPATPhotonProducer);
DEFINE_FWK_MODULE(SmearedPATMuonProducer);
DEFINE_FWK_MODULE(SmearedPATTauProducer);
