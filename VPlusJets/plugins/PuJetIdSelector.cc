// ########################################                                                                                                                                             
// ### -*- C++ -*-                                                                                                                                                                       
// ### Class:    PuJetIdSelector                                                                                                                                                  
// ### class PuJetIdSelector PuJetIdSelector.cc ElectroWeakAnalysis/VPlusJets/plugins/PuJetIdSelector.cc                                                                  
// ### Original Author:  Raffaele Angelo Gerosa                                                                                                                                          
// ### Created:  Tue Oct 08 11:57:22 CEST 2013                                                                                                                                          
// ######################################## 

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "CMGTools/External/interface/PileupJetIdentifier.h"

#include <memory>
#include <vector>
#include <sstream>
#include <cmath>

////////////////////////////////////////////////////////////////////////////////
// class definition
////////////////////////////////////////////////////////////////////////////////
template<typename T>
class PuJetIdSelector : public edm::EDProducer{

  public:
    // construction/destruction
    PuJetIdSelector(const edm::ParameterSet& iConfig);
    virtual ~PuJetIdSelector();
  
    // member functions
    void produce(edm::Event& iEvent,const edm::EventSetup& iSetup);
    void endJob();

  private:  
    // member data
    edm::InputTag  src_;
    std::string    moduleLabel_;
    std::string    idLabel_;  
    std::string    valueMapLabel_;  
    bool           applyTightID_;
    bool           applyMediumID_;
    bool           applyLooseID_;
    bool           applyMVAID_ ;

    unsigned int nTot_;
    unsigned int nPassed_;

    bool debug_ ;
};


//______________________________________________________________________________
template<typename T>
PuJetIdSelector<T>::PuJetIdSelector(const edm::ParameterSet& iConfig){

  if(iConfig.existsAs<edm::InputTag>("src"))
    src_ = iConfig.getParameter<edm::InputTag>("src");
  else src_ = edm::InputTag("selectedPatJetsPFlow");

  if (iConfig.existsAs<std::string>("@module_label"))
    moduleLabel_ = iConfig.getParameter<std::string>   ("@module_label") ;

  if(iConfig.existsAs<std::string>("idLabel"))
    idLabel_ = iConfig.getParameter<std::string>("idLabel") ;
  else idLabel_ = "loose" ;

  if(iConfig.existsAs<std::string>("valueMapLabel"))
    valueMapLabel_ = iConfig.getParameter<std::string>("valueMapLabel");
  else valueMapLabel_ = "puJetMvaChs" ;

  if(iConfig.existsAs<bool>("applyMVAID"))
    applyMVAID_ = iConfig.getParameter<bool>("applyMVAID");
  else valueMapLabel_ = true ;


  if(iConfig.existsAs<bool>("debug"))
    debug_ = iConfig.getParameter<bool>("debug");
  else debug_ = false ;


  nTot_ = 0 ; nPassed_ = 0;

  produces<std::vector<T> >();

  /// ------- Decode the ID criteria --------
  applyTightID_  = false;
  applyMediumID_ = false;
  applyLooseID_  = false;

  if( (idLabel_.compare("tight")==0) || 
      (idLabel_.compare("Tight")==0) || 
      (idLabel_.compare("TIGHT")==0) )  
    applyTightID_ = true;
  
  else if( (idLabel_.compare("medium")==0) || 
           (idLabel_.compare("Medium")==0) || 
           (idLabel_.compare("MEDIUM")==0) )  
    applyMediumID_ = true;
  
  else if( (idLabel_.compare("loose")==0) || 
           (idLabel_.compare("Loose")==0) || 
           (idLabel_.compare("LOOSE")==0) )  
    applyLooseID_ = true;
  
}


//______________________________________________________________________________
template<typename T>
PuJetIdSelector<T>::~PuJetIdSelector(){}



//______________________________________________________________________________
template<typename T>
void PuJetIdSelector<T>::produce(edm::Event& iEvent,const edm::EventSetup& iSetup){

  /////// Pileup jet identification value maps /////
  edm::Handle<edm::ValueMap<float> > puJetIdMVA;
  edm::Handle<edm::ValueMap<int> >   puJetIdFlag;

  if(applyMVAID_){

    iEvent.getByLabel(valueMapLabel_,"fullDiscriminant",puJetIdMVA);
    iEvent.getByLabel(valueMapLabel_,"fullId",puJetIdFlag);
  }
  else{ iEvent.getByLabel(valueMapLabel_,"cutbasedDiscriminant",puJetIdMVA);
        iEvent.getByLabel(valueMapLabel_,"cutbasedId",puJetIdFlag);
  }
   
  edm::Handle<edm::View<T> > jetsHandle;
  iEvent.getByLabel( src_, jetsHandle ); 

  std::auto_ptr<std::vector<T> > passingJets(new std::vector<T >);
  
  bool* isPassing = new bool[jetsHandle->size()];

  typename edm::View<T>::const_iterator itJet = jetsHandle->begin();

  for(unsigned int iJet=0; itJet!=jetsHandle->end(); iJet++,++itJet) { 

    isPassing[iJet]=false;
    
    // Jet Id
    int   idflag = (*puJetIdFlag)[jetsHandle -> refAt(iJet)];

    /// ------- Apply selection --------
    if(applyTightID_) { if(PileupJetIdentifier::passJetId( idflag, PileupJetIdentifier::kTight ) == true )  isPassing[iJet]= true; }
    if(applyMediumID_){ if(PileupJetIdentifier::passJetId( idflag, PileupJetIdentifier::kMedium) == true )  isPassing[iJet]= true; }
    if(applyLooseID_) { if(PileupJetIdentifier::passJetId( idflag, PileupJetIdentifier::kLoose ) == true )  isPassing[iJet]= true; }

    int wp = 0 ;

    if(PileupJetIdentifier::passJetId( idflag, PileupJetIdentifier::kLoose ) == true )  wp=1;
    if(PileupJetIdentifier::passJetId( idflag, PileupJetIdentifier::kMedium) == true )  wp=2;
    if(PileupJetIdentifier::passJetId( idflag, PileupJetIdentifier::kTight ) == true )  wp=3;


    const std::type_info & type = typeid(*itJet);

    if ( type == typeid(pat::Jet) ) {

      pat::Jet pfjet  = static_cast<const pat::Jet &>(*itJet);

      pfjet.addUserInt("PileUpJetIDWorkingPoint",wp);

      passingJets->push_back(pfjet);
    }

  }

   
  unsigned int counter=0;
  typename std::vector<T>::iterator tIt, endcands = passingJets->end();

  for (tIt = passingJets->begin(); tIt != endcands; ++tIt, ++counter) {
    if(!isPassing[counter]) passingJets->erase(tIt); 
  }

  nTot_    += jetsHandle->size();
  nPassed_ += passingJets->size();

  delete [] isPassing;  
  iEvent.put(passingJets);
}


//______________________________________________________________________________
template<typename T>
void PuJetIdSelector<T>::endJob(){

  if(debug_){
   std::stringstream ss;
   std::cout<<"#######################"<<std::endl;
   std::cout<<"### PuJetIdSelector ###"<<std::endl;
   std::cout<<"#######################"<<std::endl;
   ss<<"nTot="<<nTot_<<" nPassed="<<nPassed_<<" effPassed="<<100.*(nPassed_/(double)nTot_)<<"%\n";
   std::cout<<ss<<std::endl;
   std::cout<<"#######################"<<std::endl;
  }

}


////////////////////////////////////////////////////////////////////////////////
// plugin definition
////////////////////////////////////////////////////////////////////////////////
typedef PuJetIdSelector<pat::Jet>         PATPuJetIdSelector;

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PATPuJetIdSelector);

