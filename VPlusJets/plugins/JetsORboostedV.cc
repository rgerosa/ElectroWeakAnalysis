// ########################################                                                                                                                                             
// ### -*- C++ -*-                                                                                                                                                                       
// ### Class:     JetsORboostedV                                                                                                                                                  
// ### class JetsORboostedV JetsORboostedV.cc ElectroWeakAnalysis/VPlusJets/plugins/JetsORboostedV.cc                                                                  
// ### Original Author:  Raffaele Angelo Gerosa                                                                                                                                          
// ### Created:  Tue Oct 08 11:57:22 CEST 2013                                                                                                                                           
// ######################################## 

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Photon.h"


class JetsORboostedV : public edm::EDFilter {

   public:
      explicit JetsORboostedV(const edm::ParameterSet&);
      ~JetsORboostedV();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:

      virtual void beginJob() ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual bool beginRun(edm::Run&, edm::EventSetup const&);
      virtual bool endRun(edm::Run&, edm::EventSetup const&);
      virtual bool beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual bool endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

      // ----------member data ---------------------------

      edm::InputTag mInputBoson_;
      edm::InputTag mInputJets_;
      edm::InputTag mInputPhotons_;

      unsigned int minNumberJets_;
      unsigned int maxNumberJets_;
      double minVPt_;
      unsigned int minNumberPhotons_;
};

JetsORboostedV::JetsORboostedV(const edm::ParameterSet& iConfig){

  if(iConfig.existsAs<edm::InputTag>("srcVectorBoson"))
    mInputBoson_ = iConfig.getParameter<edm::InputTag>("srcVectorBoson") ; 
  else mInputBoson_ = edm::InputTag("bestWmunu") ;

  if(iConfig.existsAs<edm::InputTag>("srcJets"))
    mInputJets_  = iConfig.getParameter<edm::InputTag>("srcJets");
  else mInputJets_ = edm::InputTag("ak5PFJetsPtSkimmed") ;

  if(iConfig.existsAs<edm::InputTag>("srcPhotons"))
    mInputPhotons_  = iConfig.getParameter<edm::InputTag>("srcPhotons");
  else mInputPhotons_ = edm::InputTag("selectedPatPhotons") ;

  if(iConfig.existsAs<int>("minNumber"))
    minNumberJets_ = iConfig.getUntrackedParameter<int>("minNumber");
  else minNumberJets_ = 2 ;

  if(iConfig.existsAs<int>("maxNumberJets"))
    maxNumberJets_ = iConfig.getUntrackedParameter<int>("maxNumberJets");
  else maxNumberJets_ = 100 ;

  if(iConfig.existsAs<double>("minVPt"))
    minVPt_ = iConfig.getUntrackedParameter<double>("minVPt");
  else minVPt_ = 100. ;

  if(iConfig.existsAs<int>("minNumberPhotons"))
    minNumberPhotons_ = iConfig.getUntrackedParameter<int>("minNumberPhotons");
  else minNumberPhotons_ = 0 ;
 
}


JetsORboostedV::~JetsORboostedV(){}



// ------------ method called on each new Event  ------------
bool JetsORboostedV::filter(edm::Event& iEvent, const edm::EventSetup& iSetup){

  bool result = false;

  edm::Handle <reco::CandidateView> boson;
  iEvent.getByLabel( mInputBoson_, boson);
  if( boson->size()!=1 ) return false; // Nothing to analyze ...

  const reco::Candidate *Vboson = &((*boson)[0]); 
  if( Vboson == 0) return false;

  if(Vboson->pt()>minVPt_) result = true;


  edm::Handle<edm::View<reco::Jet> > jets;
  iEvent.getByLabel( mInputJets_, jets );

  if(jets->size() >= minNumberJets_ && jets->size() <=maxNumberJets_) result = true;
  
  edm::Handle<reco::CandidateView> photonH;
  iEvent.getByLabel(mInputPhotons_,photonH);

  if(photonH->size() < minNumberPhotons_ ) result = false;
  return result;
}

// ------------ method called once each job just before starting event loop  ------------
void JetsORboostedV::beginJob(){}

// ------------ method called once each job just after ending the event loop  ------------
void JetsORboostedV::endJob() {}

// ------------ method called when starting to processes a run  ------------
bool JetsORboostedV::beginRun(edm::Run&, edm::EventSetup const&){ return true; }

// ------------ method called when ending the processing of a run  ------------
bool JetsORboostedV::endRun(edm::Run&, edm::EventSetup const&){ return true; }

// ------------ method called when starting to processes a luminosity block  ------------
bool JetsORboostedV::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&){ return true; }

// ------------ method called when ending the processing of a luminosity block  ------------
bool JetsORboostedV::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&){ return true; }

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void JetsORboostedV::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(JetsORboostedV);
