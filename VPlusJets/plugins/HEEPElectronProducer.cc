// ########################################
// ### -*- C++ -*-                                                                                                            
// ### Class:      HEEPElectronProducer                                                                                                                   
// ### class HEEPElectronProducer HEEPElectronProducer.cc ElectroWeakAnalysis/VPlusJets/plugins/HEEPElectronProducer.cc                                                              
// ### Original Author:  Raffaele Angelo Gerosa                                                                                                                     
// ### Created:  Tue Oct 08 11:57:22 CEST 2013                                                                                                   
// ########################################  

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

//trigger                                                                                                                                                                            
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "EgammaAnalysis/ElectronTools/interface/ElectronEffectiveArea.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"


class HEEPElectronProducer : public edm::EDProducer {
 public:
  explicit HEEPElectronProducer(const edm::ParameterSet & );
  ~HEEPElectronProducer();

 private:
  virtual void beginJob() ;
  virtual void endJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);

  edm::InputTag eleLabel_;
  std::string   IdLabel_;
  double        pTCutValue_;
  bool          applyLooseID_;
  bool          applyTightID_;

  unsigned int nTot_;
  unsigned int nPassed_;

  bool debug_ ;

};

HEEPElectronProducer::HEEPElectronProducer(const edm::ParameterSet &  iConfig){

  if( iConfig.existsAs<edm::InputTag>("electronCollection") )
    eleLabel_=iConfig.getParameter< edm::InputTag >("electronCollection");
  else eleLabel_= edm::InputTag("heepPatElectrons");

  if( iConfig.existsAs<std::string>("eleIdType") )
    IdLabel_=iConfig.getParameter<std::string>("eleIdType");
  else IdLabel_= std::string("TightID");

  if( iConfig.existsAs<double>("pTCutValue") )
    pTCutValue_=iConfig.getParameter<double>("pTCutValue");
  else pTCutValue_ = 0.;

  if( iConfig.existsAs<bool>("debug") )
    debug_=iConfig.getParameter<bool>("debug");
  else debug_ = false;


  applyTightID_ = false ;
  applyLooseID_ = false ;

  if( IdLabel_ == "TightID" || IdLabel_ == "Tight" || IdLabel_ == "tightID" || IdLabel_ == "tightId" || IdLabel_ == "tightid" || IdLabel_ == "tight" ){ 
    applyTightID_ = true ; 
    if(!iConfig.existsAs<double>("pTCutValue") ) pTCutValue_ = 90. ; 
  }

  if( IdLabel_ == "LooseID" || IdLabel_ == "Loose" || IdLabel_ == "LooseID" || IdLabel_ == "looseId" || IdLabel_ == "looseid" || IdLabel_ == "loose" ){
     applyLooseID_ = true ;
     if(!iConfig.existsAs<double>("pTCutValue") ) pTCutValue_ = 20. ;
  }

  produces < pat::ElectronCollection >();

}

HEEPElectronProducer::~HEEPElectronProducer(){}

void HEEPElectronProducer::produce(edm::Event & iEvent, const edm::EventSetup & iSetup){

  edm::Handle<edm::View<pat::Electron> > eleHandle;
  iEvent.getByLabel(eleLabel_,eleHandle);
  const edm::View< pat::Electron>& eles = *(eleHandle.product());

  std::auto_ptr<pat::ElectronCollection> passingElectrons(new pat::ElectronCollection);

  bool* isPassing = new bool[eles.size()];

  for(size_t eleNr = 0; eleNr != eles.size(); ++eleNr) {

   const pat::Electron & ele = eles[eleNr];

   isPassing[eleNr]=false;
 
   int eleCutCode = ele.userInt("HEEPId");
   Double_t et    = ele.caloEnergy()*sin(ele.p4().theta());
   Double_t scEta = fabs(ele.superCluster()->eta());

   if(eleCutCode == 0 && et > pTCutValue_ && fabs(ele.eta()) < 2.5 && !(scEta > 1.4442 && scEta < 1.566) && fabs(ele.phi()) < 3.2 && applyTightID_ ) isPassing[eleNr]=true; 
   if(eleCutCode == 0 && et > pTCutValue_ && fabs(ele.eta()) < 2.5 && !(scEta > 1.4442 && scEta < 1.566) && fabs(ele.phi()) < 3.2 && applyLooseID_ ) isPassing[eleNr]=true;
  }

  unsigned int counter=0;
  edm::View<pat::Electron>::const_iterator tIt, endcands = eles.end();
  for (tIt = eles.begin(); tIt != endcands; ++tIt, ++counter) {
      if(isPassing[counter]) passingElectrons->push_back( *tIt );
  }

  nTot_  +=eles.size();
  nPassed_+=passingElectrons->size();

  delete [] isPassing;
  iEvent.put(passingElectrons);

}

// ------------ method called once each job just before starting event loop  ------------                                                                                                         
void HEEPElectronProducer::beginJob(){}

//______________________________________________________________________________                                                                                                                  
void HEEPElectronProducer::endJob(){

  if(debug_){
   std::cout<<" ############################ "<<std::endl;
   std::cout<<" ### HEEPElectronProducer ### "<<std::endl;
   std::cout<<" ############################ "<<std::endl;
   std::stringstream ss;
   ss<<"nTot="<<nTot_<<" nPassed="<<nPassed_<<" effPassed="<<100.*(nPassed_/(double)nTot_)<<"%\n";
   std::cout<<ss<<std::endl;
   std::cout<<" ############################ "<<std::endl;
  }

}

DEFINE_FWK_MODULE(HEEPElectronProducer);
