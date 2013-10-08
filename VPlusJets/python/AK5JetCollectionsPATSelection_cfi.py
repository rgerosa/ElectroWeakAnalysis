import FWCore.ParameterSet.Config as cms

from ElectroWeakAnalysis.VPlusJets.AllPassFilter_cfi import AllPassFilter
from PhysicsTools.SelectorUtils.pfJetIDSelector_cfi  import pfJetIDSelector

from RecoJets.Configuration.GenJetParticles_cff import *
from RecoJets.JetProducers.ak5GenJets_cfi       import *
from RecoMET.Configuration.GenMETParticles_cff  import *
from RecoMET.METProducers.genMetTrue_cfi        import *

##########################################################################

def AK5JetCollectionsPATSelection(process,                                  
                                  patJetCollection,
                                  patSmearedJetCollection,
                                  isPileUpJetID,
                                  useMVAPileUpJetID,
                                  useSmearedCollection,
                                  jetPtThreshold,
                                  isRequireTwoJets,
                                  isMC):


 print "                                  "
 print "##################################"
 print "## AK5 Jet collection Selection ##"
 print "##################################"
 print "                                  "

 print "Chosen Options:                           "
 print "input pat AK5 jet collection                  = %s"%patJetCollection
 print "input pat Smeared AK5 jet collection          = %s"%patSmearedJetCollection
 print "run Pile Up Jet ID                            = %d"%isPileUpJetID
 print "use the mva chs pile up jet id                = %d"%useMVAPileUpJetID 
 print "use Smeared jet collection for selections     = %d"%useSmearedCollection
 print "jet pT threshold to be applied                = %f"%jetPtThreshold
 print "jet two separated jets --> resolved analysis  = %d"%isRequireTwoJets
 print "is running on data or on MC                   = %d"%isMC
 print "                                  "

                  ### produce PileUp PF jet ID collection using cur based puJetIdChs
 process.ak5PFnoPUJets = cms.EDProducer("PATPuJetIdSelector",
                                         src = patJetCollection,
                                         idLabel = cms.string("loose"),
                                         valueMapLabel = cms.string("puJetMvaChs"),
                                         applyMVAID = cms.bool(False))

 ### in case of mva ID , change the value map
 if useMVAPileUpJetID:
    process.ak5PFnoPUJets.applyMVAID = cms.bool(True)


 ### produce loose pile up jet id collection on smeared jets, using cut based id
 process.ak5PFnoPUSmearedJets = process.ak5PFnoPUJets.clone( src = patSmearedJetCollection,
                                                             valueMapLabel = cms.string("puSmearedJetMvaChs"),
                                                             applyMVAID = cms.bool(False))


 ### in case of mva ID , change the value map
 if useMVAPileUpJetID :
    process.ak5PFnoPUSmearedJets.applyMVAID = cms.bool(True) 

 # Apply loose PF jet ID
 process.ak5PFGoodJets = cms.EDFilter("PFJetIDSelectionFunctorFilter",
                                       filterParams = pfJetIDSelector.clone(),
                                       src = cms.InputTag("ak5PFnoPUJets"),
                                       filter = cms.bool(True))

 if not isPileUpJetID:
     process.ak5PFGoodJets.src = patJetCollection

 ### apply loose jet iD on  Smeared jet collection 
 process.ak5PFGoodSmearedJets = process.ak5PFGoodJets.clone(src = cms.InputTag("ak5PFnoPUSmearedJets"))

 if not isPileUpJetID:
    process.ak5PFGoodJets.src = patSmearedJetCollection


       
 ### clean jets from electron and muons   --> by default this use the loose jet id collection, looseMuons and looseElectrons   
 process.ak5PFJetsClean = cms.EDProducer("PFPATJetCleaner",
                                          srcJets = cms.InputTag("ak5PFGoodJets"),
                                          module_label = cms.string(""),
                                          srcObjects = cms.VInputTag(cms.InputTag("looseElectrons"),cms.InputTag("looseMuons")),
                                          deltaRMin = cms.double(0.3))


 ### Cleaning of smeared jet collection 
 process.ak5PFSmearedJetsClean = process.ak5PFJetsClean.clone(srcJets = cms.InputTag("ak5PFGoodSmearedJets"))


 ### pT selection on top of the cleaned jet collection --> all eta, only central and only forward for both smeared and not smeared collection
     
 process.ak5PFJetsPtSkimmed = cms.EDFilter("PATJetRefSelector",
                                            src = cms.InputTag("ak5PFJetsClean"),
                                            cut = cms.string('pt > %f'%jetPtThreshold))


 if useSmearedCollection and isMC:
     process.ak5PFJetsPtSkimmed.src = cms.InputTag("ak5PFSmearedJetsClean")

 process.ak5PFJetsPtSkimmedCentral = cms.EDFilter("PATJetRefSelector",
                                                   src = cms.InputTag("ak5PFJetsClean"),
                                                   cut = cms.string('pt > %f && abs(eta) < 2.4'%jetPtThreshold))


 if useSmearedCollection and isMC:
     process.ak5PFJetsPtSkimmedCentral.src = cms.InputTag("ak5PFSmearedJetsClean")
 
 process.ak5PFJetsPtSkimmedForward = cms.EDFilter("PATJetRefSelector",
                                                 src = cms.InputTag("ak5PFJetsClean"),
                                                 cut = cms.string('pt > %f && abs(eta) > 2.4 && abs(eta) < 9.9'%jetPtThreshold))

 if useSmearedCollection and isMC:
     process.ak5PFJetsPtSkimmedForward.src = cms.InputTag("ak5PFSmearedJetsClean")


 ### Filter to require at least two jets in the event , to be applied for two jets topology --> after pT skim

 process.RequireTwoJets = cms.EDFilter("PATCandViewCountFilter",
                                        minNumber = cms.uint32(2),
                                        maxNumber = cms.uint32(100),
                                        src = cms.InputTag("ak5PFJetsPtSkimmed"))

 ### Define the most general sequence
 
 process.ak5PFJetPath = cms.Sequence( process.ak5PFnoPUJets*
                                      process.ak5PFnoPUSmearedJets*
                                      process.ak5PFGoodJets*
                                      process.ak5PFGoodSmearedJets*
                                      process.ak5PFJetsClean*
                                      process.ak5PFSmearedJetsClean*
                                      process.ak5PFJetsPtSkimmed*
                                      process.ak5PFJetsPtSkimmedCentral*
                                      process.ak5PFJetsPtSkimmedForward)


 if isRequireTwoJets :
   process.ak5PFJetPath += process.RequireTwoJets
                                            
 if not isPileUpJetID :
   process.ak5PFJetPath.remove(process.ak5PFnoPUJets)
   process.ak5PFJetPath.remove(process.ak5PFnoPUSmearedJets)
   

 if not useSmearedCollection :
   process.ak5PFJetPath.remove(process.ak5PFGoodSmearedJets)
   process.ak5PFJetPath.remove(process.ak5PFnoPUSmearedJets)
   process.ak5PFJetPath.remove(process.ak5PFSmearedJetsClean)

 if not isMC :
   process.ak5PFJetPath.remove(process.ak5PFGoodSmearedJets)
   process.ak5PFJetPath.remove(process.ak5PFnoPUSmearedJets)
   process.ak5PFJetPath.remove(process.ak5PFSmearedJetsClean)
            
     
   
 ################################################
 ### Gen AK5 Jets and partons --> only for MC ###
 ################################################


 process.genPartons = cms.EDProducer("PartonSelector",
                                      src = cms.InputTag("genParticles"),
                                      withLeptons = cms.bool(False))

 process.ak5flavourByRef = cms.EDProducer("JetPartonMatcher",
                                           jets = patJetCollection,
                                           coneSizeToAssociate = cms.double(0.3),
                                           partons = cms.InputTag("genPartons"))


 process.ak5SmearedflavourByRef = process.ak5flavourByRef.clone(jets = patSmearedJetCollection)

 process.ak5tagJet = cms.EDProducer("JetFlavourIdentifier",
                                     srcByReference = cms.InputTag("ak5flavourByRef"),
                                     physicsDefinition = cms.bool(False))

 process.ak5tagSmearedJet = process.ak5tagJet.clone(srcByReference = cms.InputTag("ak5SmearedflavourByRef"))
                                                    
 process.genTagJetPath = cms.Sequence( process.genPartons*
                                       process.ak5flavourByRef*
                                       process.ak5SmearedflavourByRef*
                                       process.ak5tagJet*
                                       process.ak5tagSmearedJet)

 if not useSmearedCollection :

  process.genTagJetPath.remove(process.ak5SmearedflavourByRef)
  process.genTagJetPath.remove(process.ak5tagSmearedJet)
 

 process.GenJetPath = cms.Sequence( process.genParticlesForJets*
                                    process.ak5GenJets*
                                    process.genParticlesForMETAllVisible*
                                    process.genMetTrue)

