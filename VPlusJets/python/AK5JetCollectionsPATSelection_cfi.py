import FWCore.ParameterSet.Config as cms

from ElectroWeakAnalysis.VPlusJets.AllPassFilter_cfi import AllPassFilter
from PhysicsTools.SelectorUtils.pfJetIDSelector_cfi  import pfJetIDSelector

from RecoJets.Configuration.GenJetParticles_cff import *
from RecoJets.JetProducers.ak5GenJets_cfi       import *
from RecoMET.Configuration.GenMETParticles_cff  import *
from RecoMET.METProducers.genMetTrue_cfi        import *

##########################################################################


def AK5JetCollectionsPATSelection(process,
                                  isMC,
                                  patJetCollection,
                                  patSmearedJetCollection,
                                  isPileUpJetID,
                                  useMVAPileUpJetID,
                                  useSmearedCollection,
                                  storeSmearandShiftCollections,
                                  jetPtThreshold,
                                  isRequireTwoJets):


 print "                                  "
 print "##################################"
 print "## AK5 Jet collection Selection ##"
 print "##################################"
 print "                                  "

 print "Chosen Options:                           "
 print "is running on data or on MC                   = %d"%isMC
 print "input pat AK5 jet collection                  = %s"%patJetCollection
 print "input pat Smeared AK5 jet collection          = %s"%patSmearedJetCollection
 print "run Pile Up Jet ID                            = %d"%isPileUpJetID
 print "use the mva chs pile up jet id                = %d"%useMVAPileUpJetID 
 print "use Smeared jet collection for selections     = %d"%useSmearedCollection
 print "run analysis also on shifted and smeared jets = %d"%storeSmearandShiftCollections
 print "jet pT threshold to be applied                = %f"%jetPtThreshold
 print "jet two separated jets --> resolved analysis  = %d"%isRequireTwoJets
 print "                                  "



 ### produce PileUp PF jet ID collection using cur based puJetIdChs
 process.ak5PFnoPUJets = cms.EDProducer("PATPuJetIdSelector",
                                         src = cms.InputTag(patJetCollection[0]),
                                         idLabel = cms.string("loose"),
                                         valueMapLabel = cms.string("puJetMvaChs"),
                                         applyMVAID = cms.bool(False))

 ### in case of mva ID , change the value map
 if useMVAPileUpJetID:
    process.ak5PFnoPUJets.applyMVAID = cms.bool(True)


 # Apply loose PF jet ID
 process.ak5PFGoodJets = cms.EDFilter("PFJetIDSelectionFunctorFilter",
                                       filterParams = pfJetIDSelector.clone(),
                                       src = cms.InputTag("ak5PFnoPUJets"),
                                       filter = cms.bool(True))

 if not isPileUpJetID:
     process.ak5PFGoodJets.src = patJetCollection[0]


 ### clean jets from electron and muons   --> by default this use the loose jet id collection, looseMuons and looseElectrons   
 process.ak5PFJetsClean = cms.EDProducer("PFPATJetCleaner",
                                          srcJets = cms.InputTag("ak5PFGoodJets"),
                                          module_label = cms.string(""),
                                          srcObjects = cms.VInputTag(cms.InputTag("looseElectrons"),cms.InputTag("looseMuons")),
                                          deltaRMin = cms.double(0.3))


 ### pT selection on top of the cleaned jet collection --> all eta, only central and only forward for both smeared and not smeared collection
     
 process.ak5PFJetsPtSkimmed = cms.EDFilter("PATJetRefSelector",
                                            src = cms.InputTag("ak5PFJetsClean"),
                                            cut = cms.string('pt > %f'%jetPtThreshold))


 process.ak5PFJetsPtSkimmedCentral = cms.EDFilter("PATJetRefSelector",
                                                   src = cms.InputTag("ak5PFJetsClean"),
                                                   cut = cms.string('pt > %f && abs(eta) < 2.4'%jetPtThreshold))

 process.ak5PFJetsPtSkimmedForward = cms.EDFilter("PATJetRefSelector",
                                                   src = cms.InputTag("ak5PFJetsClean"),
                                                   cut = cms.string('pt > %f && abs(eta) > 2.4'%jetPtThreshold))



 ### Filter to require at least two jets in the event , to be applied for two jets topology --> after pT skim

 process.RequireTwoJets = cms.EDFilter("PATCandViewCountFilter",
                                        minNumber = cms.uint32(2),
                                        maxNumber = cms.uint32(100),
                                        src = cms.InputTag("ak5PFJetsPtSkimmed"))

 process.RequireTwoJetsStep = AllPassFilter.clone()
 
 ### Define the most general sequence
 
 process.ak5PFJetPath = cms.Sequence( process.ak5PFnoPUJets*
                                      process.ak5PFGoodJets*
                                      process.ak5PFJetsClean*
                                      process.ak5PFJetsPtSkimmed*
                                      process.ak5PFJetsPtSkimmedCentral*
                                      process.ak5PFJetsPtSkimmedForward)

 if isRequireTwoJets :
   process.ak5PFJetPath += process.RequireTwoJets*process.RequireTwoJetsStep
                                            
 if not isPileUpJetID :
   process.ak5PFJetPath.remove(process.ak5PFnoPUJets)


 ### run the same analysis also on the other jet collection only for the MC   
 
 if isMC :

  for module in patJetCollection :

   noPUJetName              = 'ak5PFnoPUJets'
   mvaChsValueMap           = 'puJetMvaChs'
   goodJetName              = 'ak5PFGoodJets'
   cleanJetName             = 'ak5PFJetsClean'
   JetsPtSkimmedName        = 'ak5PFJetsPtSkimmed' 
   JetsPtSkimmedNameCentral = 'ak5PFJetsPtSkimmedCentral'
   JetsPtSkimmedNameForward = 'ak5PFJetsPtSkimmedForward'

   if module == patJetCollection[0] : continue

   if "EnUp" in module and "shifted" in module:

       noPUJetName              = noPUJetName+'EnUp'
       mvaChsValueMap           = mvaChsValueMap+'EnUp'
       goodJetName              = goodJetName+'EnUp'
       cleanJetName             = cleanJetName+'EnUp'
       JetsPtSkimmedName        = JetsPtSkimmedName+'EnUp'
       JetsPtSkimmedNameCentral = JetsPtSkimmedNameCentral+'EnUp'
       JetsPtSkimmedNameForward = JetsPtSkimmedNameForward+'EnUp'
       
   elif "EnDown" in module and "shifted" in module :

       noPUJetName              = noPUJetName+'EnDown'
       mvaChsValueMap           = mvaChsValueMap+'EnDown'
       goodJetName              = goodJetName+'EnDown'
       cleanJetName             = cleanJetName+'EnDown'
       JetsPtSkimmedName        = JetsPtSkimmedName+'EnDown'
       JetsPtSkimmedNameCentral = JetsPtSkimmedNameCentral+'EnDown'
       JetsPtSkimmedNameForward = JetsPtSkimmedNameForward+'EnDown'

   setattr(process, noPUJetName, process.ak5PFnoPUJets.clone(src = cms.InputTag(module),
                                                             valueMapLabel = cms.string(mvaChsValueMap)))
   setattr(process, goodJetName, process.ak5PFGoodJets.clone(src = cms.InputTag(noPUJetName)))

   if not isPileUpJetID :
    setattr(process, goodJetName, process.ak5PFGoodJets.clone(src = cms.InputTag(module)))
           
   setattr(process, cleanJetName, process.ak5PFJetsClean.clone(src = cms.InputTag(goodJetName)))
   setattr(process, JetsPtSkimmedName, process.ak5PFJetsPtSkimmed.clone(src = cms.InputTag(cleanJetName)))
   setattr(process, JetsPtSkimmedNameCentral, process.ak5PFJetsPtSkimmedCentral.clone(src = cms.InputTag(cleanJetName)))
   setattr(process, JetsPtSkimmedNameForward, process.ak5PFJetsPtSkimmedForward.clone(src = cms.InputTag(cleanJetName)))

   process.ak5PFJetPath += getattr(process,noPUJetName)
   process.ak5PFJetPath += getattr(process,goodJetName)
   process.ak5PFJetPath += getattr(process,cleanJetName)
   process.ak5PFJetPath += getattr(process,JetsPtSkimmedName)
   process.ak5PFJetPath += getattr(process,JetsPtSkimmedNameCentral)
   process.ak5PFJetPath += getattr(process,JetsPtSkimmedNameForward)           

   if not storeSmearandShiftCollections:

       process.ak5PFJetPath.remove(getattr(process,noPUJetName))
       process.ak5PFJetPath.remove(getattr(process,cleanJetName))
       process.ak5PFJetPath.remove(getattr(process,goodJetName))
       process.ak5PFJetPath.remove(getattr(process,JetsPtSkimmedName))
       process.ak5PFJetPath.remove(getattr(process,JetsPtSkimmedNameCentral))
       process.ak5PFJetPath.remove(getattr(process,JetsPtSkimmedNameForward))

   if not isPileUpJetID :
       process.ak5PFJetPath.remove(getattr(process,noPUJetName))

 
  if useSmearedCollection :
      
   for module in patSmearedJetCollection :

    noPUJetName              = 'smearedak5PFnoPUJets'
    mvaChsValueMap           = 'SmearedJetMvaChs'

    if module == patSmearedJetCollection[0] :
        mvaChsValueMap       = 'puSmearedJetMvaChs'

    goodJetName              = 'smearedak5PFGoodJets'
    cleanJetName             = 'smearedak5PFJetsClean'
    JetsPtSkimmedName        = 'smearedak5PFJetsPtSkimmed' 
    JetsPtSkimmedNameCentral = 'smearedak5PFJetsPtSkimmedCentral'
    JetsPtSkimmedNameForward = 'smearedak5PFJetsPtSkimmedForward'

    if "EnUp" in module and ( "smeared" in module or "Smeared" in module):

       noPUJetName              = noPUJetName+'EnUp'
       mvaChsValueMap           = 'pushifted'+mvaChsValueMap+'EnUp'
       goodJetName              = goodJetName+'EnUp'
       cleanJetName             = cleanJetName+'EnUp'
       JetsPtSkimmedName        = JetsPtSkimmedName+'EnUp'
       JetsPtSkimmedNameCentral = JetsPtSkimmedNameCentral+'EnUp'
       JetsPtSkimmedNameForward = JetsPtSkimmedNameForward+'EnUp'
       
    elif "EnDown" in module and  ( "smeared" in module or "Smeared" in module):

       noPUJetName              = noPUJetName+'EnDown'
       mvaChsValueMap           = 'pushifted'+mvaChsValueMap+'EnDown'
       goodJetName              = goodJetName+'EnDown'
       cleanJetName             = cleanJetName+'EnDown'
       JetsPtSkimmedName        = JetsPtSkimmedName+'EnDown'
       JetsPtSkimmedNameCentral = JetsPtSkimmedNameCentral+'EnDown'
       JetsPtSkimmedNameForward = JetsPtSkimmedNameForward+'EnDown'

    elif "ResDown" in module and ( "smeared" in module or "Smeared" in module):

       noPUJetName              = noPUJetName+'ResDown'
       goodJetName              = goodJetName+'ResDown'
       mvaChsValueMap           = 'pu'+mvaChsValueMap+'ResDown'
       cleanJetName             = cleanJetName+'ResDown'
       JetsPtSkimmedName        = JetsPtSkimmedName+'ResDown'
       JetsPtSkimmedNameCentral = JetsPtSkimmedNameCentral+'ResDown'
       JetsPtSkimmedNameForward = JetsPtSkimmedNameForward+'ResDown'

    elif "ResUp" in module and ( "smeared" in module or "Smeared" in module):

       noPUJetName              = noPUJetName+'ResUp'
       mvaChsValueMap           = 'pu'+mvaChsValueMap+'ResUp'
       goodJetName              = goodJetName+'ResUp'
       cleanJetName             = cleanJetName+'ResUp'
       JetsPtSkimmedName        = JetsPtSkimmedName+'ResUp'
       JetsPtSkimmedNameCentral = JetsPtSkimmedNameCentral+'ResUp'
       JetsPtSkimmedNameForward = JetsPtSkimmedNameForward+'ResUp'

    setattr(process, noPUJetName, process.ak5PFnoPUJets.clone(src = cms.InputTag(module),
                                                             valueMapLabel = cms.string(mvaChsValueMap)))
    setattr(process, goodJetName, process.ak5PFGoodJets.clone(src = cms.InputTag(noPUJetName)))

    if not isPileUpJetID :
     setattr(process, goodJetName, process.ak5PFGoodJets.clone(src = cms.InputTag(module)))
           
    setattr(process, cleanJetName, process.ak5PFJetsClean.clone(src = cms.InputTag(goodJetName)))
    setattr(process, JetsPtSkimmedName, process.ak5PFJetsPtSkimmed.clone(src = cms.InputTag(cleanJetName)))
    setattr(process, JetsPtSkimmedNameCentral, process.ak5PFJetsPtSkimmedCentral.clone(src = cms.InputTag(cleanJetName)))
    setattr(process, JetsPtSkimmedNameForward, process.ak5PFJetsPtSkimmedForward.clone(src = cms.InputTag(cleanJetName)))

    process.ak5PFJetPath += getattr(process,noPUJetName)
    process.ak5PFJetPath += getattr(process,goodJetName)
    process.ak5PFJetPath += getattr(process,cleanJetName)
    process.ak5PFJetPath += getattr(process,JetsPtSkimmedName)
    process.ak5PFJetPath += getattr(process,JetsPtSkimmedNameCentral)
    process.ak5PFJetPath += getattr(process,JetsPtSkimmedNameForward)           

    if not storeSmearandShiftCollections and module != patSmearedJetCollection[0] :

     process.ak5PFJetPath.remove(getattr(process,noPUJetName))
     process.ak5PFJetPath.remove(getattr(process,goodJetName))
     process.ak5PFJetPath.remove(getattr(process,cleanJetName))
     process.ak5PFJetPath.remove(getattr(process,JetsPtSkimmedName))
     process.ak5PFJetPath.remove(getattr(process,JetsPtSkimmedNameCentral))
     process.ak5PFJetPath.remove(getattr(process,JetsPtSkimmedNameForward))           
        
        
    if not isPileUpJetID :
       process.ak5PFJetPath.remove(getattr(process,noPUJetName))
   

 ################################################
 ### Gen AK5 Jets and partons --> only for MC ###
 ################################################


 process.genPartons = cms.EDProducer("PartonSelector",
                                      src = cms.InputTag("genParticles"),
                                      withLeptons = cms.bool(False))

 process.ak5flavourByRef = cms.EDProducer("JetPartonMatcher",
                                           jets = cms.InputTag(patJetCollection[0]),
                                           coneSizeToAssociate = cms.double(0.3),
                                           partons = cms.InputTag("genPartons"))


 process.ak5SmearedflavourByRef = process.ak5flavourByRef.clone(jets = cms.InputTag(patSmearedJetCollection[0]))

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

