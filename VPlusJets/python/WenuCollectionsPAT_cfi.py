import FWCore.ParameterSet.Config as cms
from ElectroWeakAnalysis.VPlusJets.AllPassFilter_cfi import AllPassFilter
from SHarper.HEEPAnalyzer.HEEPSelectionCuts_cfi import *
 

def WenuCollectionsPAT(process,
                       patMuonCollection ,
                       patElectronCollection,
                       vertexCollection,
                       isQCD,
                       isHEEPID,
                       pTCutValue,
                       pTCutLooseMuonVeto,
                       pTCutLooseElectronVeto,
                       isTransverseMassCut,
                       TransverseMassCutValue,
                       patTypeICorrectedMetSysShifted):
    

 isMuonAnalyzer = False ## here we are looking for electron channel

 print "                                         "
 print "#########################################"
 print "## Electron Selection in W->enu events ##"
 print "#########################################"
 print "                                         "

 print "Chosen Options:                          "
 print "input pat muon collection                           = %s"%patMuonCollection
 print "input pat electron collection                       = %s"%patElectronCollection
 print "input vertex collection                             = %s"%vertexCollection
 print "is running QCD isolation cut                        = %d"%isQCD
 print "run the HEEP elctron ID instead of Higgs MVA ID     = %d"%isHEEPID
 print "chosen pT threshold for the ID                      = %f"%pTCutValue
 print "chosen pT threshold for loose muon veto             = %f"%pTCutLooseMuonVeto
 print "chosen pT threshold for loose electron veto         = %f"%pTCutLooseElectronVeto
 print "apply tranverse mass cut                            = %d"%isTransverseMassCut
 print "transverse mass cut value                           = %f"%TransverseMassCutValue
 print "patTypeICorrectedMetSysShifted                      = %s"%patTypeICorrectedMetSysShifted
 print "                                         "


 tightEleIdLabel = "tight" ## id label for tight electrons
 if isQCD:
  tightEleIdLabel = "qcd"

 ## modified WP70
 if isHEEPID:
   ### heep electrons producer starting from a set of pat electrons
   process.heepPatElectronsPFlow = cms.EDProducer("HEEPAttStatusToPAT",
                                                   eleLabel = patElectronCollection,
                                                   barrelCuts = cms.PSet(heepBarrelCuts),
                                                   endcapCuts = cms.PSet(heepEndcapCuts),
                                                   applyRhoCorrToEleIsol = cms.bool(True),
                                                   eleIsolEffectiveAreas = cms.PSet (heepEffectiveAreas),
                                                   eleRhoCorrLabel = cms.InputTag("kt6PFJetsPFlow","rho"),
                                                   verticesLabel = vertexCollection)  

   ### tight heep electrons definition through plugins/HEEPElectronProducer.cc
   process.tightElectrons  = cms.EDProducer("HEEPElectronProducer",
                                             electronCollection = cms.InputTag("heepPatElectronsPFlow"),
                                             eleIdType = cms.string("TightID"),
                                             pTCutValue = cms.double(pTCutValue))

   ## requrie at least one tight electrons
   process.tightElectronFilter = cms.EDFilter("PATCandViewCountFilter",
                                               minNumber = cms.uint32(1),
                                               maxNumber = cms.uint32(999999),
                                               src = cms.InputTag("tightElectrons"))
          
   process.tightLeptonStep = AllPassFilter.clone()

   ## build the leptonic W candidate with the met and the tight electron
   process.WToEnu = cms.EDProducer("CandViewShallowCloneCombiner",
                                    decay = cms.string("tightElectrons "+patTypeICorrectedMetSysShifted[0]),
                                    cut = cms.string('daughter(0).pt >20 && daughter(1).pt >20  && sqrt(2*daughter(0).pt*daughter(1).pt*(1-cos(daughter(0).phi-daughter(1).phi)))>0'),
                                    checkCharge = cms.bool(False))

   if isTransverseMassCut :
       process.WToEnu.cut = cms.string('daughter(1).pt >20  && sqrt(2*daughter(0).pt*daughter(1).pt*(1-cos(daughter(0).phi-daughter(1).phi)))>%f'%isTransverseMassCut)


   process.bestWToEnu =cms.EDFilter("LargestPtCandViewSelector",
                                     maxNumber = cms.uint32(10),
                                     src = cms.InputTag("WToEnu"))
      
   process.bestWToLepnuStep = AllPassFilter.clone()


   ### final sequence
   process.WSequence = cms.Sequence(process.heepPatElectronsPFlow*
                                    process.tightElectrons *
                                    process.tightElectronFilter *
                                    process.tightLeptonStep *
                                    process.WToEnu *
                                    process.bestWToEnu *
                                    process.bestWToLepnuStep)

   from ElectroWeakAnalysis.VPlusJets.LooseLeptonVetoPAT_cfi import LooseLeptonVetoPAT

   ## Lepton Veto Sequence
   LooseLeptonVetoPAT(process,
                      isQCD,
                      isHEEPID,
                      isMuonAnalyzer,
                      patMuonCollection,
                      pTCutLooseMuonVeto,
                      pTCutLooseElectronVeto)

   process.WPath = cms.Sequence(process.WSequence*
                                process.VetoSequence)
      
 else:

   ### Tight electron selection according to egamma mva Higgs id  
   process.tightElectrons = cms.EDProducer("PATElectronIdSelector",
                                            src = patElectronCollection,
                                            idLabel = cms.string(tightEleIdLabel),  # refers to Higgs Cut Based or MVA WP
                                            useMVAbasedID = cms.bool(True))

   ### require at least one tight electrons  
   process.tightElectronFilter = cms.EDFilter("PATCandViewCountFilter",
        minNumber = cms.uint32(1),
        maxNumber = cms.uint32(999999),
        src = cms.InputTag("tightElectrons"))

   ### tight ele filter
   process.tightLeptonStep = AllPassFilter.clone()


   process.WToEnu = cms.EDProducer("CandViewShallowCloneCombiner",
                                    decay = cms.string("tightElectrons "+patTypeICorrectedMetSysShifted[0]),
                                    cut = cms.string('daughter(1).pt >20  && sqrt(2*daughter(0).pt*daughter(1).pt*(1-cos(daughter(0).phi-daughter(1).phi)))>0'),
                                    checkCharge = cms.bool(False))

   if isTransverseMassCut :
       process.WToEnu.cut = cms.string('daughter(1).pt >20  && sqrt(2*daughter(0).pt*daughter(1).pt*(1-cos(daughter(0).phi-daughter(1).phi)))>%f'%TransverseMassCutValue)


   process.bestWToEnu =cms.EDFilter("LargestPtCandViewSelector",
                                   maxNumber = cms.uint32(10),
                                   src = cms.InputTag("WToEnu"))
   
   process.bestWToLepnuStep = AllPassFilter.clone()

   ### final sequence
   process.WSequence = cms.Sequence(process.tightElectrons *
                         process.tightElectronFilter *
                         process.tightLeptonStep *
                         process.WToEnu *
                         process.bestWToEnu *
                         process.bestWToLepnuStep
      )


   ### call python function for loose lepton producer and veto
   LooseLeptonVetoPAT(process,
                      isQCD,
                      isHEEPID,
                      isMuonAnalyzer,
                      patMuonCollection,
                      pTCutLooseMuonVeto,
                      pTCutLooseElectronVeto)
 
   process.WPath = cms.Sequence(process.WSequence*
                                process.VetoSequence)
