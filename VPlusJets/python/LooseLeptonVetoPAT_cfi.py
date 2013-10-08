import FWCore.ParameterSet.Config as cms
from ElectroWeakAnalysis.VPlusJets.AllPassFilter_cfi import AllPassFilter
from SHarper.HEEPAnalyzer.HEEPSelectionCuts_cfi import *

def LooseLeptonVetoPAT(process,
                       isQCD,
                       isHEEPID,
                       isMuonAnalyzer,
                       patMuonCollection,
                       pTCutLooseMuonVeto,
                       patElectronCollection = cms.InputTag("selectedPatElectronsPFlow"),
                       vertexCollection = cms.InputTag("goodOfflinePrimaryVertices"),
                       looseEleIdLabel="loose"):


 print "                               "
 print "###############################"
 print "## Loose Lepton veto Moduel  ##"
 print "###############################"
 print "                               "

 print "Chosen Options:                      "
 print "input pat muon collection                                         = %s"%patMuonCollection
 print "input pat electron collection                                     = %s"%patElectronCollection
 print "input vertex collection                                           = %s"%vertexCollection
 print "is running QCD isolation cut                                      = %d"%isQCD
 print "run the High pT muon ID and HEEP electron instead of Higgs ones   = %d"%isHEEPID
 print "flag to understand we are looking for W->munu or W->enu           = %d"%isMuonAnalyzer
 print "chosen pT threshold for loose muon veto                           = %f"%pTCutLooseMuonVeto
 print "id label considered for loose electron veto                       = %s"%looseEleIdLabel
 print "                                      "

                
 if isHEEPID :

   ## produce loose muon collection
   process.looseMuons = cms.EDFilter("PATMuonRefSelector",
                                                  src = patMuonCollection,
                                                  cut = cms.string(""))

   ## set the cut to the loose Hight pT muon ID
   process.looseMuons.cut = cms.string(" isGlobalMuon && isTrackerMuon && pt>%f && abs(eta) < 2.4 && abs(phi)<3.2 && trackIso/pt < 0.1 && abs(dB) <0.2"
                                       " && globalTrack().hitPattern().numberOfValidPixelHits>0 "
                                       " && globalTrack().hitPattern().numberOfValidMuonHits>0 "
                                       " && globalTrack().hitPattern().trackerLayersWithMeasurement>5 "
                                       " && numberOfMatchedStations>1 && ptError/pt<0.3"%pTCutLooseMuonVeto)

   ## filter the event where there are loose muons --> 
   process.looseMuonFilter = cms.EDFilter("PATCandViewCountFilter",
                                           src = cms.InputTag("looseMuons"))

   process.looseMuonStep = AllPassFilter.clone()
   
   if isMuonAnalyzer:
    ## if W->munu analysis is running, we need to produce heep electron from patElectronCollection --> this is why the HEEPID flag is set to true
    process.heepPatElectronsPFlow = cms.EDProducer("HEEPAttStatusToPAT",
                                                    eleLabel = patElectronCollection,
                                                    barrelCuts = cms.PSet(heepBarrelCuts),
                                                    endcapCuts = cms.PSet(heepEndcapCuts),
                                                    applyRhoCorrToEleIsol = cms.bool(True),
                                                    eleIsolEffectiveAreas = cms.PSet (heepEffectiveAreas),
                                                    eleRhoCorrLabel = cms.InputTag("kt6PFJetsPFlow","rho"),
                                                    verticesLabel = vertexCollection)

   ## define the loose heep electron collection
   process.looseElectrons = cms.EDProducer("HEEPElectronProducer",
                                            electronCollection = patElectronCollection,
                                            eleIdType = cms.string("LooseID"),
                                            pTCutValue = cms.double(20.))

   if isMuonAnalyzer :
     process.looseElectrons.electronCollection = "heepPatElectronsPFlow"

   ## select events with a least one for the moment
   process.looseElectronFilter = cms.EDFilter("PATCandViewCountFilter",
                                               src = cms.InputTag("looseElectrons"),
                                               minNumber = cms.uint32(1),
                                               maxNumber = cms.uint32(999))


   if isMuonAnalyzer : ## if the W->munu analysis is running, we need exactly one loose muon == tight and no loose electrons
                       process.looseMuonFilter.minNumber = cms.uint32(1)
                       process.looseMuonFilter.maxNumber = cms.uint32(1)
                       process.looseElectronFilter.minNumber = cms.uint32(0)
                       process.looseElectronFilter.maxNumber = cms.uint32(0)
   else :              ## if the W->enu analysis is running, we need exactly one loose electron == tight and no loose muons
                       process.looseMuonFilter.minNumber = cms.uint32(0)
                       process.looseMuonFilter.maxNumber = cms.uint32(0)
                       process.looseElectronFilter.minNumber = cms.uint32(1)
                       process.looseElectronFilter.maxNumber = cms.uint32(1)
                      
   process.looseElectronStep = AllPassFilter.clone()

   if isMuonAnalyzer:
     process.VetoSequence = cms.Sequence(process.looseMuons*
                                         process.looseMuonFilter*
                                         process.looseMuonStep*
                                         process.heepPatElectronsPFlow*
                                         process.looseElectrons*
                                         process.looseElectronFilter*
                                         process.looseElectronStep)
                  
   else: process.VetoSequence = cms.Sequence(process.looseElectrons*
                                             process.looseElectronFilter*
                                             process.looseElectronStep*
                                             process.looseMuons*
                                             process.looseMuonFilter*
                                             process.looseMuonStep)
                                                
 else:

  ##  Define loose electron selection for veto ######
  process.looseElectrons = cms.EDProducer("PATElectronIdSelector",
                                           src = patElectronCollection,
                                           idLabel = cms.string(looseEleIdLabel),
                                           useMVAbasedID_   = cms.bool(True))

  ##  Define loose muon selection for veto ######
  process.looseMuons = cms.EDFilter("PATMuonRefSelector",
                                     src = patMuonCollection,
                                     cut = cms.string(""))

  process.looseMuons.cut = cms.string(" pt>%f &&isPFMuon && (isGlobalMuon || isTrackerMuon) && abs(eta)<2.4 && (pfIsolationR04().sumChargedHadronPt+"
                                      " max(0.,pfIsolationR04().sumNeutralHadronEt+pfIsolationR04().sumPhotonEt-0.5*pfIsolationR04().sumPUPt))/pt< 0.2"%pTCutLooseMuonVeto)

  if isMuonAnalyzer : ## for W->munu events just one tight muon and no loose electrons
                      nLooseElectron = 0 ;
                      nLooseMuon     = 1 ;
  else :              ## for W->enu events just one tight electron and no loose muons
                      nLooseElectron = 1 ;
                      nLooseMuon     = 0 ;
 
  process.looseElectronFilter = cms.EDFilter("PATCandViewCountFilter",
                                              minNumber = cms.uint32(nLooseElectron),
                                              maxNumber = cms.uint32(nLooseElectron),
                                              src = cms.InputTag("looseElectrons"))

  process.looseElectronStep = AllPassFilter.clone()


  process.looseMuonFilter = cms.EDFilter("PATCandViewCountFilter",
                                          minNumber = cms.uint32(nLooseMuon),
                                          maxNumber = cms.uint32(nLooseMuon),
                                          src = cms.InputTag("looseMuons"))

  process.looseMuonStep = AllPassFilter.clone()
   
  process.VetoSequence = cms.Sequence(process.looseElectrons*
                                      process.looseElectronFilter*
                                      process.looseElectronStep*
                                      process.looseMuons*
                                      process.looseMuonFilter*
                                      process.looseMuonStep)
