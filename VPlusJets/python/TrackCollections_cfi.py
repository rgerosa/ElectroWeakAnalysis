import FWCore.ParameterSet.Config as cms

from ElectroWeakAnalysis.VPlusJets.AllPassFilter_cfi import AllPassFilter

#--------------------------
# Counter1: All read events
AllEventsStep = AllPassFilter.clone()

##-------- Scraping veto --------
noscraping = cms.EDFilter("FilterOutScraping",
                           applyfilter = cms.untracked.bool(True),
                           debugOn = cms.untracked.bool(False),
                           numtrack = cms.untracked.uint32(10),
                           thresh = cms.untracked.double(0.2)
                          )

noscrapingStep = AllPassFilter.clone()

##---------HBHE Noise Filter ------
from CommonTools.RecoAlgos.HBHENoiseFilter_cfi import HBHENoiseFilter
HBHENoiseFilter.minIsolatedNoiseSumE = cms.double(999999.)
HBHENoiseFilter.minNumIsolatedNoiseChannels = cms.int32(999999)
HBHENoiseFilter.minIsolatedNoiseSumEt = cms.double(999999.)

HBHENoiseStep = AllPassFilter.clone()

## The CSC beam halo tight filter ____________________________________________||
from RecoMET.METAnalyzers.CSCHaloFilter_cfi import CSCBasedHaloFilter

CSCHaloStep = AllPassFilter.clone()

## The HCAL laser filter _____________________________________________________||
from RecoMET.METFilters.hcalLaserEventFilter_cfi import hcalLaserEventFilter
hcalLaserEventFilter.vetoByRunEventNumber=cms.untracked.bool(False)
hcalLaserEventFilter.vetoByHBHEOccupancy=cms.untracked.bool(True)

hcalLaserEventStep = AllPassFilter.clone()

## The ECAL dead cell trigger primitive filter _______________________________||
from RecoMET.METFilters.EcalDeadCellTriggerPrimitiveFilter_cfi import EcalDeadCellTriggerPrimitiveFilter
EcalDeadCellTriggerPrimitiveFilter.tpDigiCollection = cms.InputTag("ecalTPSkimNA")

EcalDeadCellTriggerPrimitiveStep = AllPassFilter.clone()

## The EE bad SuperCrystal filter ____________________________________________||
from RecoMET.METFilters.eeBadScFilter_cfi import eeBadScFilter

eeBadScStep = AllPassFilter.clone()

## The tracking failure filter _______________________________________________||
from RecoMET.METFilters.trackingFailureFilter_cfi import trackingFailureFilter
trackingFailureStep = AllPassFilter.clone()

from RecoMET.METFilters.trackingPOGFilters_cfi import *
trackingPOGStep     = AllPassFilter.clone()

# Tracking TOBTEC fakes filter ##
tobtecfakesfilter.filter=cms.bool(False)
tobtecfakesStep   = AllPassFilter.clone()

primaryVertexFilter = cms.EDFilter("VertexSelector",
                                    src = cms.InputTag("offlinePrimaryVertices"),
                                    cut = cms.string("!isFake & ndof > 4 & abs(z) <= 24 & position.Rho <= 2"),
                                    filter = cms.bool(True))

primaryVertexStep = AllPassFilter.clone()

TrackVtxPath = cms.Sequence(AllEventsStep*
                            noscraping*
                            noscrapingStep*
                            HBHENoiseFilter*
                            HBHENoiseStep*
                            CSCBasedHaloFilter*
                            CSCHaloStep*
                            hcalLaserEventFilter*
                            hcalLaserEventStep*
                            EcalDeadCellTriggerPrimitiveFilter*
                            hcalLaserEventStep*
                            EcalDeadCellTriggerPrimitiveFilter*
                            EcalDeadCellTriggerPrimitiveStep*
                            eeBadScFilter*
                            eeBadScStep*
                            trackingFailureFilter*
                            trackingFailureStep*
                            ~manystripclus53X *
                            ~toomanystripclus53X *
                            ~logErrorTooManyClusters *
                            ~logErrorTooManyTripletsPairs *
                            ~logErrorTooManySeeds *
                            trackingPOGStep*
                            tobtecfakesfilter*
                            tobtecfakesStep*
                            primaryVertexFilter*
                            primaryVertexFilter)
