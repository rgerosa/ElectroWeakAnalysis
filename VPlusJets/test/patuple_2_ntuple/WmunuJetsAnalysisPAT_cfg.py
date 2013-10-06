import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
import pprint

##########################
### Parsing Parameters ###
##########################

options = VarParsing ('python')

options.register ('isMC', True, VarParsing.multiplicity.singleton, VarParsing.varType.int,
                   "Run this on simulation")

options.register ('isQCD', False, VarParsing.multiplicity.singleton, VarParsing.varType.int,
                   "Run the QCD Isolation Selection")

options.register ('isTransverseMassCut', False, VarParsing.multiplicity.singleton, VarParsing.varType.int,
                   "Apply a mT cut on the leptonic leg (W) of the selected events")

options.register ('isHEEPID', True, VarParsing.multiplicity.singleton, VarParsing.varType.int,
                   "Use the HEEP electron ID instead of the MVA one")

options.register ('globalTag', '', VarParsing.multiplicity.singleton, VarParsing.varType.string,
                   "Global Tag to be used")

options.register ('numEventsToRun', -1, VarParsing.multiplicity.singleton, VarParsing.varType.int,
                   "Number of events to process: -1 means all the input events")

options.register ('OutputFileName', 'WmunuJetAnalysisntuple.root', VarParsing.multiplicity.singleton, VarParsing.varType.string,
                   "Name of the output file")

options.register ('reportEvery', 500, VarParsing.multiplicity.singleton, VarParsing.varType.int,
                   "Number of events after which print the report")

options.register ('runMetFilters', False, VarParsing.multiplicity.singleton, VarParsing.varType.int,
                   "Run the python/TrackingCollections_cfi to apply met filters and primary vertex collection skim")

options.register ('useSmearedCollection', True, VarParsing.multiplicity.singleton, VarParsing.varType.int,
                   "In case of MC analysis use smearedJets and smearedMet as default")

options.register ('hltPath', 'HLT_IsoMu24_*, HLT_Mu40_eta2p1*', VarParsing.multiplicity.singleton, VarParsing.varType.string,
                   "List of HLT path to be required when running on data. Should be separetad by ")

options.register ('didPhotonSmearing', True, VarParsing.multiplicity.singleton, VarParsing.varType.int,
                   "true if in PAT the photon smearing for resolution has been done")

options.register ('didTauSmearing', True, VarParsing.multiplicity.singleton, VarParsing.varType.int,
                   "true if in PAT the tau smearing for resolution has been done")


options.parseArguments()

############################
#### Create the Process ####
############################

process = cms.Process("demo")

##---------  Load standard Reco modules ------------
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')

##----- this config frament brings you the generator information ----
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("PhysicsTools.HepMCCandAlgos.genParticles_cfi")
process.load("Configuration.StandardSequences.Generator_cff")

##----- Detector geometry : some of these needed for b-tag -------
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi")
process.load("Geometry.CommonDetUnit.globalTrackingGeometry_cfi")
process.load("RecoMuon.DetLayers.muonDetLayerGeometry_cfi")

##----- B-tags --------------
process.load("RecoBTag.Configuration.RecoBTag_cff")

##----- Global tag: conditions database ------------
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

##----- Counter module ------------
process.load("ElectroWeakAnalysis.VPlusJets.AllPassFilter_cfi")


#############################
### Set of the Global Tag ###
#############################

if not options.isMC :
    if options.globalTag !='' :
     process.GlobalTag.globaltag = cms.string(globalTag)
    else :
     process.GlobalTag.globaltag = 'GR_P_V39_AN3::All'        
else:
    if options.globalTag !='' :
     process.GlobalTag.globaltag = cms.string(globalTag)
    else:
     process.GlobalTag.globaltag = 'START53_V7E::All'

###############################
### Input Module Definition ###
###############################

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(''))

if options.isMC:
    process.source.fileNames = cms.untracked.vstring('file:/data2/rgerosa/tlbsm_53x_v3_mc.root')
else :
    process.source.fileNames = cms.untracked.vstring('file:/data2/rgerosa/tlbsm_53x_v3_data.root')
    

################################
### Output Module Definition ###
################################     

process.maxEvents                                = cms.untracked.PSet(input = cms.untracked.int32(options.numEventsToRun))
process.MessageLogger.destinations               = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = options.reportEvery
process.options                                  = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )
#process.options                                 = cms.untracked.PSet( SkipEvent = cms.untracked.vstring('ProductNotFound'))

process.TFileService = cms.Service("TFileService", 
                                    fileName      = cms.string( options.OutputFileName ),
                                    closeFileFast = cms.untracked.bool(False))

##---------  Vertex and track Collections -----------
process.load("ElectroWeakAnalysis.VPlusJets.TrackCollections_cfi")

process.primaryVertexFilter.src = cms.InputTag("goodOfflinePrimaryVertices");
process.primaryVertexFilter.cut = cms.string(" ");

##-------- Muon events of interest --------

process.HLTMu = cms.EDFilter("HLTHighLevel",
                              TriggerResultsTag  = cms.InputTag("TriggerResults","","HLT"),
                              HLTPaths           = cms.vstring(),
                              eventSetupPathsKey = cms.string(''),
                              andOr = cms.bool(True), #----- True = OR, False = AND between the HLTPaths
                              throw = cms.bool(False) # throw exception on unknown path names
                            )

for path in options.hltPath.split(",") :
    
    process.HLTMu.HLTPaths.append(path.replace(' ',''))
    
#######################################
### MET shift correction on the fly ###
#######################################

from ElectroWeakAnalysis.VPlusJets.patMETSysShiftCorrection_cfi import *

runCorrection = cms.string("Run2012ABCD")
vertexCollection = cms.InputTag("goodOfflinePrimaryVertices")

metShiftSystematicCorrection(process,
                             options.isMC,
                             options.useSmearedCollection,
                             runCorrection,
                             vertexCollection,
                             options.didPhotonSmearing,
                             options.didTauSmearing)


##---------  W-->munu Collection ------------
#process.load("ElectroWeakAnalysis.VPlusJets.WmunuCollectionsPAT_cfi")

##---------  Jet Collection ----------------
#process.load("ElectroWeakAnalysis.VPlusJets.JetCollectionsPAT_cfi")

#
##########################################
## Filter to require at least two jets in the event
#process.RequireTwoJetsORboostedV = cms.EDFilter("JetsORboostedV",
#    minNumber = cms.untracked.int32(2),
#    maxNumber = cms.untracked.int32(100),
#    srcJets = cms.InputTag("ak5PFJetsLooseIdAll"),
#    srcVectorBoson = cms.InputTag("bestWmunu"),
#    minVpt = cms.untracked.double(100.),
#    minNumberPhotons = cms.untracked.int32(0)
#)
#process.RequireTwoJetsORboostedVStep = process.AllPassFilter.clone()


##-------- Save V+jets trees --------
#process.VplusJets = cms.EDAnalyzer("VplusJetsAnalysis", 
#    jetType = cms.string("PF"),
  #  srcPFCor = cms.InputTag("selectedPatJetsPFlow"),
#    srcPFCor = cms.InputTag("ak5PFJetsLooseId"),
#    srcPhoton = cms.InputTag("photons"),
#    IsoValPhoton = cms.VInputTag(cms.InputTag('phoPFIso:chIsoForGsfEle'),
#                                 cms.InputTag('phoPFIso:phIsoForGsfEle'),
#                                 cms.InputTag('phoPFIso:nhIsoForGsfEle'),
#                                                           ),
#    srcPFCorVBFTag = cms.InputTag("ak5PFJetsLooseIdVBFTag"), 
#    srcVectorBoson = cms.InputTag("bestWmunu"),
#    VBosonType     = cms.string('W'),
#    LeptonType     = cms.string('muon'),                          
#    TreeName    = cms.string('WJet'),
#    srcPrimaryVertex = cms.InputTag("goodOfflinePrimaryVertices"),                               
#    runningOverMC = cms.bool(isMC),			
#    runningOverAOD = cms.bool(False),			
#    srcMet = cms.InputTag("patType1CorrectedPFMet"),
#    srcMet = cms.InputTag("patMetShiftCorrected"), # type1 + shift corrections
#    srcMetMVA = cms.InputTag("pfMEtMVA"),
#    srcMetMVA = cms.InputTag("tcMet"),  
#    srcGen  = cms.InputTag("ak5GenJets"),
#    srcMuons  = cms.InputTag("selectedPatMuonsPFlow"),
#    srcBeamSpot  = cms.InputTag("offlineBeamSpot"),
#    srcCaloMet  = cms.InputTag("patMETs"),
#    srcgenMet  = cms.InputTag("genMetTrue"),
#    srcGenParticles  = cms.InputTag("genParticles"),
#    srcTcMet    = cms.InputTag("patMETsAK5TC"),
#    srcJetsforRho = cms.string("kt6PFJetsPFlow"),                               
#    srcJetsforRho_lepIso = cms.string("kt6PFJetsForIsolation"),       
#    srcJetsforRhoCHS = cms.string("kt6PFJetsChsPFlow"),
#    srcJetsforRho_lepIsoCHS = cms.string("kt6PFJetsChsForIsolationPFlow"),
#    srcFlavorByValue = cms.InputTag("ak5tagJet"),
#    bTagger=cms.string("simpleSecondaryVertexHighEffBJetTags"),
#    applyJECToGroomedJets=cms.bool(True),
#    doGroomedAK5 = cms.bool(True),
 #   doGroomedAK6 = cms.bool(True),
#    doGroomedAK7 = cms.bool(True),
#    doGroomedAK8 = cms.bool(True),
#    doGroomedAK10 = cms.bool(True),
#    doGroomedAK12 = cms.bool(True),
#    doGroomedCA8 = cms.bool(False),
#    doGroomedCA12 = cms.bool(False)
#)

#if isMC:
#    process.VplusJets.JEC_GlobalTag_forGroomedJet = cms.string("START53_V15")
#else:
#    process.VplusJets.JEC_GlobalTag_forGroomedJet = cms.string("FT_53_V10_AN3")


process.myseq = cms.Sequence( process.TrackVtxPath *
                              process.HLTMu*
                              process.metShiftSystematicCorrectionSequence)
#                              process.pfMEtSysShiftCorrSequence *
#                              process.patMetShiftCorrected )
#    process.WPath *
#    process.GenJetPath *
##    process.btagging * 
#    process.TagJetPath *
#    process.PFJetPath *
#    process.RequireTwoJetsORboostedV *
#    process.RequireTwoJetsORboostedVStep
#    )

if options.isMC:
     process.myseq.remove ( process.HLTMu)
#else:
#    process.myseq.remove ( process.GenJetPath)
#    process.myseq.remove ( process.TagJetPath)

if not options.runMetFilters :
    process.myseq.remove(process.TrackVtxPath)

process.p = cms.Path( process.myseq)
#* process.VplusJets)

############################
## Dump the output Python ##
############################

processDumpFile = open('processDump.py', 'w')
print >> processDumpFile, process.dumpPython()
   
