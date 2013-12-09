import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
import pprint
#CMSSW_BASE =os.environment()
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

options.register ('globalTag', 'GR_R_53_V10::All', VarParsing.multiplicity.singleton, VarParsing.varType.string,
                  "Global Tag to be used")

options.register ('numEventsToRun', 10 , VarParsing.multiplicity.singleton, VarParsing.varType.int,
                  "Number of events to process: -1 means all the input events")

options.register ('outputFileName', 'WmunuJetAnalysisntuple.root', VarParsing.multiplicity.singleton, VarParsing.varType.string,
                  "Name of the output file")

options.register ('reportEvery', 500, VarParsing.multiplicity.singleton, VarParsing.varType.int,
                  "Number of events after which print the report")

options.register ('runMetFilters', False, VarParsing.multiplicity.singleton, VarParsing.varType.int,
                  "Run the python/TrackingCollections_cfi to apply met filters and primary vertex collection skim")

options.register ('useSmearedCollection', True, VarParsing.multiplicity.singleton, VarParsing.varType.int,
                  "In case of MC analysis use smearedJets and smearedMet as default")

options.register ('hltPath', 'HLT_IsoMu24_*, HLT_Mu40_eta2p1*', VarParsing.multiplicity.singleton, VarParsing.varType.string,
                  "List of HLT path to be required when running on data. Should be separetad by ")

options.register ('didPhotonSystematic', True, VarParsing.multiplicity.singleton, VarParsing.varType.int,
                  "true if in PAT the photon smearing for resolution has been done")

options.register ('didTauSystematic', False, VarParsing.multiplicity.singleton, VarParsing.varType.int,
                  "true if in PAT the tau smearing for resolution has been done")

options.register ('isPileUpJetID', True, VarParsing.multiplicity.singleton, VarParsing.varType.int,
                  "true if you want to skim jet collection for the loose pile Up jet ID")

options.register ('isRequireTwoJets', False, VarParsing.multiplicity.singleton, VarParsing.varType.int,
                  "true if you want to do the two jet analysis --> resolved jet, non boosted category")

options.register ('storeSmearandShiftCollections', False, VarParsing.multiplicity.singleton, VarParsing.varType.int,
                  "true if you want to store and skim also shifted/smeared jet collections")

options.register ('skipAnalyzerAndDumpOutput', False , VarParsing.multiplicity.singleton, VarParsing.varType.int,
                  "true if you don't want to run the analyzer but dump a output file with all the collections keep*")


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
    if options.globalTag !='GR_R_53_V10::All' :  #inserted
     process.GlobalTag.globaltag = cms.string(globalTag)
    else :
     process.GlobalTag.globaltag = 'GR_P_V39_AN3::All'        
else:
    if options.globalTag !='GR_R_53_V10::All' :  #inserted
     process.GlobalTag.globaltag = cms.string(globalTag)
    else:
     process.GlobalTag.globaltag = 'START53_V7E::All'

###############################
### Input Module Definition ###
###############################

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(''))

if options.isMC:
    process.source.fileNames = cms.untracked.vstring('file:/uscms_data/d3/bmahakud/NewRaffaeleBibhuWJet/CMSSW_5_3_8_patch1/src/ElectroWeakAnalysis/VPlusJets/test/patuple_2_ntuple/tlbsm_53x_v3_mc.root')
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

if not options.skipAnalyzerAndDumpOutput :
  process.TFileService = cms.Service("TFileService", 
                                      fileName      = cms.string( options.outputFileName ),
                                      closeFileFast = cms.untracked.bool(False))
else:

  process.output = cms.OutputModule( "PoolOutputModule",
                                      fileName = cms.untracked.string(options.outputFileName),
                                      SelectEvents   = cms.untracked.PSet(SelectEvents = cms.vstring('p')),
                                      outputCommands = cms.untracked.vstring('keep *_*_*_*PAT*',
                                                                             'keep *_*_*_*demo*',
                                                                             'keep *_*_*_*SIM*'),
                                      dropMetaData   = cms.untracked.string('ALL'))


## --- to measure all the evetns processed
process.AllEventsStep = process.AllPassFilter.clone()
    
process.load("ElectroWeakAnalysis.VPlusJets.WmunuCollectionsPAT_cfi")
##---------  Vertex and track Collections -----------
process.load("ElectroWeakAnalysis.VPlusJets.TrackCollections_cfi")

#process.load("ElectroWeakAnalysis.VPlusJets.B2GJetCollectionsPAT_cfi")

process.primaryVertexFilter.src = cms.InputTag("goodOfflinePrimaryVertices");
process.primaryVertexFilter.cut = cms.string(" ");

process.MetFilterStep = process.AllPassFilter.clone()


##-------- Muon events of interest --------

process.HLTMu = cms.EDFilter("HLTHighLevel",
                              TriggerResultsTag  = cms.InputTag("TriggerResults","","HLT"),
                              HLTPaths           = cms.vstring(),
                              eventSetupPathsKey = cms.string(''),
                              andOr = cms.bool(True), #----- True = OR, False = AND between the HLTPaths
                              throw = cms.bool(False) # throw exception on unknown path names
                            )

process.HLTMuFilterStep = process.AllPassFilter.clone()

for path in options.hltPath.split(",") :    
    process.HLTMu.HLTPaths.append(path.replace(' ',''))
    
#######################################
### MET shift correction on the fly ###
#######################################

from ElectroWeakAnalysis.VPlusJets.patMETSysShiftCorrection_cfi import *

runCorrection                  = cms.string("Run2012ABCD")
vertexCollection               = cms.InputTag("goodOfflinePrimaryVertices")

patTypeIMetCorrected           = ['patMETsPFlow']
patTypeIMetCorrectedForMetUncertainty  = []

patTypeIMetCorrectedShifted = []
patTypeIMetCorrectedShiftedForMetUncertainty = []

metShiftSystematicCorrection(process,
                             options.isMC,
                             options.useSmearedCollection,
                             runCorrection,
                             vertexCollection,
                             options.didPhotonSystematic,
                             options.didTauSystematic,
                             patTypeIMetCorrected,
                             patTypeIMetCorrectedForMetUncertainty,
                             patTypeIMetCorrectedShifted,
                             patTypeIMetCorrectedShiftedForMetUncertainty)


####################################################
### Lepton Step : Muon Selection and Lepton Veto ### 
####################################################

##--------- W-->munu Collection ------------
from ElectroWeakAnalysis.VPlusJets.WmunuCollectionsPAT_cfi import*

patMuonCollection     = cms.InputTag("selectedPatMuonsPFlow")
patElectronCollection = cms.InputTag("selectedPatElectronsPFlow")

if options.isHEEPID:
    pTCutValue = 50.
    pTCutLooseMuonVeto = 20.
    pTCutLooseElectronVeto = 20.

else :
    pTCutValue = 25.
    pTCutLooseMuonVeto = 10.
    pTCutLooseElectronVeto = 10.
    
TransverseMassCutValue = 30.

WmunuCollectionsPAT(process,
                    patMuonCollection ,
                    patElectronCollection,
                    vertexCollection,
                    options.isQCD,
                    options.isHEEPID,
                    pTCutValue,
                    pTCutLooseMuonVeto,
                    pTCutLooseElectronVeto,
                    options.isTransverseMassCut,
                    TransverseMassCutValue,
                    patTypeIMetCorrectedShifted)


##############################################
### Standard AK5 jet collection selection  ### 
##############################################

from ElectroWeakAnalysis.VPlusJets.AK5JetCollectionsPATSelection_cfi import *

patJetCollection        = []
patJetCollection.append('selectedPatJetsPFlow')
patJetCollection.append('shiftedPatJetsPFlowEnUp')
patJetCollection.append('shiftedPatJetsPFlowEnDown')

patSmearedJetCollection = []
patSmearedJetCollection.append('smearedPatJetsPFlow')
patSmearedJetCollection.append('shiftedSmearedPatJetsPFlowEnUp')
patSmearedJetCollection.append('shiftedSmearedPatJetsPFlowEnDown')
patSmearedJetCollection.append('smearedPatJetsPFlowResDown')
patSmearedJetCollection.append('smearedPatJetsPFlowResUp')

jetPtThreshold          = 30.
useMVAPileUpJetID       = True

patCentralAK5Jets = []
patForwardAK5Jets = []

AK5JetCollectionsPATSelection(process,
                              options.isMC,
                              patJetCollection,
                              patSmearedJetCollection,
                              options.isPileUpJetID,
                              useMVAPileUpJetID,
                              options.useSmearedCollection,
                              options.storeSmearandShiftCollections,
                              jetPtThreshold,
                              options.isRequireTwoJets,
                              patCentralAK5Jets,
                              patForwardAK5Jets)

#############################################################################################################################
### Filter to require or at least two jets or one with hight pT --> preselection common to boosted and unboosted category ###    
#############################################################################################################################

process.RequireTwoJetsORboostedV = cms.EDFilter("JetsORboostedV",
                                                 minNumber = cms.untracked.int32(2),
                                                 maxNumber = cms.untracked.int32(100),
                                                 srcJets = cms.InputTag("ak5PFJetsPtSkimmed"),
                                                 srcVectorBoson = cms.InputTag("bestWmunu"),
                                                 srcPhotons = cms.InputTag("cleanPatPhotons"),
                                                 minVpt = cms.untracked.double(100.),
                                                 minNumberPhotons = cms.untracked.int32(0))

process.RequireTwoJetsORboostedVStep = process.AllPassFilter.clone()
print "BibhuprasadMahakud"
print patCentralAK5Jets
print patForwardAK5Jets

##-------- Save V+jets trees --------

process.VplusJets = cms.EDAnalyzer("VplusJetsAnalysis", 
    TreeName       = cms.string('WJet'),
    JetCollection  = cms.string('ak5PFJetsPtSkimmed'),
    LeptonType     = cms.string('muon'),                          
    runningOverMC  = cms.bool(options.isMC),			
    runningOverAOD = cms.bool(False),
    srcPFCorVBFTag = cms.InputTag("ak5PFJetsPtSkimmedForward"),
    srcPFCor = cms.InputTag("ak5PFJetsPtSkimmedCentral"),			
    useSmearedCollection          = cms.bool(options.useSmearedCollection),
    storeSmearandShiftCollections = cms.bool(options.storeSmearandShiftCollections),
    srcMuons  = cms.InputTag("selectedPatMuonsPFlow"), ## the real muon used for fillin the info is taken from Vboson as daughter
    jetType   = cms.string("PF"),
    srcPFPatCentralAK5Jets = cms.VInputTag(patCentralAK5Jets),
    srcPFPatForwardAK5Jets = cms.VInputTag(patForwardAK5Jets), 
    srcGenJets     = cms.InputTag("ak5GenJets"),
    srcJetsforRho  = cms.string("kt6PFJetsPFlow"),                               
    srcJetsforRho_lepIso = cms.string("kt6PFJetsForIsolation"),       
    srcFlavorByValue     = cms.InputTag("ak5tagJet"),
    VBosonType     = cms.string('W'),
    srcVectorBoson = cms.InputTag("bestWmunu"),
    srcCorrectedMet = cms.InputTag('patTypeIMetCorrectedShifted'),
    srcMetForUncertainty = cms.VInputTag(patTypeIMetCorrectedShiftedForMetUncertainty), 
    srcGenMet    = cms.InputTag("genMetTrue"),
    srcPhoton    = cms.InputTag("cleanPatPhotons"),
    IsoValPhoton = cms.VInputTag(cms.InputTag('phoPFIso:chIsoForGsfEle'),
                                 cms.InputTag('phoPFIso:phIsoForGsfEle'),
                                 cms.InputTag('phoPFIso:nhIsoForGsfEle')),
    srcPrimaryVertex = cms.InputTag("goodOfflinePrimaryVertices"),            
    srcGen  = cms.string('ak5GenJets'),
  #  srcMet = cms.InputTag("patMetShiftCorrected"),
    srcMet = cms.InputTag("pfType1CorrectedMet"),
    srcBeamSpot      = cms.InputTag("offlineBeamSpot"),
    srcGenParticles  = cms.InputTag("genParticles"),
    bTagger          = cms.string("combinedSecondaryVertexBJetTags"),
    pfCandidateForGrooming = cms.InputTag("pfNoElectronPFlow"),
    pfCandiatePileUp       = cms.InputTag("pfPileUpPFlow"),
 #   bTagger=cms.string("simpleSecondaryVertexHighEffBJetTags"),
    doGroomedAK3     = cms.bool(True),
    doGroomedAK35    = cms.bool(True),
    doGroomedAK4     = cms.bool(True),
    doGroomedAK45    = cms.bool(True),
    doGroomedAK5     = cms.bool(True),
    doGroomedAK55    = cms.bool(True),
    doGroomedAK6     = cms.bool(True),
    doGroomedAK65    = cms.bool(True),
    doGroomedAK7     = cms.bool(True),
    doGroomedAK75    = cms.bool(True),
    doGroomedAK8     = cms.bool(True),
    doGroomedAK85    = cms.bool(True),
    doGroomedAK9     = cms.bool(True),
    doGroomedAK95     = cms.bool(True),
    doGroomedAK10    = cms.bool(True),
    doGroomedAK105    = cms.bool(True),
    doGroomedAK11    = cms.bool(True),
    doGroomedAK115    = cms.bool(True),
    doGroomedAK12    = cms.bool(True),
#    doGroomedCA8     = cms.bool(True),
#    doGroomedCA12    = cms.bool(True),
    applyJECToGroomedJets =cms.bool(True)
)

 


if options.isMC:
     process.VplusJets.JEC_GlobalTag_forGroomedJet = cms.string("START53_V7G")  
else:
     process.VplusJets.JEC_GlobalTag_forGroomedJet = cms.string("FT_53_V10_AN3")


process.myseq = cms.Sequence( process.AllEventsStep*
                              process.TrackVtxPath *
                              process.MetFilterStep*
                              process.HLTMu*
                              process.HLTMuFilterStep*                              
                              process.metShiftSystematicCorrectionSequence*
                              process.WPath*
                              process.GenJetPath*
                              process.genTagJetPath*
                              process.ak5PFJetPath*
                              process.RequireTwoJetsORboostedV*
                              process.RequireTwoJetsORboostedVStep)                             


if options.isMC:
     process.myseq.remove(process.HLTMu)
     process.myseq.remove(process.HLTMuFilterStep)
else:
     process.myseq.remove(process.GenJetPath)
     process.myseq.remove(process.genTagJetPath)

if not options.runMetFilters :
    process.myseq.remove(process.TrackVtxPath)
    process.myseq.remove(process.MetFilterStep)

if options.skipAnalyzerAndDumpOutput :
    process.myseq.remove(process.AllEventsStep)
    process.myseq.remove(process.MetFilterStep)
    process.myseq.remove(process.HLTMuFilterStep)
    process.myseq.remove(process.RequireTwoJetsORboostedVStep)
    process.ak5PFJetPath.remove(process.RequireTwoJetsStep)
    process.WSequence.remove(process.tightLeptonStep)
    process.WSequence.remove(process.bestWToLepnuStep)
    process.VetoSequence.remove(process.looseMuonStep)
    process.VetoSequence.remove(process.looseElectronStep)
    
process.p = cms.Path( process.myseq  )

if options.skipAnalyzerAndDumpOutput :
 process.EndPath = cms.EndPath(process.output)
else:
    process.p += process.VplusJets
############################
## Dump the output Python ##
############################

processDumpFile = open('processDump.py', 'w')
print >> processDumpFile, process.dumpPython()
   
