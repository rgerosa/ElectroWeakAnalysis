import FWCore.ParameterSet.Config as cms

import RecoMET.METProducers.METSigParams_cfi as jetResolutions

from PhysicsTools.PatUtils.patPFMETCorrections_cff import *

from PhysicsTools.PatAlgos.producersLayer1.metProducer_cfi import patMETs


def runMetUncertainty(process,
                      useData,
                      doSmearing,
                      patJetCollection,
                      patMuonCollection,
                      patElectronCollection,
                      patTauCollection,
                      patPhotonCollection,
                      patUnclusteredCandidate,
                      TypeIpatMetcollection,
                      doTauSmearing,
                      doPhotonSmearing
                      ):

 print "                                                               "
 print "                                                               "
 print "###############################################################"
 print "################# metUncertainty Module #######################"
 print "###############################################################"
 print "                                                               "
 print "Options used : " 
 print "                                                               "
 print "useData                 = %d"%useData
 print "doSmearing              = %d"%doSmearing
 print "patJetCollection        = %s"%patJetCollection
 print "patElectronCollection   = %s"%patElectronCollection
 print "patMuonCollection       = %s"%patMuonCollection
 print "patTauCollection        = %s"%patTauCollection
 print "patPhotonCollection     = %s"%patPhotonCollection
 print "patUnclusteredCandidate = %s"%patUnclusteredCandidate
 print "TypeIpatMetcollection   = %s"%TypeIpatMetcollection


 ##### Additional Info for MC Jet Smearing  
 resolutionFile = cms.FileInPath('PhysicsTools/PatUtils/data/pfJetResolutionMCtoDataCorrLUT.root') ### file for Jet Resolution Smearing on the MC
 histoName      = cms.string('pfJetResolutionMCtoDataCorrLUT')  ## name of the histogram inside the previous file

 ##### Additional Info for JetEnergy Collection
 JECFileForJetShift   = cms.FileInPath('PhysicsTools/PatUtils/data/Summer12_V2_DATA_AK5PF_UncertaintySources.txt') ### JEC for data in txt format

 #### electron energy scale according to the HEEP prescription: 0.6% in EB and 1.5% in EE 
 EleEnergyScaleUncEB = cms.double(0.006)
 EleEnergyScaleUncEE = cms.double(0.015)

 #### electron energy resolution according to the HEEP prescription
 EleEnergyResEB      = cms.double(0.018)
 EleEnergyResEE      = cms.double(0.030)
 EleEnergyResUncEB   = cms.double(0.003)
 EleEnergyResUncEE   = cms.double(0.006)
 
 #### muon energy scale uncertainty of 0.2% as stated here:
 #### https://twiki.cern.ch/twiki/bin/viewauth/CMS/MuonReferenceResolution
 MuonEnergyScaleUnc   = cms.double(0.002)

 #### muon energy resolution of 0.6% as stated here:
 #### https://twiki.cern.ch/twiki/bin/viewauth/CMS/MuonReferenceResolution
 MuonEnergyRes      = cms.double(0.006)
 MuonEnergyResUnc   = cms.double(0.001)

 #### tau energy scale of 3% as stated in the metUncertainty tool
 TauEnergyScaleUnc = cms.double(0.030)

 ### tau energy resolution
 TauEnergyRes    = cms.double(0.030)
 TauEnergyResUnc = cms.double(0.010)

 #### Photon energy scale of 
 PhotonEnergyScaleUncEB = cms.double(0.010)
 PhotonEnergyScaleUncEE = cms.double(0.015)

 ### tau energy resolution
 PhotonEnergyResEB    = cms.double(0.018)
 PhotonEnergyResEE    = cms.double(0.030)
 PhotonEnergyResUncEB = cms.double(0.003)
 PhotonEnergyResUncEE = cms.double(0.006)

 #### unclustered energy resolution  (pfCandidate no in jet and jet pT<10 GeV) of 10% as stated in the metUncertainty tool
 UnclusteredEnergyUnc    = cms.double(0.100)

 print "                                "
 print "resolutionFile         = %s"%resolutionFile
 print "histoName              = %s"%histoName
 print "JECFileForJetShift     = %s"%JECFileForJetShift
 print "EleEnergyScaleUncEB    = %s"%EleEnergyScaleUncEB
 print "EleEnergyScaleUncEE    = %s"%EleEnergyScaleUncEE
 print "EleEnergyResEB         = %s"%EleEnergyResEB
 print "EleEnergyResEE         = %s"%EleEnergyResEE
 print "EleEnergyResUncEB      = %s"%EleEnergyResUncEB
 print "EleEnergyResUncEE      = %s"%EleEnergyResUncEE
 print "MuonEnergyScaleUnc     = %s"%MuonEnergyScaleUnc
 print "MuonEnergyRes          = %s"%MuonEnergyRes
 print "MuonEnergyResUnc       = %s"%MuonEnergyResUnc
 print "TauEnergyScaleUnc      = %s"%TauEnergyScaleUnc
 print "TauEnergyRes           = %s"%TauEnergyRes
 print "TauEnergyResUnc        = %s"%TauEnergyResUnc
 print "PhotonEnergyScaleUncEB = %s"%PhotonEnergyScaleUncEB
 print "PhotonEnergyScaleUncEE = %s"%PhotonEnergyScaleUncEE
 print "PhotonEnergyResEB      = %s"%PhotonEnergyResEB
 print "PhotonEnergyResEE      = %s"%PhotonEnergyResEE
 print "PhotonEnergyResUncEB   = %s"%PhotonEnergyResUncEB
 print "PhotonEnergyResUncEE   = %s"%PhotonEnergyResUncEE
 print "UnclusteredEnergyUnc   = %s"%UnclusteredEnergyUnc
 print "                                "

 if useData :
     print "                                "
     print "Error in parsing parameters --> useData cannot be true"
     print "MET Uncertainty Tool will not be run"
     print "                                "
     return ;
     
 ########## Smeared Jet Collection Producer taking into account different resolution between data and MC --> only for MC
 process.smearedPatJetsPFlow = cms.EDProducer("SmearedPATJetProducer",
       src = patJetCollection,
       skipCorrJetPtThreshold = cms.double(0.01),
       skipRawJetPtThreshold = cms.double(10.0),
       shiftBy = cms.double(0),
       inputFileName = resolutionFile,
       skipJetSelection = cms.string('jecSetsAvailable & abs(energy - correctedP4("Uncorrected").energy) > (5.*min(energy, correctedP4("Uncorrected").energy))'),
       dRmaxGenJetMatch = cms.string('TMath::Min(0.5, 0.1 + 0.3*TMath::Exp(-0.05*(genJetPt - 10.)))'),
       lutName = histoName,
       sigmaMaxGenJetMatch = cms.double(5.0),
       jetResolutions = cms.PSet(
                         ptresolthreshold = cms.double(10.0),
                         resolutionsEra = cms.string('Spring10'),
                         resolutionsAlgo = cms.string('AK5PF'),
                         jdpt0 = cms.vdouble(0.749, 0.829, 1.099, 1.355, 1.584, 1.807, 2.035, 2.217, 2.378, 2.591),
                         jdpt1 = cms.vdouble(0.718, 0.813, 1.133, 1.384, 1.588, 1.841, 2.115, 2.379, 2.508, 2.772),
                         jdpt2 = cms.vdouble(0.841, 0.937, 1.316, 1.605, 1.919, 2.295, 2.562, 2.722, 2.943, 3.293),
                         jdpt3 = cms.vdouble(0.929, 1.04, 1.46, 1.74, 2.042, 2.289, 2.639, 2.837, 2.946, 2.971),
                         jdpt4 = cms.vdouble(0.85, 0.961, 1.337, 1.593, 1.854, 2.005, 2.209, 2.533, 2.812, 3.047),
                         jdpt5 = cms.vdouble(1.049, 1.149, 1.607, 1.869, 2.012, 2.219, 2.289, 2.412, 2.695, 2.865),
                         jdpt6 = cms.vdouble(1.213, 1.298, 1.716, 2.015, 2.191, 2.612, 2.863, 2.879, 2.925, 2.902),
                         jdpt7 = cms.vdouble(1.094, 1.139, 1.436, 1.672, 1.831, 2.05, 2.267, 2.549, 2.785, 2.86),
                         jdpt8 = cms.vdouble(0.889, 0.939, 1.166, 1.365, 1.553, 1.805, 2.06, 2.22, 2.268, 2.247),
                         jdpt9 = cms.vdouble(0.843, 0.885, 1.245, 1.665, 1.944, 1.981, 1.972, 2.875, 3.923, 7.51),
                         jdphi0 = cms.vdouble(0.034, 0.034, 0.034, 0.034, 0.032, 0.031, 0.028, 0.027, 0.027, 0.027),
                         jdphi1 = cms.vdouble(0.034, 0.035, 0.035, 0.035, 0.035, 0.034, 0.031, 0.03, 0.029, 0.027),
                         jdphi2 = cms.vdouble(0.04, 0.04, 0.04, 0.04, 0.04,0.038, 0.036, 0.035, 0.034, 0.033),
                         jdphi3 = cms.vdouble(0.042, 0.043, 0.044, 0.043, 0.041, 0.039, 0.039, 0.036, 0.034, 0.031),
                         jdphi4 = cms.vdouble(0.042, 0.042, 0.043, 0.042, 0.038, 0.036, 0.036, 0.033, 0.031, 0.031),
                         jdphi5 = cms.vdouble(0.069, 0.069, 0.064, 0.058, 0.053,0.049, 0.049, 0.043, 0.039, 0.04),
                         jdphi6 = cms.vdouble(0.084, 0.08, 0.072, 0.065, 0.066, 0.06, 0.051, 0.049, 0.045, 0.045),
                         jdphi7 = cms.vdouble(0.077, 0.072, 0.059, 0.05, 0.045, 0.042, 0.039, 0.039, 0.037, 0.031),
                         jdphi8 = cms.vdouble(0.059, 0.057, 0.051, 0.044, 0.038, 0.035, 0.037, 0.032, 0.028, 0.028),
                         jdphi9 = cms.vdouble(0.062, 0.059, 0.053, 0.047, 0.042, 0.045, 0.036, 0.032, 0.034, 0.044),
                         HB_EtResPar = cms.vdouble(0.0, 1.22, 0.05), 
                         HE_EtResPar = cms.vdouble(0.0, 1.3, 0.05), 
                         EB_EtResPar = cms.vdouble(0.2, 0.03, 0.005),
                         EE_EtResPar = cms.vdouble(0.2, 0.03, 0.005),
                         HF_EtResPar = cms.vdouble(0.0, 1.82, 0.09),
                         HO_EtResPar = cms.vdouble(0.0, 1.3, 0.005),
                         PF_EtResType1 = cms.vdouble(0.05, 0, 0),
                         PF_EtResType2 = cms.vdouble(0.05, 0, 0),
                         PF_EtResType3 = cms.vdouble(0.05, 0, 0),
                         PF_EtResType4 = cms.vdouble(0.042, 0.1, 0.0),
                         PF_EtResType5 = cms.vdouble(0.41, 0.52, 0.25),
                         PF_EtResType6 = cms.vdouble(0.0, 1.22, 0.05),
                         PF_EtResType7 = cms.vdouble(0.0, 1.22, 0.05),
                         HB_PhiResPar = cms.vdouble(0.02511),
                         HE_PhiResPar = cms.vdouble(0.02511),
                         EB_PhiResPar = cms.vdouble(0.00502),
                         EE_PhiResPar = cms.vdouble(0.02511),
                         HF_PhiResPar = cms.vdouble(0.05022),
                         HO_PhiResPar = cms.vdouble(0.02511),
                         PF_PhiResType1 = cms.vdouble(0.002),
                         PF_PhiResType2 = cms.vdouble(0.002),
                         PF_PhiResType3 = cms.vdouble(0.002),
                         PF_PhiResType4 = cms.vdouble(0.0028, 0.0, 0.0022),
                         PF_PhiResType5 = cms.vdouble(0.1, 0.1, 0.13),
                         PF_PhiResType6 = cms.vdouble(0.02511),
                         PF_PhiResType7 = cms.vdouble(0.02511))
 )

 process.smearedPatJetsPFlowResDown = process.smearedPatJetsPFlow.clone(shiftBy = cms.double(-1.0))

 process.smearedPatJetsPFlowResUp = process.smearedPatJetsPFlow.clone(shiftBy = cms.double(1.0))

 if doSmearing and not useData:
     
  process.patJetSmearingSequence = cms.Sequence(process.smearedPatJetsPFlow*
                                                process.smearedPatJetsPFlowResDown*
                                                process.smearedPatJetsPFlowResUp)

 ########## Shift jets --> pay attention to the label                                                                           
 process.shiftedSmearedPatJetsPFlowEnUp = cms.EDProducer("ShiftedPATJetProducer",
                                                         addResidualJES = cms.bool(False),
                                                         src = cms.InputTag("smearedPatJetsPFlow"),
                                                         jetCorrInputFileName = JECFileForJetShift,
                                                         shiftBy = cms.double(1.0),
                                                         jetCorrUncertaintyTag = cms.string('SubTotalDataMC'),
                                                         jetCorrLabelUpToL3Res = cms.string('ak5PFL1FastL2L3Residual'),
                                                         jetCorrLabelUpToL3 = cms.string('ak5PFL1FastL2L3'))


 process.shiftedSmearedPatJetsPFlowEnDown  = process.shiftedSmearedPatJetsPFlowEnUp.clone(shiftBy = cms.double(-1.0))

 process.shiftedPatJetsPFlowEnUp           = process.shiftedSmearedPatJetsPFlowEnUp.clone(src = patJetCollection)

 process.shiftedPatJetsPFlowEnDown         = process.shiftedPatJetsPFlowEnUp.clone(shiftBy = cms.double(-1.0))

 if doSmearing and not useData :
     
       process.shiftJetsForMetUncertainty = cms.Sequence(process.shiftedSmearedPatJetsPFlowEnUp*
                                                         process.shiftedSmearedPatJetsPFlowEnDown*
                                                         process.shiftedPatJetsPFlowEnUp*
                                                         process.shiftedPatJetsPFlowEnDown)
 else:
       process.shiftJetsForMetUncertainty = cms.Sequence(process.shiftedPatJetsPFlowEnUp*
                                                         process.shiftedPatJetsPFlowEnDown)
       
 ########## Shift electrons 
 process.shiftedPatElectronsPFlowEnUp = cms.EDProducer("ShiftedPATElectronProducer",
                                                       src = patElectronCollection,
                                                       shiftBy = cms.double(1.0),
                                                       binning = cms.VPSet(cms.PSet(binUncertainty = EleEnergyScaleUncEB,
                                                                                    binSelection = cms.string('isEB')),
                                                                           cms.PSet(binUncertainty = EleEnergyScaleUncEE,
                                                                                    binSelection = cms.string('!isEB')))
                                                       )
   
 process.shiftedPatElectronsPFlowEnDown = process.shiftedPatElectronsPFlowEnUp.clone( shiftBy = cms.double(-1.0) )

 process.shiftElectronsForMetUncertainty = cms.Sequence(process.shiftedPatElectronsPFlowEnUp*
                                                        process.shiftedPatElectronsPFlowEnDown)

 ############ Smeared Electrons
 process.smearedPatElectronsPFlow = cms.EDProducer("SmearedPATElectronProducer",
                                                   src = patElectronCollection,
                                                   binning = cms.VPSet(cms.PSet(
                                                         binSelection = cms.string('isEB'),
                                                         binResolution = EleEnergyResEB,
                                                         binResolutionUncertainty = EleEnergyResUncEB),
                                                       cms.PSet(
                                                         binSelection = cms.string('isEB'),
                                                         binResolution = EleEnergyResEE,
                                                         binResolutionUncertainty = EleEnergyResUncEE)))

 if EleEnergyResUncEB !=0  or EleEnergyResUncEE !=0 :

     process.smearedPatElectronsPFlowResUp   = process.smearedPatElectronsPFlow.clone(shiftBy = cms.double(1.0))
     process.smearedPatElectronsPFlowResDown = process.smearedPatElectronsPFlow.clone(shiftBy = cms.double(-1.0))


 if doSmearing and not useData:
     
  process.smearedElectronsForMetUncertainty = cms.Sequence(process.smearedPatElectronsPFlow) 

  if EleEnergyResUncEB !=0  or EleEnergyResUncEE !=0 :

      process.smearedElectronsForMetUncertainty += process.smearedPatElectronsPFlowResUp*process.smearedPatElectronsPFlowResDown
   
 
 ########## Shift muons
 process.shiftedPatMuonsPFlowEnUp = cms.EDProducer("ShiftedPATMuonProducer",
                                                   src = patMuonCollection,
                                                   uncertainty = MuonEnergyScaleUnc,
                                                   shiftBy = cms.double(1.0))

 process.shiftedPatMuonsPFlowEnDown  = process.shiftedPatMuonsPFlowEnUp.clone(shiftBy = cms.double(-1.0))

 process.shiftMuonsForMetUncertainty = cms.Sequence(process.shiftedPatMuonsPFlowEnUp*
                                                    process.shiftedPatMuonsPFlowEnDown)
 
 ############ Smeared Muons
 process.smearedPatMuonsPFlow = cms.EDProducer("SmearedPATMuonProducer",
                                                src = patMuonCollection,
                                                resolution = MuonEnergyRes,
                                                resolutionUnc = MuonEnergyResUnc) 

 if MuonEnergyResUnc !=0  :
     
     process.smearedPatMuonsPFlowResUp   = process.smearedPatMuonsPFlow.clone(shiftBy = cms.double(1.0))
     process.smearedPatMuonsPFlowResDown = process.smearedPatMuonsPFlow.clone(shiftBy = cms.double(-1.0))

 if doSmearing and not useData:
     
  process.smearedMuonsForMetUncertainty = cms.Sequence(process.smearedPatMuonsPFlow) 

  if MuonEnergyResUnc !=0:

      process.smearedMuonsForMetUncertainty += process.smearedPatMuonsPFlowResUp*process.smearedPatMuonsPFlowResDown

 ############ Shift taus
 process.shiftedPatTausPFlowEnUp = cms.EDProducer("ShiftedPATTauProducer",
                                                   src = patTauCollection,
                                                   uncertainty = TauEnergyScaleUnc,
                                                   shiftBy = cms.double(1.0))

 process.shiftedPatTausPFlowEnDown  = process.shiftedPatTausPFlowEnUp.clone(shiftBy = cms.double(-1.0))

 process.shiftTausForMetUncertainty = cms.Sequence(process.shiftedPatTausPFlowEnUp*process.shiftedPatTausPFlowEnDown)

 ############ Smeared Taus
 process.smearedPatTausPFlow = cms.EDProducer("SmearedPATTauProducer",
                                               src = patTauCollection,
                                               resolution = TauEnergyRes,
                                               resolutionUnc = TauEnergyResUnc)

 if TauEnergyResUnc !=0  :
     
     process.smearedPatTausPFlowResUp   = process.smearedPatTausPFlow.clone(shiftBy = cms.double(1.0))
     process.smearedPatTausPFlowResDown = process.smearedPatTausPFlow.clone(shiftBy = cms.double(-1.0))

 if doSmearing and not useData: 

  process.smearedTausForMetUncertainty = cms.Sequence(process.smearedPatTausPFlow) 

  if TauEnergyResUnc !=0  :

   process.smearedTausForMetUncertainty += process.smearedPatTausPFlowResUp*process.smearedPatTausPFlowResDown

 ############ Shift photon
 process.shiftedPatPhotonPFlowEnUp = cms.EDProducer("ShiftedPATPhotonProducer",
                                                   src = patPhotonCollection,
                                                   binning = cms.VPSet(cms.PSet(binUncertainty = PhotonEnergyScaleUncEB,
                                                                                binSelection = cms.string('isEB')),
                                                                       cms.PSet(binUncertainty = PhotonEnergyScaleUncEE,
                                                                                binSelection = cms.string('!isEB'))),
                                                   shiftBy = cms.double(1.0))
  

 process.shiftedPatPhotonPFlowEnDown  = process.shiftedPatPhotonPFlowEnUp.clone(shiftBy = cms.double(-1.0))

 process.shiftPhotonForMetUncertainty = cms.Sequence(process.shiftedPatPhotonPFlowEnUp*process.shiftedPatPhotonPFlowEnDown)

 ############ Smeared photon
 process.smearedPatPhotonPFlow = cms.EDProducer("SmearedPATPhotonProducer",
                                                 src = patPhotonCollection,
                                                 binning = cms.VPSet(cms.PSet(
                                                         binSelection = cms.string('isEB'),
                                                         binResolution = PhotonEnergyResEB,
                                                         binResolutionUncertainty = PhotonEnergyResUncEB),
                                                       cms.PSet(
                                                         binSelection = cms.string('isEB'),
                                                         binResolution = PhotonEnergyResEE,
                                                         binResolutionUncertainty = PhotonEnergyResUncEE)))

 if PhotonEnergyResUncEB !=0 or PhotonEnergyResUncEE !=0  :
     
     process.smearedPatPhotonPFlowResUp   = process.smearedPatPhotonPFlow.clone(shiftBy = cms.double(1.0))
     process.smearedPatPhotonPFlowResDown = process.smearedPatPhotonPFlow.clone(shiftBy = cms.double(-1.0))

 if doSmearing and not useData: 

  process.smearedPhotonForMetUncertainty = cms.Sequence(process.smearedPatPhotonPFlow) 

  if PhotonEnergyResUncEB !=0  or PhotonEnergyResUncEE !=0 :

   process.smearedPhotonForMetUncertainty += process.smearedPatPhotonPFlowResUp*process.smearedPatPhotonPFlowResDown


 ##### raw met used for smearing procedure
 process.patMETforMetUncertainty = process.patMETs.clone(metSource = cms.InputTag("pfMet")) 

 ## correction to the met due to the jet smearing
 process.patPFMETcorrSmearedJet = cms.EDProducer("ShiftedParticleMETcorrInputProducer",
                                                  srcOriginal = patJetCollection,
                                                  srcShifted = cms.InputTag("smearedPatJetsPFlow"))

 process.patPFSmearedMet = cms.EDProducer("CorrectedPATMETProducer", 
                                           applyType2Corrections = cms.bool(False),
                                           srcType1Corrections = cms.VInputTag(cms.InputTag("patPFMETcorrSmearedJet")),
                                           applyType1Corrections = cms.bool(True),
                                           src = cms.InputTag("patMETforMetUncertainty"))

 process.patPFSmearedMetForMetUncertainty = cms.Sequence( process.patPFMETcorrSmearedJet*
                                                          process.patPFSmearedMet)

 ### derived the TypeI correction value for the smeared jet
 process.patPFSmearedJetMETtype1Corr  = cms.EDProducer("PATPFJetMETcorrInputProducer", 
                                                        src = cms.InputTag("smearedPatJetsPFlow"),
                                                        type1JetPtThreshold = cms.double(10.0),
                                                        skipEMfractionThreshold = cms.double(0.9),
                                                        skipEM = cms.bool(True),
                                                        offsetCorrLabel = cms.string('L1FastJet'),
                                                        skipMuons = cms.bool(True),
                                                        skipMuonSelection = cms.string('isGlobalMuon | isStandAloneMuon'),
                                                        jetCorrEtaMax = cms.double(9.9),
                                                        jetCorrLabel = cms.string('L3Absolute'))
 
 process.patType1CorrectedPFSmearedMet = cms.EDProducer("CorrectedPATMETProducer", ## should be chs independent
                                                         src = cms.InputTag("patPFSmearedMet"),
                                                         applyType1Corrections = cms.bool(True),
                                                         srcType1Corrections = cms.VInputTag(cms.InputTag("patPFSmearedJetMETtype1Corr","type1")),
                                                         srcCHSSums = cms.VInputTag(cms.InputTag("pfchsMETcorrPFlow","type0")),
                                                         applyType2Corrections = cms.bool(False),
                                                         type0Rsoft = cms.double(0.6),
                                                         applyType0Corrections = cms.bool(False))

 if doSmearing and not useData:
  process.patType1CorrectedPFSmearedMetForMetUncertainty = cms.Sequence(process.patPFSmearedJetMETtype1Corr*
                                                                        process.patType1CorrectedPFSmearedMet)

   
 ### met correction for jet up and down
 process.patPFMETcorrSmearedJetEnUp  = cms.EDProducer("ShiftedParticleMETcorrInputProducer",
                                                       srcOriginal = cms.InputTag("smearedPatJetsPFlow"),
                                                       srcShifted  = cms.InputTag("shiftedSmearedPatJetsPFlowEnUp"))

 process.patPFMETcorrSmearedJetEnDown = process.patPFMETcorrSmearedJetEnUp.clone(srcShifted = cms.InputTag("shiftedSmearedPatJetsPFlowEnDown"))

 process.patPFMETcorrJetEnUp          = process.patPFMETcorrSmearedJetEnUp.clone(srcOriginal = patJetCollection,
                                                                                 srcShifted = cms.InputTag("shiftedPatJetsPFlowEnUp"))

 process.patPFMETcorrJetEnDown        = process.patPFMETcorrJetEnUp.clone(srcShifted = cms.InputTag("shiftedPatJetsPFlowEnDown"))

 if doSmearing and not useData :
      process.patPFMETcorrJetForMetUncertainty = cms.Sequence( process.patPFMETcorrSmearedJetEnUp*
                                                               process.patPFMETcorrSmearedJetEnDown*
                                                               process.patPFMETcorrJetEnUp*
                                                               process.patPFMETcorrJetEnDown)
      
 else :
      process.patPFMETcorrJetForMetUncertainty = cms.Sequence( process.patPFMETcorrJetEnUp*
                                                               process.patPFMETcorrJetEnDown)
      
 #### propagation on the typeI met correction
 process.patType1CorrectedPFSmearedMetJetEnUp = cms.EDProducer("CorrectedPATMETProducer",
                                                                applyType2Corrections = cms.bool(False),
                                                                srcType1Corrections = cms.VInputTag(cms.InputTag("patPFMETcorrSmearedJetEnUp")),
                                                                applyType1Corrections = cms.bool(True),
                                                                src = cms.InputTag("patType1CorrectedPFSmearedMet"))

 process.patType1CorrectedPFSmearedMetJetEnDown = process.patType1CorrectedPFSmearedMetJetEnUp.clone(srcType1Corrections = cms.VInputTag(cms.InputTag("patPFMETcorrSmearedJetEnDown")))
  
 process.patType1CorrectedPFMetJetEnUp          = process.patType1CorrectedPFSmearedMetJetEnUp.clone( src = TypeIpatMetcollection,
                                                                                                      srcType1Corrections = cms.VInputTag(cms.InputTag("patPFMETcorrJetEnUp")))

 process.patType1CorrectedPFMetJetEnDown        = process.patType1CorrectedPFSmearedMetJetEnUp.clone( src = TypeIpatMetcollection,
                                                                                                      srcType1Corrections = cms.VInputTag(cms.InputTag("patPFMETcorrJetEnDown")))

 if doSmearing and not useData :
    process.patType1CorrectedPFMetJetForMetUncertainty = cms.Sequence(process.patType1CorrectedPFSmearedMetJetEnUp*
                                                                      process.patType1CorrectedPFSmearedMetJetEnDown*
                                                                      process.patType1CorrectedPFMetJetEnUp*
                                                                      process.patType1CorrectedPFMetJetEnDown)
    
 else:
    process.patType1CorrectedPFMetJetForMetUncertainty = cms.Sequence(process.patType1CorrectedPFMetJetEnUp*
                                                                      process.patType1CorrectedPFMetJetEnDown)

 ###### Resolution up and down effect on the met
 process.patPFMETcorrSmearedJetResUp = cms.EDProducer("ShiftedParticleMETcorrInputProducer",
                                                      srcOriginal = cms.InputTag("smearedPatJetsPFlow"),
                                                      srcShifted = cms.InputTag("smearedPatJetsPFlowResUp"))

 process.patPFMETcorrSmearedJetResDown = process.patPFMETcorrSmearedJetResUp.clone(srcShifted = cms.InputTag("smearedPatJetsPFlowResDown"))


 process.patType1CorrectedPFSmearedMetJetResUp = cms.EDProducer("CorrectedPATMETProducer",
                                                                 applyType2Corrections = cms.bool(False),
                                                                 srcType1Corrections = cms.VInputTag(cms.InputTag("patPFMETcorrSmearedJetResUp")),
                                                                 applyType1Corrections = cms.bool(True),
                                                                 src = cms.InputTag("patType1CorrectedPFSmearedMet"))

 process.patType1CorrectedPFSmearedMetJetResDown = process.patType1CorrectedPFSmearedMetJetResUp.clone(srcType1Corrections = cms.VInputTag(cms.InputTag("patPFMETcorrSmearedJetResDown")))

 if doSmearing and not useData:
    process.patType1CorrectedPFMetJetResForMetUncertainty =   cms.Sequence( process.patPFMETcorrSmearedJetResUp*
                                                                            process.patPFMETcorrSmearedJetResDown*
                                                                            process.patType1CorrectedPFSmearedMetJetResUp*
                                                                            process.patType1CorrectedPFSmearedMetJetResDown)


 ###### contribution due to unclustered energy
 ## correction due to pfParticle outside jets
    
 process.pfCandMETcorrUnclusteredEnUp = cms.EDProducer("ShiftedMETcorrInputProducer",
                                                        src = cms.VInputTag(patUnclusteredCandidate),
                                                        uncertainty = UnclusteredEnergyUnc,
                                                        shiftBy = cms.double(1.0))

 process.pfCandMETcorrUnclusteredEnDown = process.pfCandMETcorrUnclusteredEnUp.clone(shiftBy = cms.double(-1.0));

 ## correction due to low pT jets < 10 GeV
 process.patPFJetMETcorrUnclusteredEnUp = cms.EDProducer("ShiftedMETcorrInputProducer",
                                                          src = cms.VInputTag(cms.InputTag("pfJetMETcorrPFlow","type2"), cms.InputTag("pfJetMETcorrPFlow","offset")),
                                                          uncertainty = cms.double(0.1),
                                                          shiftBy = cms.double(1.0))

 process.patPFJetMETcorrUnclusteredEnDown = process.patPFJetMETcorrUnclusteredEnUp.clone(shiftBy = cms.double(-1.0))
 
 process.patType1CorrectedPFMetUnclusteredEnUp = cms.EDProducer("CorrectedPATMETProducer",
                                                                 applyType2Corrections = cms.bool(False),
                                                                 srcType1Corrections = cms.VInputTag(cms.InputTag("pfCandMETcorrUnclusteredEnUp"),
                                                                                                     cms.InputTag("patPFJetMETcorrUnclusteredEnUp","type2"),
                                                                                                     cms.InputTag("patPFJetMETcorrUnclusteredEnUp","offset"),
                                                                                                     cms.InputTag("patPFJetMETcorrUnclusteredEnUp","type2")),
                                                                 applyType1Corrections = cms.bool(True),
                                                                 src = TypeIpatMetcollection)

                                                                                            
 process.patType1CorrectedPFMetUnclusteredEnDown  =  process.patType1CorrectedPFMetUnclusteredEnUp.clone(
                                                                  srcType1Corrections = cms.VInputTag(cms.InputTag("pfCandMETcorrUnclusteredEnDown"),
                                                                                                      cms.InputTag("patPFJetMETcorrUnclusteredEnDown","type2"),
                                                                                                      cms.InputTag("patPFJetMETcorrUnclusteredEnDown","offset"),
                                                                                                      cms.InputTag("patPFJetMETcorrUnclusteredEnDown","type2")))

 process.patType1CorrectedPFSmearedMetUnclusteredEnUp    =  process.patType1CorrectedPFMetUnclusteredEnUp.clone(src = cms.InputTag("patType1CorrectedPFSmearedMet"))

 process.patType1CorrectedPFSmearedMetUnclusteredEnDown  =  process.patType1CorrectedPFMetUnclusteredEnDown.clone(src = cms.InputTag("patType1CorrectedPFSmearedMet"))

 if doSmearing and not useData :
      process.patType1CorrectedPFMetUnclusteredForMetUncertainty = cms.Sequence(process.pfCandMETcorrUnclusteredEnUp*
                                                                                process.pfCandMETcorrUnclusteredEnDown*
                                                                                process.patPFJetMETcorrUnclusteredEnUp*
                                                                                process.patPFJetMETcorrUnclusteredEnDown*
                                                                                process.patType1CorrectedPFMetUnclusteredEnUp*
                                                                                process.patType1CorrectedPFMetUnclusteredEnDown*
                                                                                process.patType1CorrectedPFSmearedMetUnclusteredEnUp*
                                                                                process.patType1CorrectedPFSmearedMetUnclusteredEnDown)
 else:
     process.patType1CorrectedPFMetUnclusteredForMetUncertainty = cms.Sequence(process.pfCandMETcorrUnclusteredEnUp*
                                                                               process.pfCandMETcorrUnclusteredEnDown*
                                                                               process.patPFJetMETcorrUnclusteredEnUp*
                                                                               process.patPFJetMETcorrUnclusteredEnDown*
                                                                               process.patType1CorrectedPFMetUnclusteredEnUp*
                                                                               process.patType1CorrectedPFMetUnclusteredEnDown)

      
 ############### electrons contribution to met 
 process.patPFMETcorrElectronEnUp = cms.EDProducer("ShiftedParticleMETcorrInputProducer",
                                                    srcOriginal = patElectronCollection,
                                                    srcShifted = cms.InputTag("shiftedPatElectronsPFlowEnUp"))

 process.patPFMETcorrElectronEnDown = process.patPFMETcorrElectronEnUp.clone(srcShifted = cms.InputTag("shiftedPatElectronsPFlowEnDown"))


 process.patType1CorrectedPFMetElectronEnUp = cms.EDProducer("CorrectedPATMETProducer",
                                                              applyType2Corrections = cms.bool(False),
                                                              srcType1Corrections = cms.VInputTag(cms.InputTag("patPFMETcorrElectronEnUp")),
                                                              src = TypeIpatMetcollection,
                                                              applyType1Corrections = cms.bool(True))

 process.patType1CorrectedPFMetElectronEnDown = process.patType1CorrectedPFMetElectronEnUp.clone(srcType1Corrections = cms.VInputTag(cms.InputTag("patPFMETcorrElectronEnDown")))
                                                                                          
 process.patType1CorrectedPFSmearedMetElectronEnUp    = process.patType1CorrectedPFMetElectronEnUp.clone(src = cms.InputTag("patType1CorrectedPFSmearedMet"))

 process.patType1CorrectedPFSmearedMetElectronEnDown  = process.patType1CorrectedPFMetElectronEnDown.clone(src = cms.InputTag("patType1CorrectedPFSmearedMet"))

 if doSmearing and not useData :
      process.patType1CorrectedPFMetElectronEnForMetUncertainty = cms.Sequence( process.patPFMETcorrElectronEnUp*
                                                                                process.patPFMETcorrElectronEnDown*
                                                                                process.patType1CorrectedPFMetElectronEnUp*
                                                                                process.patType1CorrectedPFMetElectronEnDown*
                                                                                process.patType1CorrectedPFSmearedMetElectronEnUp*
                                                                                process.patType1CorrectedPFSmearedMetElectronEnDown
                                                                              )  
 else:
      process.patType1CorrectedPFMetElectronEnForMetUncertainty = cms.Sequence( process.patPFMETcorrElectronEnUp*
                                                                                process.patPFMETcorrElectronEnDown*
                                                                                process.patType1CorrectedPFMetElectronEnUp*
                                                                                process.patType1CorrectedPFMetElectronEnDown
                                                                              )  

 ## electron energy resolution
 process.patPFMETcorrElectronRes = cms.EDProducer("ShiftedParticleMETcorrInputProducer",
                                                   srcOriginal = patElectronCollection,
                                                   srcShifted = cms.InputTag("smearedPatElectronsPFlow"))

      
 process.patPFMETcorrElectronResUp   = process.patPFMETcorrElectronRes.clone(srcShifted = cms.InputTag("smearedPatElectronsPFlowResUp"))
 process.patPFMETcorrElectronResDown = process.patPFMETcorrElectronRes.clone(srcShifted = cms.InputTag("smearedPatElectronsPFlowResDown"))

 if doSmearing and not useData:

   process.patPFMETcorrElectronResForMetUncertainty = cms.Sequence(process.patPFMETcorrElectronRes)

   if EleEnergyResUncEB !=0  or EleEnergyResUncEE !=0 :

    process.patPFMETcorrElectronResForMetUncertainty += process.patPFMETcorrElectronResUp*process.patPFMETcorrElectronResDown   


 process.patType1CorrectedPFMetElectronRes = cms.EDProducer("CorrectedPATMETProducer",
                                                             applyType2Corrections = cms.bool(False),
                                                             srcType1Corrections = cms.VInputTag(cms.InputTag("patPFMETcorrElectronRes")),
                                                             src = TypeIpatMetcollection,
                                                             applyType1Corrections = cms.bool(True))

 
 process.patType1CorrectedPFMetElectronResUp   = process.patType1CorrectedPFMetElectronRes.clone(srcType1Corrections = cms.VInputTag(cms.InputTag("patPFMETcorrElectronResUp")))

 process.patType1CorrectedPFMetElectronResDown = process.patType1CorrectedPFMetElectronRes.clone(srcType1Corrections = cms.VInputTag(cms.InputTag("patPFMETcorrElectronResDown")))

 process.patType1CorrectedPFSmearedMetElectronRes = process.patType1CorrectedPFMetElectronRes.clone(src = cms.InputTag("patType1CorrectedPFSmearedMet"))
                                                                              
 process.patType1CorrectedPFSmearedMetElectronResUp = process.patType1CorrectedPFSmearedMetElectronRes.clone(srcType1Corrections =
                                                                                                               cms.VInputTag(cms.InputTag("patPFMETcorrElectronResUp")))
 process.patType1CorrectedPFSmearedMetElectronResDown = process.patType1CorrectedPFSmearedMetElectronRes.clone(srcType1Corrections =
                                                                                                               cms.VInputTag(cms.InputTag("patPFMETcorrElectronResDown")))
   
 if doSmearing and not useData:
     process.patType1CorrectedPFMetElectronResForMetUncertainty = cms.Sequence(process.patType1CorrectedPFMetElectronRes*
                                                                               process.patType1CorrectedPFSmearedMetElectronRes)

     if EleEnergyResUncEB !=0  or EleEnergyResUncEE !=0 :
         process.patType1CorrectedPFMetElectronResForMetUncertainty += process.patType1CorrectedPFMetElectronResUp*process.patType1CorrectedPFMetElectronResDown
         process.patType1CorrectedPFMetElectronResForMetUncertainty += process.patType1CorrectedPFSmearedMetElectronResUp*process.patType1CorrectedPFSmearedMetElectronResDown

                  
 ######### muons contribution to met 
 process.patPFMETcorrMuonEnUp = cms.EDProducer("ShiftedParticleMETcorrInputProducer",
                                                srcOriginal = patMuonCollection,
                                                srcShifted = cms.InputTag("shiftedPatMuonsPFlowEnUp"))

 process.patPFMETcorrMuonEnDown = process.patPFMETcorrMuonEnUp.clone(srcShifted = cms.InputTag("shiftedPatMuonsPFlowEnDown"))

 process.patType1CorrectedPFSmearedMetMuonEnUp = cms.EDProducer("CorrectedPATMETProducer",
                                                      applyType2Corrections = cms.bool(False),
                                                      srcType1Corrections = cms.VInputTag(cms.InputTag("patPFMETcorrMuonEnUp")),
                                                      src = cms.InputTag("patType1CorrectedPFSmearedMet"),
                                                      applyType1Corrections = cms.bool(True))

 process.patType1CorrectedPFSmearedMetMuonEnDown  = process.patType1CorrectedPFSmearedMetMuonEnUp.clone(srcType1Corrections = cms.VInputTag(cms.InputTag("patPFMETcorrMuonEnDown")))
                                                                                          
 process.patType1CorrectedPFMetMuonEnUp           = process.patType1CorrectedPFSmearedMetMuonEnUp.clone(src = TypeIpatMetcollection)

 process.patType1CorrectedPFMetMuonEnDown         = process.patType1CorrectedPFSmearedMetMuonEnDown.clone(src = TypeIpatMetcollection)

 if doSmearing and not useData :
      process.patType1CorrectedPFMetMuonEnForMetUncertainty = cms.Sequence(process.patPFMETcorrMuonEnUp*
                                                                           process.patPFMETcorrMuonEnDown*
                                                                           process.patType1CorrectedPFSmearedMetMuonEnUp*
                                                                           process.patType1CorrectedPFSmearedMetMuonEnDown*
                                                                           process.patType1CorrectedPFMetMuonEnUp*
                                                                           process.patType1CorrectedPFMetMuonEnDown)
 else:
      process.patType1CorrectedPFMetMuonEnForMetUncertainty = cms.Sequence(process.patPFMETcorrMuonEnUp*
                                                                           process.patPFMETcorrMuonEnDown*
                                                                           process.patType1CorrectedPFMetMuonEnUp*
                                                                           process.patType1CorrectedPFMetMuonEnDown)


 ## electron energy resolution
 process.patPFMETcorrMuonRes = cms.EDProducer("ShiftedParticleMETcorrInputProducer",
                                               srcOriginal = patMuonCollection,
                                               srcShifted = cms.InputTag("smearedPatMuonsPFlow"))

 process.patPFMETcorrMuonResUp   = process.patPFMETcorrMuonRes.clone(srcShifted = cms.InputTag("smearedPatMuonsPFlowResUp"))
 process.patPFMETcorrMuonResDown = process.patPFMETcorrMuonRes.clone(srcShifted = cms.InputTag("smearedPatMuonsPFlowResDown"))

 process.patType1CorrectedPFMetMuonRes = cms.EDProducer("CorrectedPATMETProducer",
                                                         applyType2Corrections = cms.bool(False),
                                                         srcType1Corrections = cms.VInputTag(cms.InputTag("patPFMETcorrMuonRes")),
                                                         src = TypeIpatMetcollection,
                                                         applyType1Corrections = cms.bool(True))

 process.patType1CorrectedPFMetMuonResUp      = process.patType1CorrectedPFMetMuonRes.clone(srcType1Corrections = cms.VInputTag(cms.InputTag("patPFMETcorrMuonResUp")))

 process.patType1CorrectedPFMetMuonResDown    = process.patType1CorrectedPFMetMuonRes.clone(srcType1Corrections = cms.VInputTag(cms.InputTag("patPFMETcorrMuonResDown")))

 process.patType1CorrectedPFSmearedMetMuonRes = process.patType1CorrectedPFMetMuonRes.clone(src = cms.InputTag("patType1CorrectedPFSmearedMet"))

 process.patType1CorrectedPFSmearedMetMuonResUp   = process.patType1CorrectedPFSmearedMetMuonRes.clone(srcType1Corrections = cms.VInputTag(cms.InputTag("patPFMETcorrMuonResUp")))

 process.patType1CorrectedPFSmearedMetMuonResDown = process.patType1CorrectedPFSmearedMetMuonRes.clone(srcType1Corrections = cms.VInputTag(cms.InputTag("patPFMETcorrMuonResDown")))

 if doSmearing and not useData:
     process.patType1CorrectedPFMetMuonResForMetUncertainty = cms.Sequence(process.patPFMETcorrMuonRes*
                                                                           process.patType1CorrectedPFMetMuonRes*
                                                                           process.patType1CorrectedPFSmearedMetMuonRes)
     if MuonEnergyResUnc !=0:
         process.patType1CorrectedPFMetMuonResForMetUncertainty += process.patPFMETcorrMuonResUp*process.patPFMETcorrMuonResDown
         process.patType1CorrectedPFMetMuonResForMetUncertainty += process.patType1CorrectedPFMetMuonResUp*process.patType1CorrectedPFMetMuonResDown
         process.patType1CorrectedPFMetMuonResForMetUncertainty += process.patType1CorrectedPFSmearedMetMuonResUp*process.patType1CorrectedPFSmearedMetMuonResDown
         
                                                                           
 ############# taus contribution to the met
 process.patPFMETcorrTausEnUp = cms.EDProducer("ShiftedParticleMETcorrInputProducer",
                                                srcOriginal = patTauCollection,
                                                srcShifted = cms.InputTag("shiftedPatTausPFlowEnUp"))

 process.patPFMETcorrTausEnDown = process.patPFMETcorrTausEnUp.clone(srcShifted = cms.InputTag("shiftedPatTausPFlowEnDown"))

 process.patType1CorrectedPFSmearedMetTausEnUp = cms.EDProducer("CorrectedPATMETProducer",
                                                                 applyType2Corrections = cms.bool(False),
                                                                 srcType1Corrections = cms.VInputTag(cms.InputTag("patPFMETcorrTausEnUp")),
                                                                 src = cms.InputTag("patType1CorrectedPFSmearedMet"),
                                                                 applyType1Corrections = cms.bool(True))

 process.patType1CorrectedPFSmearedMetTausEnDown = process.patType1CorrectedPFSmearedMetTausEnUp.clone(srcType1Corrections = cms.VInputTag(cms.InputTag("patPFMETcorrTausEnDown")))
                                                                                          

 process.patType1CorrectedPFMetTausEnUp   = process.patType1CorrectedPFSmearedMetTausEnUp.clone(src = TypeIpatMetcollection)

 process.patType1CorrectedPFMetTausEnDown = process.patType1CorrectedPFSmearedMetTausEnDown.clone(src = TypeIpatMetcollection)

 if doSmearing and not useData :
      process.patType1CorrectedPFMetTausEnForMetUncertainty = cms.Sequence(process.patPFMETcorrTausEnUp*
                                                                           process.patPFMETcorrTausEnDown*
                                                                           process.patType1CorrectedPFSmearedMetTausEnUp*
                                                                           process.patType1CorrectedPFSmearedMetTausEnDown*
                                                                           process.patType1CorrectedPFMetTausEnUp*
                                                                           process.patType1CorrectedPFMetTausEnDown)
 else:
      process.patType1CorrectedPFMetTausEnForMetUncertainty = cms.Sequence(process.patPFMETcorrTausEnUp*
                                                                           process.patPFMETcorrTausEnDown*
                                                                           process.patType1CorrectedPFMetTausEnUp*
                                                                           process.patType1CorrectedPFMetTausEnDown)

 ## tau energy resolution
      
 process.patPFMETcorrTauRes = cms.EDProducer("ShiftedParticleMETcorrInputProducer",
                                              srcOriginal = patTauCollection,
                                              srcShifted = cms.InputTag("smearedPatTausPFlow"))

 
 process.patPFMETcorrTauResUp   = process.patPFMETcorrTauRes.clone(srcShifted = cms.InputTag("smearedPatTausPFlowResUp"))
 process.patPFMETcorrTauResDown = process.patPFMETcorrTauRes.clone(srcShifted = cms.InputTag("smearedPatTausPFlowResDown"))

     
 process.patType1CorrectedPFMetTauRes = cms.EDProducer("CorrectedPATMETProducer",
                                                        applyType2Corrections = cms.bool(False),
                                                        srcType1Corrections = cms.VInputTag(cms.InputTag("patPFMETcorrTauRes")),
                                                        src = TypeIpatMetcollection,
                                                        applyType1Corrections = cms.bool(True))

 process.patType1CorrectedPFMetTauResUp   = process.patType1CorrectedPFMetTauRes.clone(srcType1Corrections = cms.VInputTag(cms.InputTag("patPFMETcorrTauResUp")))

 process.patType1CorrectedPFMetTauResDown = process.patType1CorrectedPFMetTauRes.clone(srcType1Corrections = cms.VInputTag(cms.InputTag("patPFMETcorrTauResDown")))

 process.patType1CorrectedPFSmearedMetTauRes = process.patType1CorrectedPFMetTauRes.clone(src = cms.InputTag("patType1CorrectedPFSmearedMet"))

 process.patType1CorrectedPFSmearedMetTauResUp   = process.patType1CorrectedPFSmearedMetTauRes.clone(srcType1Corrections = cms.VInputTag(cms.InputTag("patPFMETcorrTauResUp")))

 process.patType1CorrectedPFSmearedMetTauResDown = process.patType1CorrectedPFSmearedMetTauRes.clone(srcType1Corrections = cms.VInputTag(cms.InputTag("patPFMETcorrTauResDown")))
     
 
 if doSmearing and not useData:
    
     process.patType1CorrectedPFMetTauResForMetUncertainty = cms.Sequence(process.patPFMETcorrTauRes*
                                                                          process.patType1CorrectedPFMetTauRes*
                                                                          process.patType1CorrectedPFSmearedMetTauRes)

     if TauEnergyResUnc != 0:
         
      process.patType1CorrectedPFMetTauResForMetUncertainty += process.patPFMETcorrTauResUp*process.patPFMETcorrTauResDown
      process.patType1CorrectedPFMetTauResForMetUncertainty += process.patType1CorrectedPFMetTauResUp*process.patType1CorrectedPFMetTauResDown
      process.patType1CorrectedPFMetTauResForMetUncertainty += process.patType1CorrectedPFSmearedMetTauResUp*process.patType1CorrectedPFSmearedMetTauResDown
      
      

 ############# photon contribution to the met
 process.patPFMETcorrPhotonEnUp = cms.EDProducer("ShiftedParticleMETcorrInputProducer",
                                                  srcOriginal = patPhotonCollection,
                                                  srcShifted = cms.InputTag("shiftedPatPhotonPFlowEnUp"))

 process.patPFMETcorrPhotonEnDown = process.patPFMETcorrPhotonEnUp.clone(srcShifted = cms.InputTag("shiftedPatPhotonPFlowEnDown"))

 process.patType1CorrectedPFSmearedMetPhotonEnUp = cms.EDProducer("CorrectedPATMETProducer",
                                                                 applyType2Corrections = cms.bool(False),
                                                                 srcType1Corrections = cms.VInputTag(cms.InputTag("patPFMETcorrPhotonEnUp")),
                                                                 src = cms.InputTag("patType1CorrectedPFSmearedMet"),
                                                                 applyType1Corrections = cms.bool(True))

 process.patType1CorrectedPFSmearedMetPhotonEnDown = process.patType1CorrectedPFSmearedMetPhotonEnUp.clone(srcType1Corrections = cms.VInputTag(cms.InputTag("patPFMETcorrPhotonEnDown")))
                                                                                          

 process.patType1CorrectedPFMetPhotonEnUp   = process.patType1CorrectedPFSmearedMetPhotonEnUp.clone(src = TypeIpatMetcollection)

 process.patType1CorrectedPFMetPhotonEnDown = process.patType1CorrectedPFSmearedMetPhotonEnDown.clone(src = TypeIpatMetcollection)

 if doSmearing and not useData :
      process.patType1CorrectedPFMetPhotonEnForMetUncertainty = cms.Sequence(process.patPFMETcorrPhotonEnUp*
                                                                             process.patPFMETcorrPhotonEnDown*
                                                                             process.patType1CorrectedPFSmearedMetPhotonEnUp*
                                                                             process.patType1CorrectedPFSmearedMetPhotonEnDown*
                                                                             process.patType1CorrectedPFMetPhotonEnUp*
                                                                             process.patType1CorrectedPFMetPhotonEnDown)
 else:
      process.patType1CorrectedPFMetPhotonEnForMetUncertainty = cms.Sequence(process.patPFMETcorrPhotonEnUp*
                                                                             process.patPFMETcorrPhotonEnDown*
                                                                             process.patType1CorrectedPFMetPhotonEnUp*
                                                                             process.patType1CorrectedPFMetPhotonEnDown)

 ## photon energy resolution      
 process.patPFMETcorrPhotonRes = cms.EDProducer("ShiftedParticleMETcorrInputProducer",
                                                srcOriginal = patPhotonCollection,
                                                srcShifted = cms.InputTag("smearedPatPhotonPFlow"))

 
 process.patPFMETcorrPhotonResUp   = process.patPFMETcorrPhotonRes.clone(srcShifted = cms.InputTag("smearedPatPhotonPFlowResUp"))
 process.patPFMETcorrPhotonResDown = process.patPFMETcorrPhotonRes.clone(srcShifted = cms.InputTag("smearedPatPhotonPFlowResDown"))

     
 process.patType1CorrectedPFMetPhotonRes = cms.EDProducer("CorrectedPATMETProducer",
                                                        applyType2Corrections = cms.bool(False),
                                                        srcType1Corrections = cms.VInputTag(cms.InputTag("patPFMETcorrPhotonRes")),
                                                        src = TypeIpatMetcollection,
                                                        applyType1Corrections = cms.bool(True))

 process.patType1CorrectedPFSmearedMetPhotonRes = process.patType1CorrectedPFMetPhotonRes.clone(src = cms.InputTag("patType1CorrectedPFSmearedMet"))

 process.patType1CorrectedPFMetPhotonResUp   = process.patType1CorrectedPFMetPhotonRes.clone(srcType1Corrections = cms.VInputTag(cms.InputTag("patPFMETcorrPhotonResUp")))
 process.patType1CorrectedPFMetPhotonResDown = process.patType1CorrectedPFMetPhotonRes.clone(srcType1Corrections = cms.VInputTag(cms.InputTag("patPFMETcorrPhotonResDown")))

 process.patType1CorrectedPFSmearedMetPhotonResUp   = process.patType1CorrectedPFSmearedMetPhotonRes.clone(srcType1Corrections = cms.VInputTag(cms.InputTag("patPFMETcorrPhotonResUp")))
 process.patType1CorrectedPFSmearedMetPhotonResDown = process.patType1CorrectedPFSmearedMetPhotonRes.clone(srcType1Corrections = cms.VInputTag(cms.InputTag("patPFMETcorrPhotonResDown")))
     
 
 if doSmearing and not useData:
    
     process.patType1CorrectedPFMetPhotonResForMetUncertainty = cms.Sequence(process.patPFMETcorrPhotonRes*
                                                                             process.patType1CorrectedPFMetPhotonRes*
                                                                             process.patType1CorrectedPFSmearedMetPhotonRes)

     if PhotonEnergyResUncEB != 0 or PhotonEnergyResUncEE != 0 :
      process.patType1CorrectedPFMetPhotonResForMetUncertainty += process.patPFMETcorrPhotonResUp*process.patPFMETcorrPhotonResDown
      process.patType1CorrectedPFMetPhotonResForMetUncertainty += process.patType1CorrectedPFMetPhotonResUp*process.patType1CorrectedPFMetPhotonResDown
      process.patType1CorrectedPFMetPhotonResForMetUncertainty += process.patType1CorrectedPFSmearedMetPhotonResUp*process.patType1CorrectedPFSmearedMetPhotonResDown
      

               
 #### Final Sequence        
 if doSmearing and not useData :
     process.metUncertainty = cms.Sequence(process.patJetSmearingSequence*
                                           process.shiftJetsForMetUncertainty*
                                           process.shiftElectronsForMetUncertainty*
                                           process.smearedElectronsForMetUncertainty*
                                           process.shiftMuonsForMetUncertainty*
                                           process.smearedMuonsForMetUncertainty*
                                           process.shiftTausForMetUncertainty*
                                           process.smearedTausForMetUncertainty*
                                           process.shiftPhotonForMetUncertainty*
                                           process.smearedPhotonForMetUncertainty*
                                           process.patMETforMetUncertainty*
                                           process.patPFSmearedMetForMetUncertainty*
                                           process.patType1CorrectedPFSmearedMetForMetUncertainty*
                                           process.patPFMETcorrJetForMetUncertainty*
                                           process.patType1CorrectedPFMetJetForMetUncertainty*
                                           process.patType1CorrectedPFMetJetResForMetUncertainty*
                                           process.patType1CorrectedPFMetUnclusteredForMetUncertainty*
                                           process.patType1CorrectedPFMetElectronEnForMetUncertainty*
                                           process.patPFMETcorrElectronResForMetUncertainty*
                                           process.patType1CorrectedPFMetElectronResForMetUncertainty*
                                           process.patType1CorrectedPFMetMuonEnForMetUncertainty*
                                           process.patType1CorrectedPFMetMuonResForMetUncertainty*
                                           process.patType1CorrectedPFMetTausEnForMetUncertainty*
                                           process.patType1CorrectedPFMetTauResForMetUncertainty*
                                           process.patType1CorrectedPFMetPhotonEnForMetUncertainty*
                                           process.patType1CorrectedPFMetPhotonResForMetUncertainty                                           
                                           )
     if not doTauSmearing :
        process.metUncertainty.remove(process.smearedTausForMetUncertainty)
        process.metUncertainty.remove(process.patType1CorrectedPFMetTauResForMetUncertainty)

     if not doPhotonSmearing :
        process.metUncertainty.remove(process.smearedPhotonForMetUncertainty)
        process.metUncertainty.remove(process.patType1CorrectedPFMetPhotonResForMetUncertainty)
         
        
 else : 
     process.metUncertainty = cms.Sequence( process.shiftJetsForMetUncertainty*
                                            process.shiftElectronsForMetUncertainty*
                                            process.shiftMuonsForMetUncertainty*
                                            process.shiftTausForMetUncertainty*
                                            process.shiftPhotonForMetUncertainty*
                                            process.patMETforMetUncertainty*
                                            process.patPFMETcorrJetForMetUncertainty*
                                            process.patType1CorrectedPFMetJetForMetUncertainty*
                                            process.patType1CorrectedPFMetUnclusteredForMetUncertainty*
                                            process.patType1CorrectedPFMetElectronEnForMetUncertainty*
                                            process.patType1CorrectedPFMetMuonEnForMetUncertainty*
                                            process.patType1CorrectedPFMetTausEnForMetUncertainty*
                                            process.patType1CorrectedPFMetPhotonEnForMetUncertainty
                                           )


 process.patseq += process.metUncertainty   

 ## Add in the output

 if doSmearing and not useData :
     process.out.outputCommands += [
                                    'keep *_smearedPatJetsPFlow*_*_*',
                                    'keep *_shiftedPatJetsPFlow*_*_*',
                                    'keep *_shiftedPatElectronsPFlowEn*_*_*',
                                    'keep *_shiftedPatMuonsPFlowEn*_*_*',
                                    'keep *_shiftedPatTausPFlowEn*_*_*',
                                    'keep *_shiftedPatPhotonPFlowEn*_*_*',
                                    'keep *_smearedPatMuonsPFlow*_*_*',
                                    'keep *_smearedPatElectronsPFlow*_*_*',
                                    'keep *_smearedPatTausPFlow*_*_*',
                                    'keep *_smearedPatPhotonPFlow*_*_*',                                    
                                    'keep *_patType1CorrectedPFSmearedMet*_*_*',
                                    'keep *_patType1CorrectedPF*MetJetEnUp*_*_*',
                                    'keep *_patType1CorrectedPF*MetJetEnDown*_*_*',
                                    'keep *_patType1CorrectedPF*MetJetRes*_*_*',
                                    'keep *_patType1CorrectedPF*MetUnclusteredEn*_*_*',
                                    'keep *_patType1CorrectedPF*MetElectronEn*_*_*',
                                    'keep *_patType1CorrectedPF*MetMuonEn*_*_*',
                                    'keep *_patType1CorrectedPF*MetTausEn*_*_*',
                                    'keep *_patType1CorrectedPF*MetPhotonEn*_*_*',
                                    'keep *_patType1CorrectedPF*MetElectronRes*_*_*',
                                    'keep *_patType1CorrectedPF*MetMuonRes*_*_*',
                                    'keep *_patType1CorrectedPF*MetTauRes*_*_*',
                                    'keep *_patType1CorrectedPF*MetPhotonRes*_*_*',
                                    ] 
 else :
     process.out.outputCommands += ['keep *_shiftedPatJetsPFlow*_*_*',
                                    'keep *_shiftedPatElectronsPFlowEn*_*_*',
                                    'keep *_shiftedPatMuonsPFlowEn*_*_*',
                                    'keep *_shiftedPatTausPFlowEn*_*_*',
                                    'keep *_shiftedPatPhotonPFlowEn*_*_*',
                                    'keep *_patType1CorrectedPF*MetJetEnUp*_*_*',
                                    'keep *_patType1CorrectedPF*MetJetEnDown*_*_*',
                                    'keep *_patType1CorrectedPF*MetUnclusteredEn*_*_*',
                                    'keep *_patType1CorrectedPF*MetElectronEn*_*_*',
                                    'keep *_patType1CorrectedPF*MetMuonEn*_*_*',
                                    'keep *_patType1CorrectedPF*MetTausEn*_*_*',
                                    'keep *_patType1CorrectedPF*MetPhotonEn*_*_*',
                                    ]
