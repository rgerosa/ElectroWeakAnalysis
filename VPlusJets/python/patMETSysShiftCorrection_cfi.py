import FWCore.ParameterSet.Config as cms

def metShiftSystematicCorrection(process,
                                 isMC,
                                 useSmearedCollection,
                                 runCorrection,
                                 vertexCollection,
                                 didPhotonSmearing,
                                 didTauSmearing):



 print "                                      "
 print "######################################"
 print "Met Shift Systematic Correction Module"
 print "######################################"
 print "                                      "

 print "Chosen Options:                      "
 print "is running on MC                   = %d"%isMC
 print "useSmeardCollection                = %d"%useSmearedCollection
 print "RunCorrection to be applied        = %s"%runCorrection
 print "Vertex collection to be considered = %s"%vertexCollection
 print "Photon smearing has been done      = %d"%didPhotonSmearing
 print "Tau smearing has been done         = %d"%didTauSmearing

 process.load("JetMETCorrections.Type1MET.pfMETsysShiftCorrections_cfi")

 ######## Some basic parameters:
 patTypeIMetCorrected = 'patMETsPFlow'

 patTypeIMetCorrectedNoSmeared = []
 patTypeIMetCorrectedNoSmeared.append('patType1CorrectedPFMetElectronEnDown')
 patTypeIMetCorrectedNoSmeared.append('patType1CorrectedPFMetElectronEnUp')
 patTypeIMetCorrectedNoSmeared.append('patType1CorrectedPFMetElectronRes')
 patTypeIMetCorrectedNoSmeared.append('patType1CorrectedPFMetElectronResDown')
 patTypeIMetCorrectedNoSmeared.append('patType1CorrectedPFMetElectronResUp')
 patTypeIMetCorrectedNoSmeared.append('patType1CorrectedPFMetJetEnDown')
 patTypeIMetCorrectedNoSmeared.append('patType1CorrectedPFMetJetEnUp')
 patTypeIMetCorrectedNoSmeared.append('patType1CorrectedPFMetMuonEnDown')
 patTypeIMetCorrectedNoSmeared.append('patType1CorrectedPFMetMuonEnUp')
 patTypeIMetCorrectedNoSmeared.append('patType1CorrectedPFMetMuonRes')
 patTypeIMetCorrectedNoSmeared.append('patType1CorrectedPFMetMuonResDown') 
 patTypeIMetCorrectedNoSmeared.append('patType1CorrectedPFMetMuonResUp')
 patTypeIMetCorrectedNoSmeared.append('patType1CorrectedPFMetPhotonEnDown')
 patTypeIMetCorrectedNoSmeared.append('patType1CorrectedPFMetPhotonEnUp')
 patTypeIMetCorrectedNoSmeared.append('patType1CorrectedPFMetPhotonRes')
 patTypeIMetCorrectedNoSmeared.append('patType1CorrectedPFMetPhotonResDown')
 patTypeIMetCorrectedNoSmeared.append('patType1CorrectedPFMetPhotonResUp')
 patTypeIMetCorrectedNoSmeared.append('patType1CorrectedPFMetTauRes')
 patTypeIMetCorrectedNoSmeared.append('patType1CorrectedPFMetTauResDown')
 patTypeIMetCorrectedNoSmeared.append('patType1CorrectedPFMetTauResUp')
 patTypeIMetCorrectedNoSmeared.append('patType1CorrectedPFMetTausEnDown')
 patTypeIMetCorrectedNoSmeared.append('patType1CorrectedPFMetTausEnUp')
 patTypeIMetCorrectedNoSmeared.append('patType1CorrectedPFMetUnclusteredEnDown')
 patTypeIMetCorrectedNoSmeared.append('patType1CorrectedPFMetUnclusteredEnUp')

 patTypeIMetCorrectedSmeared = []
 patTypeIMetCorrectedSmeared.append('patType1CorrectedPFSmearedMet')
 patTypeIMetCorrectedSmeared.append('patType1CorrectedPFSmearedMetElectronEnDown')
 patTypeIMetCorrectedSmeared.append('patType1CorrectedPFSmearedMetElectronEnUp')
 patTypeIMetCorrectedSmeared.append('patType1CorrectedPFSmearedMetElectronRes')
 patTypeIMetCorrectedSmeared.append('patType1CorrectedPFSmearedMetElectronResDown')
 patTypeIMetCorrectedSmeared.append('patType1CorrectedPFSmearedMetElectronResUp')
 patTypeIMetCorrectedSmeared.append('patType1CorrectedPFSmearedMetJetEnDown')
 patTypeIMetCorrectedSmeared.append('patType1CorrectedPFSmearedMetJetEnUp')
 patTypeIMetCorrectedSmeared.append('patType1CorrectedPFSmearedMetJetResDown')
 patTypeIMetCorrectedSmeared.append('patType1CorrectedPFSmearedMetJetResUp')
 patTypeIMetCorrectedSmeared.append('patType1CorrectedPFSmearedMetMuonEnDown')
 patTypeIMetCorrectedSmeared.append('patType1CorrectedPFSmearedMetMuonEnUp')
 patTypeIMetCorrectedSmeared.append('patType1CorrectedPFSmearedMetMuonRes')
 patTypeIMetCorrectedSmeared.append('patType1CorrectedPFSmearedMetMuonResDown')
 patTypeIMetCorrectedSmeared.append('patType1CorrectedPFSmearedMetMuonResUp')
 patTypeIMetCorrectedSmeared.append('patType1CorrectedPFSmearedMetPhotonEnDown')
 patTypeIMetCorrectedSmeared.append('patType1CorrectedPFSmearedMetPhotonEnUp')
 patTypeIMetCorrectedSmeared.append('patType1CorrectedPFSmearedMetPhotonRes')
 patTypeIMetCorrectedSmeared.append('patType1CorrectedPFSmearedMetPhotonResDown')
 patTypeIMetCorrectedSmeared.append('patType1CorrectedPFSmearedMetPhotonResUp')
 patTypeIMetCorrectedSmeared.append('patType1CorrectedPFSmearedMetTauRes')
 patTypeIMetCorrectedSmeared.append('patType1CorrectedPFSmearedMetTauResDown')
 patTypeIMetCorrectedSmeared.append('patType1CorrectedPFSmearedMetTauResUp')
 patTypeIMetCorrectedSmeared.append('patType1CorrectedPFSmearedMetTausEnDown')
 patTypeIMetCorrectedSmeared.append('patType1CorrectedPFSmearedMetTausEnUp')                                
 patTypeIMetCorrectedSmeared.append('patType1CorrectedPFSmearedMetUnclusteredEnDown')                                
 patTypeIMetCorrectedSmeared.append(B'patType1CorrectedPFSmearedMetUnclusteredEnUp')

 ## basic declaration of the sys shift producer
 process.pfMEtSysShiftCorr = cms.EDProducer("SysShiftMETcorrInputProducer",
                                             src = cms.InputTag('pfMet'),
                                             srcVertices = vertexCollection,
                                             parameter   = process.pfMEtSysShiftCorrParameters_2012runAvsNvtx_data )

 process.pfMetShiftCorrected = cms.EDProducer("CorrectedPATMETProducer",
                                               src = cms.InputTag('pfMet'),
                                               applyType1Corrections = cms.bool(True),
                                               srcType1Corrections = cms.VInputTag(cms.InputTag('pfMEtSysShiftCorr')),
                                               applyType2Corrections = cms.bool(False))
 
 ### taking the runCorrection as reference, apply the following corrections, most of them are taken from JetMETCorrections.Type1MET.pfMETsysShiftCorrections_cfi
 if (runCorrection == "A" or runCorrection == "Run2012A" or runCorrection == "2012A") and not isMC :
    process.pfMEtSysShiftCorr.parameter =  process.pfMEtSysShiftCorrParameters_2012runAvsNvtx_data

 elif (runCorrection == "A" or runCorrection == "Run2012A" or runCorrection == "2012A") and isMC :
    process.pfMEtSysShiftCorr.parameter =  process.pfMEtSysShiftCorrParameters_2012runAvsNvtx_mc

 elif (runCorrection == "B" or runCorrection == "Run2012B" or runCorrection == "2012B") and not isMC :
    process.pfMEtSysShiftCorr.parameter =  process.pfMEtSysShiftCorrParameters_2012runBvsNvtx_data

 elif (runCorrection == "B" or runCorrection == "Run2012B" or runCorrection == "2012B") and  isMC :
    process.pfMEtSysShiftCorr.parameter =  process.pfMEtSysShiftCorrParameters_2012runBvsNvtx_mc

 elif (runCorrection == "ABCD" or runCorrection == "Run2012ABCD" or runCorrection == "2012ABCD") and not isMC :
     process.pfMEtSysShiftCorr.parameter =  cms.PSet( # CV: ReReco data + Summer'13 JEC
                                                      px = cms.string("+4.83642e-02 + 2.48870e-01*Nvtx"),
                                                      py = cms.string("-1.50135e-01 - 8.27917e-02*Nvtx"))

 elif (runCorrection == "ABCD" or runCorrection == "Run2012ABCD" or runCorrection == "2012ABCD") and isMC :
     process.pfMEtSysShiftCorr.parameter =      cms.PSet( # CV: Summer'12 MC + Summer'13 JEC
                                                          px = cms.string("+1.62861e-01 - 2.38517e-02*Nvtx"),
                                                          py = cms.string("+3.60860e-01 - 1.30335e-01*Nvtx"))


 ### input met to be corrected:

 process.metShiftSystematicCorrectionSequence = cms.Sequence()

 if not isMC :

  nameShifter = patTypeIMetCorrected+'MEtSysShiftCorr'
  nameCorrMet = patTypeIMetCorrected+'SysShifted'
  
  setattr(process, nameShifter, process.pfMEtSysShiftCorr.clone(src = cms.InputTag(patTypeIMetCorrected)))
  setattr(process, nameCorrMet, process.pfMetShiftCorrected.clone(src = cms.InputTag(patTypeIMetCorrected),
                                                                  srcType1Corrections = cms.VInputTag(cms.InputTag(nameShifter))))

  process.metShiftSystematicCorrectionSequence += getattr(process,nameShifter)* getattr(process,nameCorrMet)
  
 else :

    nameShifter = patTypeIMetCorrected+'MEtSysShiftCorr'
    nameCorrMet = patTypeIMetCorrected+'SysShifted'

    setattr(process, nameShifter, process.pfMEtSysShiftCorr.clone(src = cms.InputTag(patTypeIMetCorrected)))
    setattr(process, nameCorrMet, process.pfMetShiftCorrected.clone(src = cms.InputTag(patTypeIMetCorrected),
                                                                    srcType1Corrections = cms.VInputTag(cms.InputTag(nameShifter))))


    process.metShiftSystematicCorrectionSequence += getattr(process,nameShifter)* getattr(process,nameCorrMet)

    if not useSmearedCollection:

        for module in patTypeIMetCorrectedNoSmeared :

            if not didPhotonSmearing and "Photon" in module : continue
            if not didTauSmearing and "Tau" in module : continue
            
            nameShifter = module+'MEtSysShiftCorr'
            nameCorrMet = module+'SysShifted'

            setattr(process,name, process.pfMEtSysShiftCorr.clone(src = cms.InputTag(module)))
            setattr(process,nameCorrMet, process.pfMetShiftCorrected.clone(src = cms.InputTag(module),
                                                                           srcType1Corrections = cms.VInputTag(cms.InputTag(nameShifter))))

            process.metShiftSystematicCorrectionSequence += getattr(process,nameShifter)*getattr(process,nameCorrMet)

    else : ## apply to both smeared and not smeared collection

        for module in patTypeIMetCorrectedNoSmeared :

            if not didPhotonSmearing and "Photon" in module : continue
            if not didTauSmearing and "Tau" in module : continue

            nameShifter = module+'MEtSysShiftCorr'
            nameCorrMet = module+'SysShifted'

            setattr(process,nameShifter, process.pfMEtSysShiftCorr.clone(src = cms.InputTag(module)))
            setattr(process,nameCorrMet, process.pfMetShiftCorrected.clone(src = cms.InputTag(module),
                                                                           srcType1Corrections = cms.VInputTag(cms.InputTag(nameShifter))))

            process.metShiftSystematicCorrectionSequence += getattr(process,nameShifter)* getattr(process,nameCorrMet)


        for module in patTypeIMetCorrectedSmeared :

            if not didPhotonSmearing and "Photon" in module : continue
            if not didTauSmearing and "Tau" in module : continue

            nameShifter = module+'MEtSysShiftCorr'
            nameCorrMet = module+'SysShifted'

            setattr(process,nameShifter, process.pfMEtSysShiftCorr.clone(src = cms.InputTag(module)))
            setattr(process,nameCorrMet, process.pfMetShiftCorrected.clone(src = cms.InputTag(module),
                                                                           srcType1Corrections = cms.VInputTag(cms.InputTag(nameShifter))))

            process.metShiftSystematicCorrectionSequence += getattr(process,nameShifter)*getattr(process,nameCorrMet)
