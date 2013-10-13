import FWCore.ParameterSet.Config as cms

def metShiftSystematicCorrection(process,
                                 isMC,
                                 useSmearedCollection,
                                 runCorrection,
                                 vertexCollection,
                                 didPhotonSystematic,
                                 didTauSystematic,
                                 patTypeIMetCorrected,
                                 patTypeIMetCorrectedForMetUncertainty,
                                 patTypeIMetCorrectedShifted,
                                 patTypeIMetCorrectedShiftedForMetUncertainty):



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
 print "Photon systematic has been done    = %d"%didPhotonSystematic
 print "Tau systematic has been done       = %d"%didTauSystematic

 process.load("JetMETCorrections.Type1MET.pfMETsysShiftCorrections_cfi")

 ######## Some basic parameters:

 if not useSmearedCollection:
     
  patTypeIMetCorrectedForMetUncertainty.append('patType1CorrectedPFMetElectronEnDown')
  patTypeIMetCorrectedForMetUncertainty.append('patType1CorrectedPFMetElectronEnUp')
  patTypeIMetCorrectedForMetUncertainty.append('patType1CorrectedPFMetElectronRes')
  patTypeIMetCorrectedForMetUncertainty.append('patType1CorrectedPFMetElectronResDown')
  patTypeIMetCorrectedForMetUncertainty.append('patType1CorrectedPFMetElectronResUp')
  patTypeIMetCorrectedForMetUncertainty.append('patType1CorrectedPFMetJetEnDown')
  patTypeIMetCorrectedForMetUncertainty.append('patType1CorrectedPFMetJetEnUp')
  patTypeIMetCorrectedForMetUncertainty.append('patType1CorrectedPFMetMuonEnDown')
  patTypeIMetCorrectedForMetUncertainty.append('patType1CorrectedPFMetMuonEnUp')
  patTypeIMetCorrectedForMetUncertainty.append('patType1CorrectedPFMetMuonRes')
  patTypeIMetCorrectedForMetUncertainty.append('patType1CorrectedPFMetMuonResDown') 
  patTypeIMetCorrectedForMetUncertainty.append('patType1CorrectedPFMetMuonResUp')
  patTypeIMetCorrectedForMetUncertainty.append('patType1CorrectedPFMetPhotonEnDown')
  patTypeIMetCorrectedForMetUncertainty.append('patType1CorrectedPFMetPhotonEnUp')
  patTypeIMetCorrectedForMetUncertainty.append('patType1CorrectedPFMetPhotonRes')
  patTypeIMetCorrectedForMetUncertainty.append('patType1CorrectedPFMetPhotonResDown')
  patTypeIMetCorrectedForMetUncertainty.append('patType1CorrectedPFMetPhotonResUp')
  patTypeIMetCorrectedForMetUncertainty.append('patType1CorrectedPFMetTauRes')
  patTypeIMetCorrectedForMetUncertainty.append('patType1CorrectedPFMetTauResDown')
  patTypeIMetCorrectedForMetUncertainty.append('patType1CorrectedPFMetTauResUp')
  patTypeIMetCorrectedForMetUncertainty.append('patType1CorrectedPFMetTausEnDown')
  patTypeIMetCorrectedForMetUncertainty.append('patType1CorrectedPFMetTausEnUp')
  patTypeIMetCorrectedForMetUncertainty.append('patType1CorrectedPFMetUnclusteredEnDown')
  patTypeIMetCorrectedForMetUncertainty.append('patType1CorrectedPFMetUnclusteredEnUp')

 else :

  patTypeIMetCorrectedForMetUncertainty.append('patType1CorrectedPFSmearedMet')
  patTypeIMetCorrectedForMetUncertainty.append('patType1CorrectedPFSmearedMetElectronEnDown')
  patTypeIMetCorrectedForMetUncertainty.append('patType1CorrectedPFSmearedMetElectronEnUp')
  patTypeIMetCorrectedForMetUncertainty.append('patType1CorrectedPFSmearedMetElectronRes')
  patTypeIMetCorrectedForMetUncertainty.append('patType1CorrectedPFSmearedMetElectronResDown')
  patTypeIMetCorrectedForMetUncertainty.append('patType1CorrectedPFSmearedMetElectronResUp')
  patTypeIMetCorrectedForMetUncertainty.append('patType1CorrectedPFSmearedMetJetEnDown')
  patTypeIMetCorrectedForMetUncertainty.append('patType1CorrectedPFSmearedMetJetEnUp')
  patTypeIMetCorrectedForMetUncertainty.append('patType1CorrectedPFSmearedMetJetResDown')
  patTypeIMetCorrectedForMetUncertainty.append('patType1CorrectedPFSmearedMetJetResUp')
  patTypeIMetCorrectedForMetUncertainty.append('patType1CorrectedPFSmearedMetMuonEnDown')
  patTypeIMetCorrectedForMetUncertainty.append('patType1CorrectedPFSmearedMetMuonEnUp')
  patTypeIMetCorrectedForMetUncertainty.append('patType1CorrectedPFSmearedMetMuonRes')
  patTypeIMetCorrectedForMetUncertainty.append('patType1CorrectedPFSmearedMetMuonResDown')
  patTypeIMetCorrectedForMetUncertainty.append('patType1CorrectedPFSmearedMetMuonResUp')
  patTypeIMetCorrectedForMetUncertainty.append('patType1CorrectedPFSmearedMetPhotonEnDown')
  patTypeIMetCorrectedForMetUncertainty.append('patType1CorrectedPFSmearedMetPhotonEnUp')
  patTypeIMetCorrectedForMetUncertainty.append('patType1CorrectedPFSmearedMetPhotonRes')
  patTypeIMetCorrectedForMetUncertainty.append('patType1CorrectedPFSmearedMetPhotonResDown')
  patTypeIMetCorrectedForMetUncertainty.append('patType1CorrectedPFSmearedMetPhotonResUp')
  patTypeIMetCorrectedForMetUncertainty.append('patType1CorrectedPFSmearedMetTauRes')
  patTypeIMetCorrectedForMetUncertainty.append('patType1CorrectedPFSmearedMetTauResDown')
  patTypeIMetCorrectedForMetUncertainty.append('patType1CorrectedPFSmearedMetTauResUp')
  patTypeIMetCorrectedForMetUncertainty.append('patType1CorrectedPFSmearedMetTausEnDown')
  patTypeIMetCorrectedForMetUncertainty.append('patType1CorrectedPFSmearedMetTausEnUp')                                
  patTypeIMetCorrectedForMetUncertainty.append('patType1CorrectedPFSmearedMetUnclusteredEnDown')                                
  patTypeIMetCorrectedForMetUncertainty.append('patType1CorrectedPFSmearedMetUnclusteredEnUp')

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

  nameShifter = patTypeIMetCorrected[0]+'MEtSysShiftCorr'
  nameCorrMet = patTypeIMetCorrected[0]+'SysShifted'
  
  setattr(process, nameShifter, process.pfMEtSysShiftCorr.clone(src = cms.InputTag(patTypeIMetCorrected[0])))
  setattr(process, nameCorrMet, process.pfMetShiftCorrected.clone(src = cms.InputTag(patTypeIMetCorrected[0]),
                                                                  srcType1Corrections = cms.VInputTag(cms.InputTag(nameShifter))))

  process.metShiftSystematicCorrectionSequence += getattr(process,nameShifter)*getattr(process,nameCorrMet)

  patTypeIMetCorrectedShifted.append(nameCorrMet)  
  
 else :

    nameShifter = patTypeIMetCorrected[0]+'MEtSysShiftCorr'
    nameCorrMet = patTypeIMetCorrected[0]+'SysShifted'

    setattr(process, nameShifter, process.pfMEtSysShiftCorr.clone(src = cms.InputTag(patTypeIMetCorrected[0])))
    setattr(process, nameCorrMet, process.pfMetShiftCorrected.clone(src = cms.InputTag(patTypeIMetCorrected[0]),
                                                                    srcType1Corrections = cms.VInputTag(cms.InputTag(nameShifter))))


    process.metShiftSystematicCorrectionSequence += getattr(process,nameShifter)* getattr(process,nameCorrMet)

    patTypeIMetCorrectedShifted.append(nameCorrMet)  

    if not useSmearedCollection:

        for module in patTypeIMetCorrectedForMetUncertainty :

            if not didPhotonSystematic and "Photon" in module : continue
            if not didTauSystematic    and "Tau"    in module : continue
            
            nameShifter = module+'MEtSysShiftCorr'
            nameCorrMet = module+'SysShifted'

            setattr(process,nameShifter, process.pfMEtSysShiftCorr.clone(src = cms.InputTag(module)))
            setattr(process,nameCorrMet, process.pfMetShiftCorrected.clone(src = cms.InputTag(module),
                                                                           srcType1Corrections = cms.VInputTag(cms.InputTag(nameShifter))))

            process.metShiftSystematicCorrectionSequence += getattr(process,nameShifter)*getattr(process,nameCorrMet)

            patTypeIMetCorrectedShiftedForMetUncertainty.append(nameCorrMet)

    else : ## apply to both smeared and not smeared collection

        for module in patTypeIMetCorrectedForMetUncertainty :

            if not didPhotonSystematic and "Photon" in module : continue
            if not didTauSystematic    and "Tau"    in module : continue

            nameShifter = module+'MEtSysShiftCorr'
            nameCorrMet = module+'SysShifted'

            setattr(process,nameShifter, process.pfMEtSysShiftCorr.clone(src = cms.InputTag(module)))
            setattr(process,nameCorrMet, process.pfMetShiftCorrected.clone(src = cms.InputTag(module),
                                                                           srcType1Corrections = cms.VInputTag(cms.InputTag(nameShifter))))

            process.metShiftSystematicCorrectionSequence += getattr(process,nameShifter)*getattr(process,nameCorrMet)

            patTypeIMetCorrectedShiftedForMetUncertainty.append(nameCorrMet)


 if useSmearedCollection: ## the shifted smeared met should be the basic one, the not smear can be used to assign a uncertainty in case

  NameTemp = patTypeIMetCorrectedShifted[0]
  patTypeIMetCorrectedShifted[0] = patTypeIMetCorrectedShiftedForMetUncertainty[0]
  patTypeIMetCorrectedShiftedForMetUncertainty[0] = NameTemp


 print "                         "
 print "Corrected MET to be used as reference %s       "%patTypeIMetCorrectedShifted
 print "Corrected MET to be used for MET Uncertainty %s"%patTypeIMetCorrectedShiftedForMetUncertainty
 print "                         "
