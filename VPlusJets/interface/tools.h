/*
 * =====================================================================================
 * 
 *       Filename:  tools.h
 * 
 *    Description:  
 * 
 *        Version:  1.0
 *        Created:  08/07/13 02:04:15 CDT
 *       Revision:  none
 *       Compiler:  gcc, root
 * 
 *         Author:  Zijun Xu, xuzijun123@gmail.com
 *        Company:  School of Physics, Peking Univ.
 * 
 * =====================================================================================
 */

#ifndef  TOOLS_INC
#define  TOOLS_INC


// system include files
#include <memory>
#include <string>
#include <iostream>
#include <map>
#include <fstream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h" 

#include "TFile.h"
#include "TTree.h"

//#include "ElectroWeakAnalysis/VPlusJets/interface/JetTreeFiller.h"
//#include "ElectroWeakAnalysis/VPlusJets/interface/GroomedJetFiller.h"
//#include "ElectroWeakAnalysis/VPlusJets/interface/PhotonTreeFiller.h"
//#include "ElectroWeakAnalysis/VPlusJets/interface/VtoElectronTreeFiller.h"
//#include "ElectroWeakAnalysis/VPlusJets/interface/VtoMuonTreeFiller.h"
//#include "ElectroWeakAnalysis/VPlusJets/interface/MCTreeFiller.h"
#include "ElectroWeakAnalysis/VPlusJets/interface/tools.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include <fastjet/ClusterSequence.hh>
#include <fastjet/GhostedAreaSpec.hh>
#include <fastjet/ClusterSequenceArea.hh>
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"
#include "fastjet/tools/MassDropTagger.hh"
#include "fastjet/GhostedAreaSpec.hh"

#include "JetSubstructure/SubstructureTools/interface/PseudoJetUserInfo.h"
#include "JetSubstructure/SubstructureTools/interface/JetSubstructureTools.h"



void print_p4(fastjet::PseudoJet tmpJ, std::string tmpName="",bool extra_info=0);
void BREAK(std::string info="");

#endif   /* ----- #ifndef TOOLS_INC  ----- */

