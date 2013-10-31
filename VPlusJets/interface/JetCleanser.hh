// JetCleanser Package
// Questions/Comments? dkrohn@physics.harvard.edu mattlow@uchicago.edu schwartz@physics.harvard.edu liantaow@uchicago.edu
//
// Copyright (c) 2013
// David Krohn, Matthew Low, Matthew Schwartz, and Lian-Tao Wang
//
// $Id$
//
//----------------------------------------------------------------------
// This file is part of FastJet contrib.
//
// It is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2 of the License, or (at
// your option) any later version.
//
// It is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
// License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this code. If not, see <http://www.gnu.org/licenses/>.
//----------------------------------------------------------------------

#ifndef __FASTJET_CONTRIB_JETCLEANSER_HH__
#define __FASTJET_CONTRIB_JETCLEANSER_HH__

#include <fastjet/internal/base.hh>

#include "fastjet/ClusterSequence.hh"
#include "fastjet/Error.hh"
#include "fastjet/JetDefinition.hh"
//#include "fastjet/FunctionOfPseudoJet.hh"

#include <map>
#include <sstream>
#include <string>

//FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

//namespace contrib{

//------------------------------------------------------------------------
/// \class JetCleanser
/// This class implements jet cleansing which is a substructure technique
/// designed to correct for the effects of pileup.
///
/// Jet cleansing takes as input either:
///  (A) "calorimeter" or "particle" four-vectors which are charged and neutral 
///       particles all grouped together, charged tracks coming from the primary 
///       interaction, and charged tracks coming from pileup.
///  (B) neutral particles, charged tracks coming from the primary
///       interaction, and charged tracks coming from pileup.
///
/// The output is a jet for which the constituents themselves have been 
/// corrected so that the output can be used to calculate both 
/// kinematic quantities and shape variables.
///
class JetCleanser {

public:
  enum cleansing_mode {
    jvf_cleansing,      /// gamma0 = gamma1
    linear_cleansing,   /// approximate gamma0 as constant
    gaussian_cleansing  /// approximate gamma's as truncated gaussians and maximize likelihood
  };

  enum input_mode {
    input_nc_together,  /// input is all charged and neutrals together (as particles or calo cells),
                        ///          charged primary interaction particles, and charged pileup particles
    input_nc_separate   /// input is all neutral particles, charged primary interaction particles,
                        ///          and charged pileup particles
  };

public:
  // constructors
  JetCleanser(fastjet::JetDefinition subjet_def, cleansing_mode cmode, input_mode imode);
  JetCleanser(double rsub, cleansing_mode cmode, input_mode imode);

  // destructor
  ~JetCleanser(){}

  // settings
  void SetGroomingParameters(double fcut, int nsjmin);
  inline void SetTrimming(double fcut) { SetGroomingParameters(fcut,0); }
  inline void SetFiltering(int nsj) { SetGroomingParameters(1.0,nsj); }
  void SetLinearParameters(double g0_mean=0.67);
  void SetGaussianParameters(double g0_mean=0.67, double g1_mean=0.67, double g0_width=0.15, double g1_width=0.25);

  // standard usage
  std::string description() const;
  fastjet::PseudoJet operator()(const fastjet::PseudoJet & jet,                             // For use with input_nc_together mode,
                       const std::vector<fastjet::PseudoJet> & tracks_lv,          // it takes a plain jet, charged tracks from
                       const std::vector<fastjet::PseudoJet> & tracks_pu) const;   // the leading vertex, and charged tracks
                                                                          // from pileup.
  fastjet::PseudoJet operator()(const std::vector<fastjet::PseudoJet> & neutrals_all,       // For use with input_nc_separate mode,
                       const std::vector<fastjet::PseudoJet> & tracks_lv,          // it takes all neutral particles, charged 
                       const std::vector<fastjet::PseudoJet> & tracks_pu) const;   // tracks from the leading vertex, and charged
                                                                          // tracks from pileup.
  void   _RunTests();


private:
  const static double jc_zero;// = 1.0e-6;

  double _rsub;
  double _fcut;
  double _nsjmin;
  fastjet::JetDefinition _subjet_def;

  cleansing_mode _cleansing_mode;
  input_mode _input_mode;

  double _linear_gamma0_mean;
  double _gaussian_gamma0_mean;
  double _gaussian_gamma0_width;
  double _gaussian_gamma1_mean;
  double _gaussian_gamma1_width;

  // helper functions
  void   _SetDefaults();
  void   _CheckRescalingValues(double & pt_all, const double & ptc_lv, const double & ptc_pu) const;
  double _GetSubjetRescaling_nctogether(double pt_all, double ptc_lv, double ptc_pu) const;
  double _GetSubjetRescaling_ncseparate(double ptn_all, double ptc_lv, double ptc_pu) const;
  double _GaussianGetMinimizedGamma0(double pt_all, double ptc_lv, double ptc_pu) const;
  double _GaussianGetGamma1(double gamma0, double pt_all, double ptc_lv, double ptc_pu) const;
  double _GaussianFunction(double x, void * params) const;
  void   _RunTestRescaling(double pt_all, double ptc_lv, double ptc_pu) const;
};
//const double JetCleanser::jc_zero = 1.0e-6;


// helper function
std::vector<fastjet::PseudoJet> RescalePseudoJetVector(const std::vector<fastjet::PseudoJet> & jets, const double s_factor);

// helper function
fastjet::PseudoJet RescalePseudoJetConstituents(const fastjet::PseudoJet & jet, const double s_factor);

// helper function
std::vector< std::vector<fastjet::PseudoJet> > ClusterSets(const fastjet::JetDefinition & jet_def, 
                                                  const std::vector<fastjet::PseudoJet> & cluster_set,
                                                  const std::vector< std::vector<fastjet::PseudoJet> > & follow_sets,
                                                  const double &ptmin=0.0);
        
        

//------------------------------------------------------------------------
/// \class FollowSetGhostInfo
/// This class keeps tracks of which set of PseudoJets and the index where
/// the associated ghost fastjet::PseudoJet came from.
/// 
class FollowSetGhostInfo : public fastjet::PseudoJet::UserInfoBase {
  public:
    inline FollowSetGhostInfo(int set_id, int ind_id) { _set_id = set_id; _ind_id = ind_id; }
    inline int GetSetId() { return _set_id; }
    inline int GetIndId() { return _ind_id; }

  private:
    int _set_id;
    int _ind_id;
};

//} // namespace contrib

//FASTJET_END_NAMESPACE

#endif  // __FASTJET_CONTRIB_JETCLEANSER_HH__
