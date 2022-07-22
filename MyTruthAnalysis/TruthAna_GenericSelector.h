#ifndef MYTRUTHANALYSIS_TRUTHANA_GENERICSELECTOR_H
#define MYTRUTHANALYSIS_TRUTHANA_GENERICSELECTOR_H

#include <vector>
#include <TLorentzVector.h>

#include "xAODTruth/TruthParticleContainer.h"
#include "xAODJet/JetContainer.h"

class TruthAna_GenericSelector
{
 public: 
  TruthAna_GenericSelector();
  ~TruthAna_GenericSelector();

  //super cool c++11 feature
  enum class CUT_PDG {YES,NO,ABS};
  enum class CUT_STATUS {YES,NO};

  //Selector methods

  std::vector<TLorentzVector> GenericParticle(const xAOD::TruthParticleContainer *cont, float pt, float eta, CUT_PDG cut_pdg, int pdgid, CUT_STATUS cut_status, int status);

  std::vector<TLorentzVector> GenericJet(const xAOD::JetContainer *cont, float pt, float eta);
  
  std::vector< std::pair< TLorentzVector,int > > GetAllPartons(const xAOD::TruthParticleContainer *cont, CUT_STATUS cut_status, int status);
  
  std::vector<TLorentzVector> RetrieveParton(std::vector< std::pair< TLorentzVector,int > > generic_stable_parton, CUT_PDG cut_pdg, int pdgid);

  float CalculateHT(std::vector<TLorentzVector> cont_jet, 
		    std::vector<TLorentzVector> cont_electron, std::vector<TLorentzVector> cont_muon, std::vector<TLorentzVector> cont_tau, 
		    std::vector<TLorentzVector> cont_neutrino);

};


#endif 
