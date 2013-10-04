// $Id: FunComptonTheory.h 236 2013-08-23 17:59:16Z mushkar $
#ifndef __FunComptonTheory_h__
#define __FunComptonTheory_h__

#include <sstream>

#include "TRandom3.h"

#include "FitTypes.h"
#include "TFit.h"
#include "TFitFunc.h"
#include "TFitFuncLinApprox.h"

namespace Fit {

//_____________________________________________
// This should be transformed into some plugin

//_____________________________________________
// Calculate Sigma2x in a given point using Pasquini's program
Double_t Sigma2xPasquini(const Fit::Point_t &p)
{
  //
  stringstream ss;
  ss << "sh $Fit/scripts/runRCSDRTOT2011_Obs7_FitAllSpinPol ";
  
  // Set parameters
  ss << "2 ";
  // Eg in GeV
  ss << p.at(1)*0.001 << " ";
  // Theta_min, Theta_max
  ss << p.at(0) << " " << p.at(0) << " ";
  // Theta_step
  ss << p.at(0)+1. << " ";
  // alpha + beta
  ss << p.at(2)+p.at(3) << " ";
  // alpha - beta
  ss << p.at(2)-p.at(3) << " ";
  // gammas
  for (UInt_t i=4; i<p.size(); ++i) ss << p.at(i) << " ";
  // tmp dir name
  stringstream ssdir;
  ssdir << gRandom->Integer(1000000);
  TString dir = ssdir.str();
  ss  << dir << " > /dev/null 2>&1";
  
  // Run the script
  system(ss.str().c_str());
  
  // Read data
  Fit::Data_t data;
  Fit::TFit::ReadDataFile("tmp/"+dir+"/output.dat",data);
  
  // Remove tmp dir
  system("rm -rf /dev/shm/fit-tmp/"+dir);
  
  // Calculate Sigma2x
  if (data.empty()) return kUndefined;
  if (data.at(0).size() != 4) return kUndefined;
  return (data.at(0).at(2)-data.at(0).at(3))/(data.at(0).at(2)+data.at(0).at(3));
}

//_____________________________________________
// Calculate Sigma3 in a given point using Pasquini's program
Double_t Sigma3Pasquini(const Fit::Point_t &p)
{
  //
  stringstream ss;
  ss << "sh $Fit/scripts/runRCSDRTOT2011_Obs10_FitAllSpinPol ";
  
  // Set parameters
  // Eg in GeV
  ss << p.at(1)*0.001 << " ";
  // Theta_min, Theta_max
  ss << p.at(0) << " " << p.at(0) << " ";
  // Theta_step
  ss << p.at(0)+1. << " ";
  // alpha + beta
  ss << p.at(2)+p.at(3) << " ";
  // alpha - beta
  ss << p.at(2)-p.at(3) << " ";
  // gammas
  for (UInt_t i=4; i<p.size(); ++i) ss << p.at(i) << " ";
  // tmp dir name
  stringstream ssdir;
  ssdir << gRandom->Integer(1000000);
  TString dir = ssdir.str();
  ss  << dir << " > /dev/null 2>&1";
  
  // Run the script
  system(ss.str().c_str());
  
  // Read data
  Fit::Data_t data;
  Fit::TFit::ReadDataFile("tmp/"+dir+"/output.dat",data);
  
  // Remove tmp dir
  system("rm -rf /dev/shm/fit-tmp/"+dir);
  
  // Sigma3
  if (data.empty()) return kUndefined;
  if (data.at(0).size() != 3) return kUndefined;
  return data.at(0).at(1);
}

//_____________________________________________
// Calculate diff c.s. in a given point using Pasquini's program
Double_t CrossSectionPasquini(const Fit::Point_t &p)
{
  //
  stringstream ss;
  ss << "sh $Fit/scripts/runRCSDRTOT2011_Obs1_FitAllSpinPol ";
  
  // Set parameters
  // Theta
  ss << p.at(0) << " ";
  // e_min, e_max
  ss << p.at(1)*0.001 << " " << p.at(1)*0.001 << " ";
  // e_step
  ss << p.at(1) << " ";
  // alpha + beta
  ss << p.at(2)+p.at(3) << " ";
  // alpha - beta
  ss << p.at(2)-p.at(3) << " ";
  // gammas
  for (UInt_t i=4; i<p.size(); ++i) ss << p.at(i) << " ";
  // tmp dir name
  stringstream ssdir;
  ssdir << gRandom->Integer(1000000);
  TString dir = ssdir.str();
  ss  << dir << " > /dev/null 2>&1";
  
  // Run the script
  system(ss.str().c_str());
  
  // Read data
  Fit::Data_t data;
  Fit::TFit::ReadDataFile("tmp/"+dir+"/output.dat",data);
  
  // Remove tmp dir
  system("rm -rf /dev/shm/fit-tmp/"+dir);
  
  // Cross section
  if (data.empty()) return kUndefined;
  if (data.at(0).size() < 2) return kUndefined;
  return data.at(0).at(1);
}

//_____________________________________________
// Calculate Sigma2x in a given point using Pascalutsa's program
Double_t Sigma2xPascalutsa(const Fit::Point_t &p)
{
  //
  stringstream ss;
  ss << "sh $Fit/scripts/runComptonEFT ";
  
  // Set parameters
  for (UInt_t i=0; i<p.size(); ++i) ss << p.at(i) << " ";
  
  // tmp dir name
  stringstream ssdir;
  ssdir << gRandom->Integer(1000000);
  TString dir = ssdir.str();
  ss << dir << " > /dev/null 2>&1";
  
  // Run the script
  system(ss.str().c_str());
  
  // Read data
  Fit::Data_t data;
  Fit::TFit::ReadDataFile("tmp/"+dir+"/xSigma2x_Adistribution.output",data);
  
  // Remove tmp dir
  system("rm -rf /dev/shm/fit-tmp/"+dir);
  
  // Sigma2x
  if (data.empty()) return kUndefined;
  if (data.at(0).size() < 2) return kUndefined;
  return data.at(0).at(1);
}

//_____________________________________________
// Calculate Sigma3 in a given point using Pascalutsa's program
Double_t Sigma3Pascalutsa(const Fit::Point_t &p)
{
  //
  stringstream ss;
  ss << "sh $Fit/scripts/runComptonEFT ";
  
  // Set parameters
  for (UInt_t i=0; i<p.size(); ++i) ss << p.at(i) << " ";
  
  // tmp dir name
  stringstream ssdir;
  ssdir << gRandom->Integer(1000000);
  TString dir = ssdir.str();
  ss << dir << " > /dev/null 2>&1";
  
  // Run the script
  system(ss.str().c_str());
  
  // Read data
  Fit::Data_t data;
  Fit::TFit::ReadDataFile("tmp/"+dir+"/xBa_Adistribution.output",data);
  
  // Remove tmp dir
  system("rm -rf /dev/shm/fit-tmp/"+dir);
  
  // Sigma3
  if (data.empty()) return kUndefined;
  if (data.at(0).size() < 2) return kUndefined;
  return -0.01*data.at(0).at(1);
}

//_____________________________________________
// Calculate diff. cross section in a given point using Pascalutsa's program
Double_t CrossSectionPascalutsa(const Fit::Point_t &p)
{
  //
  stringstream ss;
  ss << "sh $Fit/scripts/runComptonEFT ";
  
  // Set parameters
  for (UInt_t i=0; i<p.size(); ++i) ss << p.at(i) << " ";
  
  // tmp dir name
  stringstream ssdir;
  ssdir << gRandom->Integer(1000000);
  TString dir = ssdir.str();
  ss << dir << " > /dev/null 2>&1";
  
  // Run the script
  system(ss.str().c_str());
  
  // Read data
  Fit::Data_t data;
  Fit::TFit::ReadDataFile("tmp/"+dir+"/xcs_Adistribution.output",data);
  
  // Remove tmp dir
  system("rm -rf /dev/shm/fit-tmp/"+dir);
  
  // diff cross section
  if (data.empty()) return kUndefined;
  if (data.at(0).size() < 3) return kUndefined;
  return data.at(0).at(1);
}

//_____________________________________________
// Run Lvov's program to get cross-section and Sigma3 as a function eg, theta, alpha and beta
void RunLvov(const Fit::Point_t &p, Fit::Data_t &data)
{
  // Set parameters
  stringstream ss;
  for (UInt_t i=0; i<5; ++i) ss << p.at(i) << " ";
  
  // tmp dir name
  stringstream ssdir;
  ssdir << gRandom->Integer(1000000);
  TString dir = ssdir.str();
  ss << dir << " > /dev/null 2>&1";
  
  // mkdir tmp
  system("mkdir -p /dev/shm/fit-tmp/"+dir);
  
  // Run
  const TString command = "cd $ComptonLvov/gngn_code; ./calc_cs_asymm "+ss.str();
  system(command);
  
  // Read data
  Fit::TFit::ReadDataFile("tmp/"+dir+"/lvov_cs_asymm.out",data);
  
  // Remove tmp dir
  system("rm -rf /dev/shm/fit-tmp/"+dir);
}

//_____________________________________________
// Cross section Lvov's
inline Double_t CrossSectionLvov(const Fit::Point_t &p)
{
  //
  Fit::Data_t data;
  RunLvov(p,data);
  
  // Cross section
  if (data.empty()) return kUndefined;
  if (data.at(0).size() < 4) return kUndefined;
  return data.at(0).at(2);
}

//_____________________________________________
// Lvov's Sigma3
inline Double_t Sigma3Lvov(const Fit::Point_t &p)
{
  //
  Fit::Data_t data;
  RunLvov(p,data);
  
  // Sigma3
  if (data.empty()) return kUndefined;
  if (data.at(0).size() < 4) return kUndefined;
  return data.at(0).at(3);
}

//_____________________________________________
// The Baldin sum rule constraint on the alpha and beta
inline Double_t ConstraintBaldinSumRule(const Double_t *p)
{
  //
  return Fit::TFit::Chi2(p[0]+p[1],0.4,13.8);
}

//_____________________________________________
// alpha-beta from McGovern, Phillips, Grießhammer, EPJ A (2013) 49, p.12,
inline Double_t ConstraintAlphaMinusBetaMPG(const Double_t *p)
{
  //
  return Fit::TFit::Chi2(p[0]-p[1],1.7,7.6);
}

//_____________________________________________
// Constraint on the alpha from the Pascalutsa's theory
inline Double_t ConstraintAlphaPascalutsa(const Double_t *p)
{
  //
  return Fit::TFit::Chi2(p[0],0.7,11.2461);
}

//_____________________________________________
// Constraint on the beta from the Pascalutsa's theory
inline Double_t ConstraintBetaPascalutsa(const Double_t *p)
{
  //
  return Fit::TFit::Chi2(p[1],0.7,3.87054);
}

//_____________________________________________
// Constraint on the alpha from the PDG
inline Double_t ConstraintAlphaPDG(const Double_t *p)
{
  //
  return Fit::TFit::Chi2(p[0],0.6,12.);
}

//_____________________________________________
// Constraint on the beta from the PDG
inline Double_t ConstraintBetaPDG(const Double_t *p)
{
  //
  return Fit::TFit::Chi2(p[1],0.5,1.9);
}

//_____________________________________________
// Constraint on the gamma0 from the previous experiments
inline Double_t ConstraintGamma0(const Double_t *p)
{
  //
  Double_t gamma0 = - p[2] - p[3] - p[4] - p[5];
  return Fit::TFit::Chi2(gamma0,0.18,-1.01);
}

//_____________________________________________
// Constraint on the gamma_pi from the previous experiments
inline Double_t ConstraintGammaPi(const Double_t *p)
{
  //
  Double_t gammaPi = - p[2] + p[3] - p[4] + p[5];
  return Fit::TFit::Chi2(gammaPi,1.8,8.);
}

//_____________________________________________
// Convert from apb amb ge1e1 gm1m1 g0 gpi => to the standard set
inline void ConvertPolariz1(const Double_t *pin, Double_t *pout)
{
  pout[0] = (pin[0]+pin[1])/2.; // alpha
  pout[1] = (pin[0]-pin[1])/2.; // beta
  pout[2] = pin[2]; // e1e1
  pout[3] = pin[3]; // m1m1
  pout[4] = (-pin[4]-pin[5])/2. - pout[2]; // e1m2
  pout[5] = (-pin[4]+pin[5])/2. - pout[3]; // m1e2

}

//_____________________________________________
// Convert from apb amb gpi => a b gpi
inline void ConvertPolariz2(const Double_t *pin, Double_t *pout)
{
  pout[0] = (pin[0]+pin[1])/2.; // alpha
  pout[1] = (pin[0]-pin[1])/2.; // beta
  pout[2] = pin[2]; // gpi
}

//_____________________________________________
// Convert from a b ge1e1 gm1m1 g0 gpi => to the standard set
inline void ConvertPolariz3(const Double_t *pin, Double_t *pout)
{
  pout[0] = pin[0]; // alpha
  pout[1] = pin[1]; // beta
  pout[2] = pin[2]; // e1e1
  pout[3] = pin[3]; // m1m1
  pout[4] = (-pin[4]-pin[5])/2. - pout[2]; // e1m2
  pout[5] = (-pin[4]+pin[5])/2. - pout[3]; // m1e2
}

//_____________________________________________
// Adding the Compton function to the TFit::fListFitFun and fListConstraints
inline void InitComptonTheory()
{  
  // Check if it's already initialized
  static Bool_t isInit = kFALSE;
  if (isInit) return;
  isInit = kTRUE;
  
  // Init random gen.
  gRandom->SetSeed(0);
  
  // Create the fitting functions
  
  // Pascalutsa
  static Fit::TFitFunc funPascalutsaSigma2x("PascalutsaSigma2x","Fitting function PascalutsaSigma2x",Fit::Sigma2xPascalutsa);
  static Fit::TFitFunc funPascalutsaSigma3("PascalutsaSigma3","Fitting function PascalutsaSigma3",Fit::Sigma3Pascalutsa);
  static Fit::TFitFunc funPascalutsaCS("PascalutsaCS","Fitting function PascalutsaCS",Fit::CrossSectionPascalutsa);
  //
  static Fit::TFitFuncLinApprox funPascalutsaSigma2xLinApprox("PascalutsaSigma2xLinApprox","Fitting function PascalutsaSigma2x (lin. approx)");
  funPascalutsaSigma2xLinApprox.SetFitFunc(&funPascalutsaSigma2x);
  static Fit::TFitFuncLinApprox funPascalutsaSigma3LinApprox("PascalutsaSigma3LinApprox","Fitting function PascalutsaSigma3 (lin. approx)");
  funPascalutsaSigma3LinApprox.SetFitFunc(&funPascalutsaSigma3);
  static Fit::TFitFuncLinApprox funPascalutsaCSLinApprox("PascalutsaCSLinApprox","Fitting function PascalutsaCS (lin. approx)");
  funPascalutsaCSLinApprox.SetFitFunc(&funPascalutsaCS);
  
  // Pasquini
  static Fit::TFitFunc funPasquiniSigma2x("PasquiniSigma2x","Fitting function PasquiniSigma2x",Fit::Sigma2xPasquini);
  static Fit::TFitFunc funPasquiniSigma3("PasquiniSigma3","Fitting function PasquniSigma3",Fit::Sigma3Pasquini);
  static Fit::TFitFunc funPasquiniCS("PasquiniCS","Fitting function PasquniCS",Fit::CrossSectionPasquini);
  //
  static Fit::TFitFuncLinApprox funPasquiniSigma2xLinApprox("PasquiniSigma2xLinApprox","Fitting function PasquiniSigma2x (lin. approx)");
  funPasquiniSigma2xLinApprox.SetFitFunc(&funPasquiniSigma2x);
  static Fit::TFitFuncLinApprox funPasquiniSigma3LinApprox("PasquiniSigma3LinApprox","Fitting function PasquiniSigma3 (lin. approx)");
  funPasquiniSigma3LinApprox.SetFitFunc(&funPasquiniSigma3);
  static Fit::TFitFuncLinApprox funPasquiniCSLinApprox("PasquiniCSLinApprox","Fitting function PasquiniCS (lin. approx)");
  funPasquiniCSLinApprox.SetFitFunc(&funPasquiniCS);
  
  // Lvov
  static Fit::TFitFunc funLvovSigma3("LvovSigma3","Fitting function LvovSigma3",Fit::Sigma3Lvov);
  static Fit::TFitFunc funLvovCS("LvovCS","Fitting function LvovCS",Fit::CrossSectionLvov);
  //
  static Fit::TFitFuncLinApprox funLvovSigma3LinApprox("LvovSigma3LinApprox","Fitting function LvovSigma3 (lin. approx)");
  funLvovSigma3LinApprox.SetFitFunc(&funLvovSigma3);
  static Fit::TFitFuncLinApprox funLvovCSLinApprox("LvovCSLinApprox","Fitting function LvovCS (lin. approx)");
  funLvovCSLinApprox.SetFitFunc(&funLvovCS);
      
  // Add possible fitting functions
  TFit::AddFitFunc(&funPascalutsaSigma2x);
  TFit::AddFitFunc(&funPascalutsaSigma3);
  TFit::AddFitFunc(&funPascalutsaCS);
  TFit::AddFitFunc(&funPasquiniSigma2x);
  TFit::AddFitFunc(&funPasquiniSigma3);
  TFit::AddFitFunc(&funPasquiniCS);
  TFit::AddFitFunc(&funLvovSigma3);
  TFit::AddFitFunc(&funLvovCS);
  //
  TFit::AddFitFunc(&funPascalutsaSigma2xLinApprox);
  TFit::AddFitFunc(&funPascalutsaSigma3LinApprox);
  TFit::AddFitFunc(&funPascalutsaCSLinApprox);
  TFit::AddFitFunc(&funPasquiniSigma2xLinApprox);
  TFit::AddFitFunc(&funPasquiniSigma3LinApprox);
  TFit::AddFitFunc(&funPasquiniCSLinApprox);
  TFit::AddFitFunc(&funLvovSigma3LinApprox);
  TFit::AddFitFunc(&funLvovCSLinApprox);
  
  // Add possible constraints
  TFit::AddConstraint("BaldinSumRule",&Fit::ConstraintBaldinSumRule);
  TFit::AddConstraint("AlphaPascalutsa",&Fit::ConstraintAlphaPascalutsa);
  TFit::AddConstraint("BetaPascalutsa",&Fit::ConstraintBetaPascalutsa);
  TFit::AddConstraint("AlphaPDG",&Fit::ConstraintAlphaPDG);
  TFit::AddConstraint("BetaPDG",&Fit::ConstraintBetaPDG);
  TFit::AddConstraint("Gamma0",&Fit::ConstraintGamma0);
  TFit::AddConstraint("GammaPi",&Fit::ConstraintGammaPi);
  TFit::AddConstraint("AlphaMinusBeta_MPG",&Fit::ConstraintAlphaMinusBetaMPG);
  
  // Add converters
  TFit::AddConverter("apb_amb_ge1e1_gm1m1_g0_gpi",&Fit::ConvertPolariz1);
  TFit::AddConverter("apb_amb_gpi",&Fit::ConvertPolariz2);
  TFit::AddConverter("a_b_ge1e1_gm1m1_g0_gpi",&Fit::ConvertPolariz3);
  
}

}

#endif
