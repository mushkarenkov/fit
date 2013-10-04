// SVN: $Id: PartDerivMain.cxx 244 2013-09-25 13:22:38Z mushkar $

// STD
#include <getopt.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
using namespace std;
// OpenMP
// #include <omp.h>
// ROOT
#include <TString.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TFrame.h>
#include <TGraphErrors.h>
#include <TLatex.h>
#include <TApplication.h>
// Fit
#include "TFit.h"
#include "FitTypes.h"
#include "TFitFunc.h"
#include "FunComptonTheory.h"

void Usage(const char*);
void ParseCommandLineArguments(int, char**,TString&);
void PrintPartialDerivatives(Fit::TFit&,const TString&,const TString&,const TString&,Fit::Point_t&,Double_t*);

void CalcPartDeriv(Fit::TFit&);

static const Double_t kInfiniti = 10e+10;

//_____________________________________________
int main(int argc, char** argv)
{
  // Parse the command line arguments
  TString nameFileParam;
  ParseCommandLineArguments(argc,argv,nameFileParam);
  
  // Create the fitting functions and constraints from the compton theory
  Fit::InitComptonTheory();
  
  // Create a fitter
  Fit::TFit fit(nameFileParam);
  
  //
  CalcPartDeriv(fit);
  
  //
  return 0;
}

//_____________________________________________
void Usage(const char *nameProgram)
{
  // Print help
  cout << endl;
  cout << "Usage: " << nameProgram << " --param=file [--help]" << endl;
  cout << endl;
  cout << "Options: " << endl;
  cout << "\t-h  --help  \t print this help and exit" << endl;
  cout << "\t-p  --param \t a file with parameters" << endl;
  cout << endl;
  
}

//_____________________________________________
void ParseCommandLineArguments(int argc, char** argv, TString &nameFileParam)
{
  // Parse the command line arguments
  
  // When no option is provided => exit
  if (argc == 1)
  {
    cout << ">>> ERROR! No option provided!" << endl;
    Usage(argv[0]);
    exit(-1);
  }
  
  // Define options
  const char *optsShort = "hp:";
  const struct option optsLong[] = {
    {"help",  no_argument      ,NULL,'h'},
    {"param", required_argument,NULL,'p'},
    {NULL,   0                 ,NULL, 0 }
  };
  
  // Parse options
  Int_t rez, iOpt = -1;
  while ( (rez=getopt_long(argc,argv,optsShort,optsLong,&iOpt)) != -1 )
  {
    switch (rez)
    {
      case 'h':
	Usage(argv[0]);
	exit(0);
	break;
      case 'p':
	nameFileParam = optarg;
	break;
      case '?':
      default:
	exit(-1);
    }
  }
  
}

//_____________________________________________
//_____________________________________________
//_____________________________________________
void PrintPartialDerivatives(Fit::TFit &fit, const TString &obs, const TString &expset, const TString &suffix, Fit::Point_t &point)
{
  //
  ofstream fileOut(obs+"_"+expset+"_"+suffix+".dat");
  //
  Fit::DataSetExp_t::const_iterator exp = fit.GetDataExp(obs,expset);
  if (exp == fit.GetDataExp()->end())
  {
    cout << ">>> ERROR: " << obs << " " << expset << " not found!" << endl;
    return;
  }
  const Fit::Data_t &dataExp = exp->fData;
  const Fit::TFitFunc *dataTh = exp->fTh;

  //
  fileOut << "Alpha = " << point[2] << endl;
  fileOut << "Beta = " << point[3] << endl;
  fileOut << "gE1E1 = " << point[4] << endl;
  fileOut << "gM1M1 = " << point[5] << endl;
  fileOut << "gE1M2 = " << point[6] << endl;
  fileOut << "gM1E2 = " << point[7] << endl;
  fileOut << endl;
//   fileOut << "sigma_g0 = " <<  "0.18 x10^-4 fm^4" << endl;
//   fileOut << "sigma_gpi = " <<  "1.8 x10^-4 fm^4" << endl;
  fileOut << endl;
  fileOut << "A = " << obs << endl;
  fileOut << endl;
  fileOut << "Eg\ttheta\texp\texp_err\t   A\t\tdA/dAlpha\tdA/dBeta\tdA/dgE1E1\tdA/dgM1M1\tdA/dgE1M2\tdA/dgM1E2" << endl;
  //
  Int_t iPoint = 1;
  for (Fit::Data_t::const_iterator iDataExp = dataExp.begin(); iDataExp != dataExp.end(); ++iDataExp)
  {
    Int_t i;
    for (i=0; i<2; ++i) point[i] = iDataExp->at(i);
    Double_t observ    = iDataExp->at(i);   // value
    Double_t observErr = iDataExp->at(i+1); // err_value
    cout << "exp point "<< iPoint << ":\t" << "e: "<< point[1] << "\ttheta: "<< point[0] << 
    "\ta: " << point[2] << 
    "\tb: "<< point[3] << 
    "\tee: "<< point[4] << 
    "\tmm: "<< point[5] << 
    "\tem: "<< point[6] << 
    "\tme: "<< point[7];
    //
//     point[1] = 50.; // MeV
    //
    Double_t observTh = dataTh->GetValue(point);
    // deriv
    Double_t deriv[6];
    for (Int_t i=0; i<6;++i)
    {
      Fit::Point_t point2 = point;
      Int_t ii = i+2;
      point2[ii] += 0.1;
      Double_t observTh2 = dataTh->GetValue(point2);
      deriv[i] = (observTh2-observTh)/(point2[ii]-point[ii]);
    }
    //
    cout << "\texp: "<< observ << "\tth: " << observTh << "\texp_err: " << observErr << endl;
    cout << "\td_alpha: "<< deriv[0] << "\td_beta: " << deriv[1] << endl;
    ++iPoint; // only for TEST
    //
    fileOut <<
    point[1] << "\t" << point[0] <<
    "\t" << observ << "\t" << observErr <<
    "\t" << observTh <<
    "\t" << deriv[0] <<
    "\t" << deriv[1] <<
    "\t" << deriv[2] <<
    "\t" << deriv[3] <<
    "\t" << deriv[4] <<
    "\t" << deriv[5] << endl;
  }
  
  //
  fileOut.close();
}

//_____________________________________________
void CalcPartDeriv(Fit::TFit &fit)
{
  //
  Fit::Point_t point1(8);
  Double_t point[8], point_conv[8];
  
  //
  point[0] =  13.99;
  point[1] =  7.37;
  point[2] = -3.12;
  point[3] =  3.39;
  point[4] = -1.05;
  point[5] =  9.01;
  
  //
  Fit::ConvertPolariz1(point,point_conv);
  for (Int_t i=0;i<6;++i) point1[i+2] = point_conv[i];
  
  //
//   PrintPartialDerivatives(fit,"Sigma2x_Pascalutsa","mar","1",point1);
//   PrintPartialDerivatives(fit,"Sigma2x_Pascalutsa","mar","2",point2);
//   PrintPartialDerivatives(fit,"Sigma2x_Pascalutsa","mar","3",point3);
//   PrintPartialDerivatives(fit,"Sigma3_Pascalutsa","bla01","1",point1);
//   PrintPartialDerivatives(fit,"Sigma3_Pascalutsa","bla01","2",point2);
//   PrintPartialDerivatives(fit,"Sigma3_Pascalutsa","bla01","3",point3);
	
  PrintPartialDerivatives(fit,"Sigma2x_Pasquini","mar","1",point1);
  PrintPartialDerivatives(fit,"Sigma3_Pasquini","bla01","1",point1);
//   PrintPartialDerivatives(fit,"Sigma2x_Pasquini","mar","2",point2);
//   PrintPartialDerivatives(fit,"Sigma3_Pasquini","bla01","2",point2);
//   PrintPartialDerivatives(fit,"Sigma2x_Pasquini","mar","3",point3);
//   PrintPartialDerivatives(fit,"Sigma3_Pasquini","bla01","3",point3);

}
