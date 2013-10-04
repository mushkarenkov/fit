// SVN: $Id: TFit.cxx 224 2013-07-27 23:25:18Z mushkar $
#include <fstream>
#include <iostream>
#include <sstream>

#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TFrame.h>
#include <TVectorD.h>

#include "TFit.h"
#include "TFitFuncLinApprox.h"
using namespace Fit;

ClassImp(TFit)

//__________________________________________________________________
//
// Describes object that does the fit
//

// Initialize static data members
vector<TFitFunc*> TFit::fListFitFun;
map<TString,Constraint_t> TFit::fListConstraints;
map<TString,Converter_t> TFit::fListConverters;

//__________________________________________________________________
// TFit::TFit(const char* appClassName, Int_t* argc, char** argv, const TString &nameFileParam) : TRint(appClassName,argc,argv)
TFit::TFit(const TString &nameFileParam) : fNparam(0)
{
  //
  
  // Init data-members
  InitDataMembers();  
  
  // Read parameters
  ReadParams(nameFileParam);

}

//__________________________________________________________________
TFit::TFit() : fNparam(0)
{
  //
  InitDataMembers();
}

//__________________________________________________________________
void TFit::InitDataMembers()
{
  //
  fMinimizer = NULL;
  fResult = new TFitResult("Result", "Result of the minimization");
  fInputPar.SetNameTitle("InputPar","");
  fConvertParams = NULL;
}

//__________________________________________________________________
TFit::~TFit()
{
  // Destructor
  
  //
  delete fResult;
}

//__________________________________________________________________
void TFit::ReadParams(const TString &nameFileParam)
{
  // Reads an input configuration file
  
  // Open the configuration file
  ifstream fileIn(nameFileParam);
  
  // Check if the file is opened
  if (!fileIn.is_open())
  {
    cout << ">>> ERROR! TFit::ReadParams(): Can't open the file " << nameFileParam << endl;
    return;
  }
  
  // Reading
  string line;
  TString field, keyword;
  vector<TString> record;
  vector<Double_t> parInit, parStep;
  vector<TString> parFix;
  while ( getline(fileIn,line) )
  {
    stringstream ss(line);
    record.clear();
    while ( ss.good() )
    {
      ss >> field;
      record.push_back(field);
    }
    // Check keyword
    keyword = record.at(0);
    //
    if (keyword == "Exp")
    {
      // Reading an experimental data file
      ReadDataExp(record);
      // Store input parameters
      StoreInputPar(record);
    }
    //
    else if (keyword == "Th")
    {
      // Reading an theoretical data file
      ReadDataTh(record);
      // Store input parameters
      StoreInputPar(record);
    }
    //
    else if (keyword == "Par")
    {
      // Add a new parameter
      Fit::Par_t p =
      {
	record.at(1),
	(record.at(2) == "f"),
	record.at(3).Atof(),
	record.at(4).Atof(),
	record.at(5).Atof(),
	record.at(6).Atof()
      };
      fParams.push_back(p);
      // Store input parameters
      StoreInputPar(record);
    }
    //
    else if (keyword == "ParConverter")
    {
      map<TString,Converter_t>::iterator iter = fListConverters.find(record.at(1));
      if (iter == fListConverters.end())
      {
	cout << ">>> WARNING! ParConverter \"" << record.at(1) << "\" not found!" << endl;
      }
      else
      {
	cout << "Parameters converter: \"" << record.at(1) << "\"" << endl;
	fConvertParams = iter->second;
      }
      // Store input parameters
      StoreInputPar(record);
    }
    //
    else if (keyword == "Constraint")
    {
      map<TString,Constraint_t>::iterator iter = fListConstraints.find(record.at(1));
      if (iter == fListConstraints.end())
      {
	cout << ">>> WARNING! Constraint \"" << record.at(1) << "\" not found!" << endl;
      }
      else
      {
	cout << "Adding constraint \"" << record.at(1) << "\"" << endl;
	fConstraints.push_back(iter->second);
      }
      // Store input parameters
      StoreInputPar(record);
    }
    //
    else if (keyword == "Output")
    {
      fNameFileResult = record.at(1);
      // Store input parameters
      StoreInputPar(record);
    }
    //
    else if (keyword == "Minimizer")
    {
      TString name = record.at(1);
      TString algo = record.at(2);
      if (name == "0") name = "Minuit2";
      if (algo == "0") algo = "";
      fMinimizer = ROOT::Math::Factory::CreateMinimizer(name.Data(), algo.Data());
      fMinimizer->SetMaxFunctionCalls(record.at(3).Atoi()); // for Minuit/Minuit2
      fMinimizer->SetMaxIterations(record.at(4).Atoi());  // for GSL
      fMinimizer->SetTolerance(record.at(5).Atof());
      fMinimizer->SetPrintLevel(record.at(6).Atoi());
      //
      fMinimize = &TFit::MinimizeNormally;
      if (record.size()>7 && record.at(7) == "LinApprox") fMinimize = &TFit::MinimizeLinApprox;
      // Store input parameters
      StoreInputPar(record);
    }
  }
  
  // Check
  if (fDataExp.empty())
  {
    cout << ">>> ERROR! No 'Exp' data is set!" << endl;
    exit(1);
  }
  if (fDataTh.empty())
  {
    cout << ">>> ERROR! No 'Th' data is set!" << endl;
    exit(1);
  }
  
  // Associate exp -> theory
  for (DataSetExp_t::iterator iDataFile = fDataExp.begin(); iDataFile != fDataExp.end(); ++iDataFile) // For each exp data file
  {
    // Get a pointer to the respective theor. data
    DataSetTh_t::iterator iDataTh = fDataTh.find(iDataFile->fKeyObserv);
    if (iDataTh == fDataTh.end())
    {
      cout << ">>> WARNING! '" << iDataFile->fKeyObserv << "' from '" << iDataFile->fKeyData << "' does not have its respective theoretical data to fit!" << endl;
      iDataFile->fTh = NULL; // no associated theoretical data
      continue;
    }
    // Set the pointer
    iDataFile->fTh = iDataTh->second;
  }
  
  // Set N params
  // (NOTE: It's assumed that all exp and theor functions depends on the same variables. Nparam = NvarTheor - NvarExp)
  fNvar = fDataExp.begin()->fData.begin()->size() - 2;
  if (fNvar < 0)
  {
    cout << "ERROR! fNvar < 0" << endl;
    exit(1);
  }
//   fNparam = fDataTh.begin()->second.fX.size() - fNvar; // TODO
  fNparam = fParams.size();
  if (fNparam < 0 )
  {
    cout << "ERROR! fNparam < 0" << endl;
    exit(1);
  }
  cout << "Nvar = " << fNvar << "\t" << "Nparam = " << fNparam << endl;
  
}

//__________________________________________________________________
void TFit::ReadDataExp(const vector<TString> &record)
{
  // Read an experimental data file
  
  Data_t data;
  TString keyObserv = record.at(1);
  TString keyData   = record.at(2);
  TString nameFile  = record.at(3);
  Double_t syserr   = record.at(4).Atof();
  if ( !ReadDataFile(nameFile,data) ) return;
  
  // Check
  for (DataSetExp_t::const_iterator iData = fDataExp.begin(); iData != fDataExp.end(); ++iData)
  {
//     if (iData->fKeyObserv != keyObserv) continue; // NOTE: At the moment, all 'Exp' data depend on the same variables
    if (iData->fData.at(0).size() != data.at(0).size())
    {
      cout << ">>> ERROR! '" << iData->fKeyObserv << "' from '" << iData->fKeyData << "' and '" << keyData << "' have different number of variables!" << endl;
      exit(1);
    }
    else
    {
      break;
    }
  }
  
  // Adding to the experimental data list
  DataExp_t dataExp = {keyObserv,keyData,data,NULL,syserr};
  fDataExp.push_back( dataExp );
}

//__________________________________________________________________
void TFit::ReadDataTh(const vector<TString> &record)
{
  // Read a theoretical data with discrete function to fit
  
  // Reading
  Data_t data;
  TString nameFile = record.at(3);
  TString key = record.at(1);
  TString type = record.at(2);
  if (type == "func")
  {
    for (UInt_t i=0; i<fListFitFun.size(); ++i)
    {
      if (static_cast<TString>(fListFitFun.at(i)->GetName()) == nameFile)
      {
	fDataTh[key] = fListFitFun.at(i);
      }
      else
      {
// 	cout << ">>> ERROR:FitFun " << nameFile << "not fount!" << endl;
      }
    }
  }
  else
  {
//     if ( !ReadDataFile(nameFile,data) ) return;
//     
//     // Check
//     for (DataSetTh_t::const_iterator iData = fDataTh.begin(); iData != fDataTh.end(); ++iData)
//     {
//       if (iData->first == key)
//       {
// 	cout << ">>> WARNING! Theoretical data for '" << key << "' is already set. '" << nameFile << "' is ignored!" << endl;
// 	return;
//       }
//     }
//     // Check that all the theoretical data (functions) have the same number of variables
//     // NOTE: At the moment, all 'Th' data depend on the same variables
//     if (!fDataTh.empty() && fDataTh.begin()->second.fX.size() != data.size()-1 )
//     {
//       cout << ">>> ERROR! '" << key << "' from '" << nameFile << "' has different number of variables!" << endl;
//       exit(1);
//     }
//     
//     // Adding to the theoretical data
//     if (key == "Sigma2xx")
//     {
//       // NOTE: Temporary solution for Sigma2x
//       TFitDiscreteFunc f;
//       f.SetData2(data);
//       fDataTh[key] = f;
//     }
//     else
//     {
//       fDataTh[key] = new TFitDiscreteFunc(data); // TODO cleanup memory
//     }
  }
}

//__________________________________________________________________
Bool_t TFit::ReadDataFile(const TString &nameFile, Data_t &data)
{
  // Read a data file
  
  // Open the file
  ifstream fileIn(TString(getenv("FitData"))+"/"+nameFile);
  
  // Check if the file is opened
  if (!fileIn.is_open())
  {
    cout << ">>> ERROR! TFit::ReadDataFile(): Can't open the file " << nameFile << endl;
    return kFALSE;
  }
  
  // Read the file
  Double_t field;
  Point_t point;
  data.clear();
  string line;
  while ( getline(fileIn,line) )
  {
    stringstream ss(line); // TODO move stringstream ss outside the loop
    point.clear();
    while ( ss.good() )
    {
      ss >> field;
      point.push_back(field);
    }
    data.push_back(point);
  }
  
  // Close the file
  fileIn.close();
  
  return kTRUE;
}

//__________________________________________________________________
void TFit::MinimizeLinApprox(const Char_t *minName, const Char_t *algoName)
{
  // Newton's method...
  
  // Convert params to standard ones
  Double_t pin[fNparam], parstd[fNparam];
  for (Int_t i=0; i<fNparam; ++i) pin[i] = fParams[i].fInit;
  ConvertParams(pin,parstd);
  // Set exp coord
  set<TFitFuncLinApprox*> setTheor;
  for (UInt_t i=0; i<fDataExp.size(); ++i) // for each set
  {
     // Get pointer to the respective theor. data
    TFitFuncLinApprox *dataTh = static_cast<TFitFuncLinApprox*>(fDataExp.at(i).fTh);
    // Check if there is corresponding theory to fit
    if (!dataTh)
    {
      cout << "ERROR >>> TFit::Minimize2(): \'" << fDataExp.at(i).fKeyObserv << "\'\t" << fDataExp.at(i).fKeyData << " does not have assosiated TFitFuncLinApprox theory!" << endl;
      return;
    }
    // Add to the set of uniq theories used in the fit, if added set params
    if (setTheor.insert(dataTh).second)
    {
      dataTh->SetNvar(fNvar);
      for (Int_t i=0; i<fNparam; ++i) dataTh->AddPar(parstd[i]);
    }
    // Get pointer to the exp. data
    const Data_t *dataExp = &(fDataExp.at(i).fData);
    // Calc the approximation parameters
    for (UInt_t j=0; j<dataExp->size(); ++j) dataTh->AddApproxAt(dataExp->at(j));
  }
  
  // Iterative minimization
  vector<vector<Double_t> > params(fNparam); // stores history of params
  vector<Double_t> chi2; // history of chi2
  for (Int_t i=0; i<fNparam; ++i) params[i].push_back(fParams[i].fInit);
  while (!IsStableMinimum(params,0.01))
  {  
    // minimize
    MinimizeNormally(minName,algoName);
    
    // chi2 history
    chi2.push_back(fMinimizer->MinValue());
    for (UInt_t i=0; i<chi2.size(); ++i) cout << chi2[i] << " ";
    cout << endl;
    
    // if is not converged => return
    if (fMinimizer->Status() != 0)
    {
      cout << "Minimization status: " << fMinimizer->Status() << " => Cannot be minimized!" << endl;
      break;
    }
    
    // Assign the new best guess to the init params values and add to history
    for (Int_t i=0; i<fNparam; ++i)
    {
      pin[i] = fParams[i].fInit = fParams[i].fValue;
      params[i].push_back(pin[i]);
    }
    
    // Convert params to the std
    ConvertParams(pin,parstd);
    
    // Set new best guess for each teor func
    for (set<TFitFuncLinApprox*>::iterator iTh = setTheor.begin(); iTh != setTheor.end(); ++iTh)
    {
      // Set new params
      (*iTh)->SetPar(fNparam,parstd);
      // Re-calculate the liniar approximation
      (*iTh)->ReCalcLinApprox();
    }
  }
  
}

//__________________________________________________________________
Bool_t TFit::IsStableMinimum(const vector<vector<Double_t> > &params, const Double_t &tol) const
{
  // Checks stability of parameters during minimization to find a stable minimum
  
  for (Int_t ipar=0; ipar<fNparam; ++ipar)
  {
    Int_t np = params[ipar].size();
    if (np<3) return kFALSE;
    for (Int_t i=1; i<3; ++i)
    {
      Double_t ddd = TMath::Abs(params[ipar].back()-params[ipar][np-i-1]);
      cout << "Check tol: " << ipar << ":\tiiter: " << np-i-1 << "\t" << ddd << " < " << tol << endl;
      if (ddd > tol) return kFALSE;
    }
  }
  
  //
  return kTRUE;
}

//__________________________________________________________________
void TFit::SetMinimizer(const Char_t *minName, const Char_t *algoName)
{
  // Set parameters of theminimizer (create it if wos not created before)
  
  // Create minimizer
  if (fMinimizer == NULL)
  {
    ROOT::Math::Minimizer *fMinimizer = ROOT::Math::Factory::CreateMinimizer(minName, algoName);
    // set tolerance , etc...                                                                                                                                                                                       
    fMinimizer->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
    fMinimizer->SetMaxIterations(10000);  // for GSL
    fMinimizer->SetTolerance(0.1);
    fMinimizer->SetPrintLevel(1);
  }
  
  // create funciton wrapper for minmizer a IMultiGenFunction type
  Int_t npar = fNparam+fDataExp.size();
  static ROOT::Math::Functor f(&(*this),&Fit::TFit::GlobalChi2,npar);
  fMinimizer->SetFunction(f);
  
  // Set the params to be minimized
  for (Int_t i=0; i<fNparam; ++i)
  {
    if (fParams[i].fFix)
      fMinimizer->SetFixedVariable(i,fParams.at(i).fName.Data(),fParams.at(i).fInit);
    else
      fMinimizer->SetLimitedVariable(i,fParams.at(i).fName.Data(),fParams.at(i).fInit,fParams.at(i).fStep,fParams.at(i).fMin,fParams.at(i).fMax);
  }
  
  // Scale factors
  fResult->fNpointsTot = 0;
  fNsys = 0;
  for (UInt_t i=0; i<fDataExp.size(); ++i)
  {
    //
    fResult->fNpointsTot += fDataExp.at(i).fData.size();
    //
    TString sname = "scale_"+fDataExp.at(i).fKeyData+"_"+fDataExp.at(i).fKeyObserv;
    if (fDataExp.at(i).fSysErr != 0.)
    {
      ++fNsys;
      cout << i+fNparam << endl;
      fMinimizer->SetVariable(i+fNparam,sname.Data(),1., 0.1);
      cout << "1111111" << endl;
    }
    else
    {
      cout << i+fNparam << endl;
      fMinimizer->SetFixedVariable(i+fNparam,sname.Data(),1.);
      cout << "2222222" << endl;
    }
  }
  
}

//__________________________________________________________________
void TFit::MinimizeNormally(const Char_t *minName, const Char_t *algoName)
{
  // Make minimization
  
  // Set minimizer parameters
  SetMinimizer(minName,algoName);
  
  // do the minimization
  fMinimizer->Minimize();
  
  // Canvas for contours
//   fCnvCI = new TCanvas("CnvCI", "Confidence intervals");
//   // NOTE: These 2 strange lines are necessary to create more than 2 canvases
//   TCanvas ccc("ccc", "Dummy");
//   ccc.Delete();
//   //
  
  //   cout << "Minimum: f(" << *(min->X()) << "): " << min->MinValue()  << endl;
  //   
  //   // expected minimum is 0
  //   if ( min->MinValue()  < 1.E-4  && f(xs) < 1.E-4) 
  //     std::cout << "Minimizer " << minName << " - " << algoName 
  //     << "   converged to the right minimum" << std::endl;
  //   else {
    //     std::cout << "Minimizer " << minName << " - " << algoName 
  //     << "   failed to converge !!!" << std::endl;
  //     Error("NumericalMinimization","fail to converge");
  //   }
  
  Bool_t minOK = kTRUE;
  if (minOK)
  {
    SaveResult();
  }
  else
  {
    cout << "Minimization is not OK: => No result!" << endl;
  }
  
}

//__________________________________________________________________
void TFit::SaveResult()
{
  // Output file
  TFile fileResult(fNameFileResult,"RECREATE","TFit output");
  
  // Fill fFitResult
  for (Int_t i=0; i<=fNparam; ++i)
  {
    fParams[i].fValue = fMinimizer->X()[i];
    fParams[i].fError = fMinimizer->Errors()[i];
  }
  fResult->fChi2 = fMinimizer->MinValue();
  fResult->fParams = fParams;
  fResult->fNparam = fNparam;
  fResult->fNconstraints = fConstraints.size();
  fResult->fNparamFree = fMinimizer->NFree();
  
  // Contours
//   const UInt_t kNp = 25;
//   Double_t xContour[kNp], yContour[kNp];
//   for (Int_t i=0; i<fNparam-1; ++i)
//   {
//     for (Int_t j=i+1; j<fNparam; ++j)
//     {
//       UInt_t np = kNp; // Need for ROOT::Math::Minimizer::Contour
//       if (!min->Contour(i,j,np,xContour,yContour)) continue;
//       vector<Double_t> x, y;
//       x.assign(xContour,xContour+kNp);
//       y.assign(yContour,yContour+kNp);
//       Contour_t contour = {i,j,x,y};
//       fResult->fContours.push_back(contour);
//     }
//   }
  
  // Write the result to the output file
  fileResult.cd();
  fResult->Write();
  
  // Write the input parameters to the output file
  fInputPar.SetName("InputPar");
  fInputPar.Write();
  
  // Close the output file
  fileResult.Close();
 
}

//__________________________________________________________________
void TFit::ConvertParams(const Double_t *pin, Double_t *pout)
{
  // Convert fitting parameters to the standard ones
  
  if (fConvertParams)
    fConvertParams(pin,pout);
  else
    for (Int_t i=0; i<fNparam; ++i) pout[i] = pin[i];
    
}

//__________________________________________________________________
Double_t TFit::GlobalChi2(const Double_t *par)
{
  // Global fit
  // Print params to the std output
  for (Int_t ii=0; ii<fNparam; ++ii) cout << fMinimizer->VariableName(ii) << ": " << par[ii] << "\t";
  cout << endl;
  
  // Init chi2
  Double_t chi2 = 0.;
  
  // Convert fitting parameters to the standard ones
  Double_t p[fMinimizer->NDim()];
  ConvertParams(par,p);
  // for the scale factors
  for (UInt_t i=fNparam; i<fMinimizer->NDim(); ++i) p[i] = par[i];
  
  // TEST
  for (Int_t i=0; i<fNparam; ++i) cout << p[i] << "\t";
  cout << endl;
  
  // Constraints
  for (UInt_t i=0; i<fConstraints.size(); ++i)
  {
    Double_t chi2tmp = fConstraints.at(i)(p);
    cout << "constraint" << i << " chi2tmp: " << chi2tmp << endl;
    chi2 += chi2tmp;
  }
  cout << "chi2 constraints: " << chi2 << endl;
  
  // Make the point
  Point_t point(fNvar+fNparam);
  for (Int_t i=0; i<fNparam; ++i)
  {
    point[i+fNvar] = p[i];
  }
  
  // Loop over each exp data file
  for (UInt_t i=0; i<fDataExp.size(); ++i)
  {
    // Check if there is corresponding theory to fit
    if (!fDataExp.at(i).fTh)
    {
      chi2 = kUndefined;
      break;
    }
    // Get pointer to the respective theor. data
    const TFitFunc *dataTh = fDataExp.at(i).fTh;
    // Get pointer to the exp. data
    const Data_t *dataExp = &(fDataExp.at(i).fData);
    
    // Sys. errors factor
    Double_t syserr = fDataExp.at(i).fSysErr;
    Double_t factor = p[fNparam+i];
    Double_t chi2tmp = 0.;
    if (syserr != 0.) chi2 += (chi2tmp = Chi2(1.,fDataExp.at(i).fSysErr,1.,1./factor));
    cout << "Sys errors term: chi2tmp = " << chi2tmp << " factor" << i << " = " << factor << " chi2: " << chi2 << endl;
    
    // Chi2 for each data file
    chi2 += (chi2tmp = Chi2DataSet(dataExp,dataTh,point,factor));
    cout << fDataExp.at(i).fKeyData << " chi2: " << chi2tmp << " chi2: " << chi2 << endl;
    if (chi2tmp == kUndefined)
    {
      chi2 = kUndefined;
      break;
    }
  }
  
  //
  Int_t ndf = fResult->fNpointsTot+fConstraints.size()-fMinimizer->NFree()+fNsys;
  cout << "================ Global Chi2: " << chi2 << "/(" << fResult->fNpointsTot << "+" << fConstraints.size() << "-" << fMinimizer->NFree() << "+" << fNsys << ") = " << chi2/ndf << endl;
  return chi2;
}

//__________________________________________________________________
Double_t TFit::Chi2DataSet(const Data_t *dataExp, const TFitFunc *dataTh, Point_t &point, const Double_t &factor)
{
  //
  Double_t chi2 = 0.;
  Double_t chi2tmp;
  Int_t iPoint = 1; // only to print
  for (Data_t::const_iterator iDataExp = dataExp->begin(); iDataExp != dataExp->end(); ++iDataExp)
  {
    // Calc chi2
    Int_t i;
    for (i=0; i<fNvar; ++i)
    {
      point[i] = iDataExp->at(i);
    }
    Double_t observ    = iDataExp->at(i);   // value
    Double_t observErr = iDataExp->at(i+1); // err_value
    Double_t observTh = dataTh->GetValue(point);
    chi2 += (chi2tmp = Chi2(observ,observErr,observTh,factor));
    
    // Print to std out
    cout << "exp point "<< iPoint << ":\t" << "e: "<< point[1] << "\ttheta: "<< point[0];
    cout << "\texp: "<< observ << "\tth: " << observTh << "\tchi2tmp: " << chi2tmp << "\tchi2: " << chi2 << endl;
    ++iPoint; // only for TEST
    
    //
    if (chi2tmp == kUndefined) return kUndefined;
  }
  
  //
  return chi2;
}

//__________________________________________________________________
Double_t TFit::Chi2(const Double_t &x, const Double_t &xErr, const Double_t &x0, const Double_t &f)
{
  //
  
  Double_t chi2 = kUndefined;
  if (x == kUndefined) return chi2;
  if (x0 == kUndefined) return chi2;
  if (xErr == kUndefined || xErr == 0.) return chi2;
  
//   return TMath::Power((x-f*x0)/xErr,2);
  return TMath::Power((f*x-x0)/(f*xErr),2);
  // TODO
  // Probably should be
  // TMath::Power((f*x-x0)/(f*xErr),2);
}

//__________________________________________________________________
Double_t TFit::Chi2(const Double_t &x, const Double_t &xErr, const Double_t &x0)
{
  //
  return Chi2(x,xErr,x0,1.);
}

//__________________________________________________________________
vector<TGraphErrors*> TFit::GetGraphsExp(const TString &keyObserv,
				    const Int_t ix, const Int_t iex,
				    const Int_t iy, const Int_t iey,
				    const Double_t *min, const Double_t *max, const Int_t n) const
{
  // Added to be able to plot from CINT
  
  vector<pair<Double_t,Double_t> > lim;
  for (Int_t i=0; i<n; ++i)
  {
    lim.push_back( make_pair(min[i],max[i]) );
  }
  
  return GetGraphsExp(keyObserv,ix,iex,iy,iey,lim);
}

//__________________________________________________________________
vector<TGraphErrors*> TFit::GetGraphsExp(const TString &keyObserv,
				    const Int_t ix, const Int_t iex,
				    const Int_t iy, const Int_t iey,
				    const vector<pair<Double_t,Double_t> > &lim) const
{
  // Plot all available exp data
  
  vector<TGraphErrors*> graphs;
  for (DataSetExp_t::const_iterator iData = fDataExp.begin(); iData != fDataExp.end(); ++iData)
  {
    if (iData->fKeyObserv != keyObserv) continue;
    graphs.push_back(iData->GetGraph(ix,iex,iy,iey,lim));
  }
  
  return graphs;
}

//__________________________________________________________________
TGraphErrors* TFit::GetGraphExp(const TString &keyObserv, const TString &keyData,
			    const Int_t ix, const Int_t iex,
			    const Int_t iy, const Int_t iey,
			    const vector<pair<Double_t,Double_t> > &lim) const
{
  //
  
  for (DataSetExp_t::const_iterator iData = fDataExp.begin(); iData != fDataExp.end(); ++iData)
  {
    if (iData->fKeyObserv != keyObserv) continue;
    if (iData->fKeyData != keyData) continue;
    return iData->GetGraph(ix,iex,iy,iey,lim);
  }
  
  return NULL;
}

//__________________________________________________________________
TGraphErrors* TFit::GetGraphExp(const TString &keyObserv, const TString &keyData,
			    const Int_t ix, const Int_t iex,
			    const Int_t iy, const Int_t iey,
			    const Double_t *min, const Double_t *max, const Int_t n) const
{
  //
  
  vector<pair<Double_t,Double_t> > lim;
  for (Int_t i=0; i<n; ++i)
  {
    lim.push_back( make_pair(min[i],max[i]) );
  }
  
  return GetGraphExp(keyObserv,keyData,ix,iex,iy,iey,lim);
  
}

//__________________________________________________________________
TGraph *TFit::GetGraphTh(const TString &keyObserv, const Int_t ix, const Int_t nPoints, const Double_t &xmin, const Double_t &xmax, const Point_t &point) const
{
  //
  
  map<TString,TFitFunc*>::const_iterator iData = fDataTh.find(keyObserv);
  if (iData == fDataTh.end()) return NULL;
  
  //
  return iData->second->GetGraph(ix,nPoints,xmin,xmax,point);
}

//__________________________________________________________________
TGraph *TFit::GetGraphFitFunc(const TString &name, const Int_t ix, const Int_t nPoints, const Double_t &xmin, const Double_t &xmax, const Point_t &point) const
{
  //
  
  for (vector<TFitFunc*>::iterator iFunc = fListFitFun.begin(); iFunc != fListFitFun.end(); ++iFunc)
  {
    if ((*iFunc)->GetName() == name) return (*iFunc)->GetGraph(ix,nPoints,xmin,xmax,point);
  }
  
  //
  return NULL;
}

//__________________________________________________________________
TFitFunc *TFit::GetFitFunc(const TString &name)
{
  //
  
  for (UInt_t i=0; i<fListFitFun.size(); ++i)
  {
    if (fListFitFun[i]->GetName() == name) return fListFitFun[i];
  }
  
  //
  return NULL;
}

//__________________________________________________________________
void TFit::StoreInputPar(const vector<TString> &record)
{
  // Store input parameters in a TString object
  TString s = "";
  for (UInt_t i=0; i<record.size(); ++i)
  {
//     fInputPar.SetString(fInputPar.GetString()+record.at(i)+" ");
    s += record.at(i)+" ";
  }
  fInputPar.SetTitle(fInputPar.GetTitle()+s+"\n");
}
