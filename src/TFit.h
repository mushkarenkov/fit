// SVN: $Id: TFit.h 209 2013-07-23 04:41:01Z mushkar $
#ifndef __TFit_h__
#define __TFit_h__

#include <vector>
#include <set>
#include <map>
using namespace std;

#include <TMath.h>
#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <Math/Functor.h>
#include <TNamed.h>
#include <TString.h>
#include <TCanvas.h>
// #include <TRint.h>
#include <TFile.h>
class TGraphErrors;

#include "FitTypes.h"
#include "TFitFunc.h"

namespace Fit
{
  //
  typedef map<TString,TFitFunc*> DataSetTh_t;
  struct Par_t
  {
    TString  fName;
    Bool_t   fFix;
    Double_t fInit;
    Double_t fStep;
    Double_t fMin;
    Double_t fMax;
    Double_t fValue;
    Double_t fError;
  };
  
  
//_____________________________________________________________________________
  class TFitResult : public TNamed
  {
  public:
    TFitResult(const TString &name, const TString &title) : TNamed(name,title) {}
    vector<Par_t> fParams;
    Int_t fNparam;
    Int_t fNparamFree;
    Int_t fNconstraints;
    Int_t fNpointsTot;
//     map<TString,Double_t> fChi2Constraints;
    Double_t fChi2;
    set<Double_t> fChi2Point;
    set<Double_t> fChi2Set;
    vector<Contour_t> fContours;
    //
    ClassDef(TFitResult,1) // The class describes an object that contains a fit result
  };
  
  // TFit
//   class TFit : public TRint
  
//_____________________________________________________________________________  
  class TFit
  {
  private:
    ROOT::Math::Minimizer *fMinimizer;
    vector<Par_t> fParams;
    Int_t fNparam; // # of parameters to fit
    Int_t fNvar; // # of variables
    Int_t fNsys;
    DataSetExp_t fDataExp; // List of experimental data
    DataSetTh_t  fDataTh; // List of theoretical data
    vector<Constraint_t> fConstraints;
    //
    void (TFit::*fMinimize)(const Char_t*, const Char_t*);
    //
    Converter_t fConvertParams;
    //
    static vector<TFitFunc*> fListFitFun;
    static map<TString,Constraint_t> fListConstraints;
    static map<TString,Converter_t> fListConverters;
    // Result
//     TString fNameFileResult;
    TFitResult *fResult;
    TNamed fInputPar; // input parameters to store in the output file
    //
    virtual void InitDataMembers();
    virtual void StoreInputPar(const vector<TString>&);
    virtual void SaveResult();
  public:
//     TFit(const char*, Int_t*, char**, const TString&);
    TFit(const TString&);
    TFit();
    virtual ~TFit();
    virtual TFitResult *GetFitResult() const { return fResult; }
    const virtual DataSetTh_t *GetDataTh() const { return &fDataTh; }
    virtual TFitFunc *GetFitFunc(const TString&);
    const virtual DataSetExp_t *GetDataExp() const { return &fDataExp; }
    const virtual DataSetExp_t::const_iterator GetDataExp(const TString &keyObserv, const TString &keyData) const
    {
      //
      DataSetExp_t::const_iterator iData;
      for (iData = fDataExp.begin(); iData != fDataExp.end(); ++iData)
      {
	if (iData->fKeyObserv != keyObserv) continue;
	if (iData->fKeyData != keyData) continue;
	break;
      }
      return iData;
    }
    virtual void ReadParams(const TString&);
    virtual void ReadDataExp(const vector<TString>&);
    virtual void ReadDataTh(const vector<TString>&);
    static Bool_t ReadDataFile(const TString&, Data_t&);
    virtual void MinimizeNormally(const Char_t *minName = "Minuit2", const Char_t *algoName = "");
    virtual void MinimizeLinApprox(const Char_t *minName = "Minuit2", const Char_t *algoName = "");
    virtual void SetMinimizer(const Char_t *minName = "Minuit2", const Char_t *algoName = "");
    virtual void Minimize(const Char_t *minName = "Minuit2", const Char_t *algoName = "") { (this->*fMinimize)(minName,algoName); }
    virtual Bool_t IsStableMinimum(const vector<vector<Double_t> >&, const Double_t&) const;
    virtual vector<TGraphErrors*> GetGraphsExp(const TString&, const Int_t, const Int_t, const Int_t, const Int_t, const vector<pair<Double_t,Double_t> >&) const;
    virtual vector<TGraphErrors*> GetGraphsExp(const TString&, const Int_t, const Int_t, const Int_t, const Int_t, const Double_t*, const Double_t*, const Int_t) const;
    virtual TGraphErrors *GetGraphExp(const TString&, const TString&, const Int_t, const Int_t, const Int_t, const Int_t, const vector<pair<Double_t,Double_t> >&) const;
    virtual TGraphErrors *GetGraphExp(const TString&, const TString&, const Int_t, const Int_t, const Int_t, const Int_t, const Double_t*, const Double_t*, const Int_t) const;
    virtual TGraph *GetGraphTh(const TString&, const Int_t, const Int_t, const Double_t&, const Double_t&, const Point_t&) const;
    virtual TGraph *GetGraphFitFunc(const TString&, const Int_t, const Int_t, const Double_t&, const Double_t&, const Point_t&) const;
    //
    static void AddFitFunc(TFitFunc *fun) { fListFitFun.push_back(fun); }
    static void AddConstraint(const TString &name, const Constraint_t &c) { fListConstraints[name] = c; }
    static void AddConverter(const TString &name, const Converter_t &c) { fListConverters[name] = c; }
    //
    static Double_t Chi2(const Double_t&, const Double_t&, const Double_t&);
    static Double_t Chi2(const Double_t&, const Double_t&, const Double_t&, const Double_t&);
    virtual Double_t Chi2DataSet(const Data_t*, const TFitFunc*, Point_t&, const Double_t&);
    virtual Double_t GlobalChi2(const Double_t*);
    //
    virtual void ConvertParams(const Double_t*, Double_t*);
    //
    TString fNameFileResult;
    
    ClassDef(TFit,1) // The class describes an object that does a fit
  };
}

#endif
