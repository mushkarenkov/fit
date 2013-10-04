// $Id: TFitFuncLinApprox.h 211 2013-07-24 03:25:34Z mushkar $
#ifndef __TFitFuncLinApprox_h__
#define __TFitFuncLinApprox_h__

#include "FitTypes.h"
#include "TFitFunc.h"

namespace Fit
{
  class TFitFuncLinApprox : public TFitFunc
  {
  private:
    const TFitFunc *fFitFunc;
    UInt_t fNvar;
    UInt_t fNpar;
    struct LinApprox_t
    {
      Point_t  fX;
      Double_t fFunc;
      vector<Double_t> fDeriv;
    };
    vector<LinApprox_t> fApprox;
    vector<Double_t> fPar; // best guessed param
    vector<Bool_t> fIsFixedPar; // if true than fDeriv set 0 to save time for calculation
  private:
    virtual Int_t GetIndexOfApproxAt(const Point_t&) const;
    virtual LinApprox_t CalcLinApproxAt(const Point_t&) const;
  public:
    TFitFuncLinApprox() {}
    TFitFuncLinApprox(const TString&, const TString&, const TFitFunc*);
    TFitFuncLinApprox(const TString&, const TString&);
    ~TFitFuncLinApprox() {}
    virtual Bool_t AddApproxAt(const Point_t&);
    virtual void AddPar(const Double_t &par)
    {
      fPar.push_back(par);
      fNpar = fPar.size();
      fIsFixedPar.push_back(kFALSE);
    }
    virtual void AddParFixed(const Double_t &par)
    {
      AddPar(par);
      fIsFixedPar[fNpar-1] = kTRUE;
    }
    virtual void SetPar(const Int_t i, const Double_t &par)
    {
      if (static_cast<Int_t>(fPar.size()) < i) fPar.resize(i+1);
      fPar[i] = par;
      fNpar = fPar.size();
    }
    virtual void SetParFixed(const Int_t i, const Bool_t &type)
    {
      if (static_cast<Int_t>(fIsFixedPar.size()) < i)
      {
	cout << "WARNING! >>> i > fIsFixedPar.size()" << endl;
	return;
      }
      fIsFixedPar[i] = type;
    }
    virtual void SetPar(const Int_t n, const Double_t *par)
    {
      for (Int_t i=0; i<n; ++i)	SetPar(i,par[i]);
    }
    virtual void SetNvar(const Int_t n) { fNvar = n; }
    virtual void SetFitFunc(const TFitFunc *fun) { fFitFunc = fun; }
    virtual void ReCalcLinApprox();
    virtual void ReCalcLinApproxAt(const Point_t&);
    virtual Double_t GetValue(const Point_t&) const;
    virtual void LoadLinApprox(const TString&) {/* TODO Read approx from file */}
    
    ClassDef(TFitFuncLinApprox,1) // The class describes an abstract fitting function
  };
}

#endif
