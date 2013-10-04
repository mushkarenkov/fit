// $Id: TFitFuncLinApprox.cxx 211 2013-07-24 03:25:34Z mushkar $

#include "TFitFuncLinApprox.h"
using namespace Fit;

ClassImp(TFitFuncLinApprox)

//______________________________________________________________
TFitFuncLinApprox::TFitFuncLinApprox(const TString &name, const TString &title, const TFitFunc *fun)
{
  // Constructor
  fName = name;
  fTitle = title;
  fFitFunc = fun;
  fNvar = 0;
  fNpar = 0;
}

//______________________________________________________________
TFitFuncLinApprox::TFitFuncLinApprox(const TString &name, const TString &title)
{
  // Constructor
  fName = name;
  fTitle = title;
  fNvar = 0;
  fNpar = 0;
}

//______________________________________________________________
Bool_t TFitFuncLinApprox::AddApproxAt(const Point_t &p)
{
  // Add a new element to fApprox array calculated at p
  
  //
  if (p.size()<fNvar) return kFALSE;
  if (!fFitFunc) return kFALSE;
  
  // Create a point with size = fNvar
  Point_t x(fNvar);
  for (UInt_t j=0; j<fNvar; ++j) x[j] = p[j];
  
  // Don't calculate and add anything if it's already calculated at this point
  if (GetIndexOfApproxAt(x)>=0) return kFALSE;
  
  // Make a new lin approximation at p
  fApprox.push_back(CalcLinApproxAt(x));
  
  //
  return kTRUE;
}

//______________________________________________________________
Int_t TFitFuncLinApprox::GetIndexOfApproxAt(const Point_t &x) const
{
  // Check if an approximation at x already exist
  
  for (UInt_t i=0; i<fApprox.size(); ++i)
  {
    if (x == fApprox[i].fX) return i;
  }
  
  //
  return -1;
}

//______________________________________________________________
TFitFuncLinApprox::LinApprox_t TFitFuncLinApprox::CalcLinApproxAt(const Point_t &x) const
{
  // Make lin approx at x
  
  LinApprox_t approx;
  approx.fX = x;
  Point_t p(fNvar+fNpar);
  for (UInt_t i=0; i<fNvar; ++i) p[i] = approx.fX[i];
  for (UInt_t i=0; i<fNpar; ++i) p[fNvar+i] = fPar[i];
  approx.fFunc = fFitFunc->GetValue(p);
  approx.fDeriv.resize(fNpar);
  for (UInt_t i=0; i<fNpar; ++i)
  {
    if (fIsFixedPar[i])
    {
      approx.fDeriv[i] = 0.;
      continue;
    }
    Point_t p1 = p;
    Int_t ii = i+fNvar;
    p1[ii] += 0.1; // TODO Do I need a castomizable step?
    Double_t f1 = fFitFunc->GetValue(p1);
    approx.fDeriv[i] = (f1-approx.fFunc)/(p1[ii]-p[ii]);
  }
  
  //
  return approx;
}

//______________________________________________________________
void TFitFuncLinApprox::ReCalcLinApprox()
{
  // Re-calculate all lin approximations in the fApprox array
  
  //
  if (!fFitFunc) return;
  for (UInt_t i=0; i<fApprox.size(); ++i) fApprox[i] = CalcLinApproxAt(fApprox[i].fX);
}

//______________________________________________________________
void TFitFuncLinApprox::ReCalcLinApproxAt(const Point_t &p)
{
  //
  
  if (p.size()<fNvar) return;
  if (!fFitFunc) return;
  
  // Create a point with size = fNvar
  Point_t x(fNvar);
  for (UInt_t j=0; j<fNvar; ++j) x[j] = p[j];
  
  // Don't calculate and add anything if it's already calculated at this point
  Int_t iApprox = GetIndexOfApproxAt(x);
  if (iApprox<0) return;
  
  // Recalc the new lin approximation at x
  fApprox[iApprox] = CalcLinApproxAt(x);
}

//______________________________________________________________
Double_t TFitFuncLinApprox::GetValue(const Point_t &p) const
{
  // Return a linearly approximated value
  
  //
  if (p.size()<fNvar) return kUndefined;
  
  // Find corresponding approximation in fApprox
  Point_t x(fNvar);
  for (UInt_t i=0; i<fNvar; ++i) x[i] = p[i];
  Int_t iApprox = GetIndexOfApproxAt(x);
  if (iApprox<0) return kUndefined;
  
  // Calc the linearly approximated value
  Double_t f = fApprox[iApprox].fFunc;
  for (UInt_t i=0; i<fNpar; ++i) f += fApprox[iApprox].fDeriv[i]*(p[i+fNvar]-fPar[i]);
    
  //
  return f;
  
}
