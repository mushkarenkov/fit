// $Id: TFitFunc.cxx 204 2013-07-19 03:33:55Z mushkar $

#include "TFitFunc.h"
using namespace Fit;

ClassImp(TFitFunc)

//______________________________________________________________________
Bool_t TFitFunc::GetAlong(const Int_t ix, const Int_t nx, const Double_t &xmin, const Double_t &xmax, Point_t p, vector<Double_t> &x, vector<Double_t> &y) const
{
  //
  
  if ( (xmax == xmin) || (nx <= 1) ) return kFALSE;
  //
  Double_t dx = (xmax-xmin)/(nx-1);
  for (Int_t i=0; i<nx; ++i)
  {
    Double_t xx = xmin + dx*i;
    p[ix] = xx;
    Double_t yy = GetValue(p);
    if (yy == kUndefined) continue;
    x.push_back(xx);
    y.push_back(yy);
  }
  //
  if (x.empty()) return kFALSE;
  return kTRUE;
}

//______________________________________________________________________
TGraph *TFitFunc::GetGraph(const Int_t ix, const Int_t nPoints, const Double_t &xmin, const Double_t &xmax, const Point_t &point) const
{
  //
  
  vector<Double_t> xx,yy;
  if (!GetAlong(ix,nPoints,xmin,xmax,point,xx,yy)) return NULL;
  //
  Int_t n = xx.size();
  Double_t *x = new Double_t[n];
  Double_t *y = new Double_t[n];
  for (Int_t i=0; i<n; ++i)
  {
    x[i] = xx.at(i);
    y[i] = yy.at(i);
  }
  TGraph *gr = new TGraph(n,x,y);
  delete x;
  delete y;
  //
  return gr;
  
}
