// $Id: TFitFunc.h 194 2013-07-16 08:17:34Z mushkar $
#ifndef __TFitFunc_h__
#define __TFitFunc_h__

#include <TNamed.h>

#include "FitTypes.h"

namespace Fit
{
  class TFitFunc : public TNamed
  {
  protected:
    Double_t (*fGetValue)(const Point_t&);
  public:
    TFitFunc() {}
    TFitFunc(const TString &name, const TString &title, Double_t (*fun)(const Fit::Point_t&)) : TNamed(name,title)
    {
      fGetValue = fun;
    }
    ~TFitFunc() {}
    virtual Bool_t GetAlong(const Int_t, const Int_t, const Double_t&, const Double_t&, Point_t, vector<Double_t>&, vector<Double_t>&) const;
    virtual Double_t GetValue(const Point_t &p) const
    {
      return fGetValue(p);
    }
    virtual TGraph *GetGraph(const Int_t, const Int_t, const Double_t&, const Double_t&, const Point_t&) const;
    
    ClassDef(TFitFunc,1) // The class describes an abstract fitting function
  };
}

#endif
