// SVN: $Id: TFitDiscreteFunc.h 94 2013-04-24 23:49:11Z mushkar $
#ifndef __TFitDiscreteFunc_h__
#define __TFitDiscreteFunc_h__

#include <set>
#include <iostream>
#include <map>
using namespace std;

#include "TFitFunc.h"
#include "FitTypes.h"

namespace Fit
{
  class TFitDiscreteFunc : public TFitFunc
  {
  private:
    Int_t fNpoints;
    virtual map<Point_t,Int_t>::const_iterator FindNextCoord(const Point_t&, const Int_t);
  public:
    vector<set<Double_t> > fX; // [Nvar][Npoint]
    map<Point_t,Int_t> fCoord; // <point,point_index>
    vector<Double_t> fValue;
    vector<vector<Double_t> > fDeriv; // [i_point][Nvar]
    //
    TFitDiscreteFunc() {}
    TFitDiscreteFunc(const Data_t &data)
    {
      SetData(data);
    }
    ~TFitDiscreteFunc() {}
    virtual Bool_t BelongToDomain(const Point_t&) const;
    virtual Bool_t GetAlong(const Int_t, Point_t, vector<Double_t>&, vector<Double_t>&) const;
    virtual Double_t Interpolate(const Point_t&) const;
    virtual Double_t Interpolate2(Point_t) const;
    virtual Double_t GetValue(const Point_t &point) const {return Interpolate2(point); }
    virtual void     SetData(const Data_t&);
    virtual void     SetData2(const Data_t &data);
    
    ClassDef(TFitDiscreteFunc,1) // The class describes a discrete function of n variables
  };
}

#endif
