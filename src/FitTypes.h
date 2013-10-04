// SVN: $Id: FitTypes.h 207 2013-07-20 03:27:41Z mushkar $
#ifndef __FitTypes_h__
#define __FitTypes_h__

#include <vector>
#include <iostream>
#include <map>
using namespace std;

#include <TString.h>
#include <TGraphErrors.h>


namespace Fit
{
  class TFitFunc;
  
  //____________________________________________________________
  typedef vector<Double_t> Point_t;
  typedef vector<Point_t> Data_t;
  
  //____________________________________________________________
  class DataExp_t
  {
  public:
    TString   fKeyObserv;
    TString   fKeyData;
    Data_t    fData;
    TFitFunc *fTh;
    Double_t  fSysErr;
//     TGraphErrors *GetGraphCINT(const Int_t ix, const Int_t iex, const Int_t iy, const Int_t iey, const Double_t *min, const Double_t *max, const Int_t n) const
//     {
//       vector<pair<Double_t,Double_t> > lim;
//       for (Int_t i=0; i<n; ++i)
//       {
// 	lim.push_back( make_pair(min[i],max[i]) );
//       }
//       return GetGraph(ix,iex,iy,iey,lim);
//     }
    TGraphErrors *GetGraph(const Int_t ix, const Int_t iex, const Int_t iy, const Int_t iey, const vector<pair<Double_t,Double_t> > &lim) const
    {
      // Get graph to plot the data
      Int_t nPoints = 0;
      const Int_t kNmaxNpoints = 10000;
      Double_t x[kNmaxNpoints], ex[kNmaxNpoints], y[kNmaxNpoints], ey[kNmaxNpoints];
      for (Data_t::const_iterator iPoint = fData.begin(); iPoint != fData.end(); ++iPoint)
      {
	if (iPoint->size() != lim.size()) { cout << "ERROR: iPoint->size() != lim.size()" << endl; return NULL;	}
	Bool_t reject = kFALSE;
	for (UInt_t i=0; i<iPoint->size(); ++i)
	{
	  if (iPoint->at(i) < lim.at(i).first || iPoint->at(i) > lim.at(i).second) { reject = kTRUE; break; }
	}
	//
	if (reject) continue;
	// x
	x[nPoints] = iPoint->at(ix);
	if (iex < 0 || static_cast<UInt_t>(iex) >= iPoint->size())
	  ex[nPoints] = 0.;
	else
	  ex[nPoints] = iPoint->at(iex);
	// y
	y[nPoints] = iPoint->at(iy);
	if (iey < 0 || static_cast<UInt_t>(iey) >= iPoint->size())
	  ey[nPoints] = 0.;
	else
	  ey[nPoints] = iPoint->at(iey);
	++nPoints;
      }
      //
      if (!nPoints) return NULL;
      return new TGraphErrors(nPoints,x,y,ex,ey);
    }
  };
  //____________________________________________________________
  typedef vector<DataExp_t> DataSetExp_t;
  
  //____________________________________________________________
  struct Contour_t
  {
    Int_t ix;
    Int_t iy;
    vector<Double_t> x;
    vector<Double_t> y;
  };
  
  //____________________________________________________________
  typedef Double_t (*Constraint_t)(const Double_t*);
  typedef void (*Converter_t)(const Double_t*, Double_t*);
  
  //
  const Double_t kUndefined = 1.0e+10;
  
}

#endif
