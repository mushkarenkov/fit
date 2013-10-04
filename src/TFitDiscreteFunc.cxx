// SVN: $Id: TFitDiscreteFunc.cxx 94 2013-04-24 23:49:11Z mushkar $

#include <Math/Interpolator.h>

#include "TFitDiscreteFunc.h"
using namespace Fit;

ClassImp(TFitDiscreteFunc)

//___________________________________________________________________
map<Point_t,Int_t>::const_iterator TFitDiscreteFunc::FindNextCoord(const Point_t &point, const Int_t iX)
{
  // Find the next point greater then "point" wrt iX coordinate
  map<Point_t,Int_t>::const_iterator iCoordNext = fCoord.end(); // default
  Point_t pointNext;
  for (set<Double_t>::const_iterator ix = fX.at(iX).upper_bound(point.at(iX)); (ix != fX.at(iX).end() && iCoordNext == fCoord.end()); ++ix)
  {
    pointNext = point;
    pointNext[iX] = *ix;
    iCoordNext = fCoord.find(pointNext);
  }
  return iCoordNext;
}

//___________________________________________________________________
Bool_t TFitDiscreteFunc::GetAlong(const Int_t i, Point_t p, vector<Double_t> &x, vector<Double_t> &y) const
{
  // returns tha arrays y and x (y=y(x)) along the ix axis
  
  for (set<Double_t>::const_iterator ix = fX.at(i).begin(); ix != fX.at(i).end(); ++ix)
  {
    p[i] = *ix;
    map<Point_t,Int_t>::const_iterator iy = fCoord.find(p);
    if (iy == fCoord.end()) continue;
    x.push_back(p.at(i));
    y.push_back(fValue.at(iy->second));
  }
  
  //
  if (x.empty()) return kFALSE;
  //
  return kTRUE;
  
}

//___________________________________________________________________
Bool_t TFitDiscreteFunc::BelongToDomain(const Point_t &p) const
{
  // Checks if p belongs to the function domain
  
  for (UInt_t i=0; i<p.size(); ++i)
  {
    if ( p.at(i) < *(fX.at(i).begin()) || p.at(i) > *(fX.at(i).rbegin()) ) return kFALSE;
  }
  
  //
  return kTRUE;
}

//___________________________________________________________________
Double_t TFitDiscreteFunc::Interpolate2(Point_t point) const
{
  // Interpolate N-dim grid function
  
  // empty point
  if (point.empty()) return kUndefined;
  
  // if # of vars in the point != # of vars in the function -> return kUndefined
  if (point.size() != fX.size()) return kUndefined;
  
  if (fX.empty())
  {
    // In case of a constant function -> return the constant
    if (!fValue.empty()) return fValue.front();
    else return kUndefined;
  }
  
  //
  if (!BelongToDomain(point)) return kUndefined;
  
  // Get a coordinate point in the function (can be taken any point)
  Point_t p = fCoord.begin()->first;
  
  // Interpolation algorithm
//   static const ROOT::Math::Interpolation::Type kTypeInterp =  ROOT::Math::Interpolation::kLINEAR;
  static const ROOT::Math::Interpolation::Type kTypeInterp =  ROOT::Math::Interpolation::kCSPLINE;
//   static const ROOT::Math::Interpolation::Type kTypeInterp =  ROOT::Math::Interpolation::kPOLYNOMIAL;
//   static const ROOT::Math::Interpolation::Type kTypeInterp =  ROOT::Math::Interpolation::kCSPLINE_PERIODIC;
//   static const ROOT::Math::Interpolation::Type kTypeInterp = ROOT::Math::Interpolation::kAKIMA;
  
  // 1-D case
  if (fX.size() == 1)
  {
    vector<Double_t> x, y;
    GetAlong(0,p,x,y);
    ROOT::Math::Interpolator interp(x, y, kTypeInterp);
    // Final return (recursion should come here)
    return interp.Eval(point.at(0));
  }
  
  // N-Dim case
  Data_t data;
  for (map<Point_t,Int_t>::const_iterator iCoord = fCoord.begin(); iCoord != fCoord.end(); ++iCoord)
  {
    if (iCoord->first.back() != p.back()) continue;
    Point_t pNew = iCoord->first;
    vector<Double_t> x, y;
    GetAlong(fX.size()-1,pNew,x,y);
    ROOT::Math::Interpolator interp(x, y, kTypeInterp);
    pNew[fX.size()-1] = interp.Eval(point.back());
    data.push_back(pNew);
  }
  
  //
  point.pop_back();
  return TFitDiscreteFunc(data).Interpolate2(point);
}


//___________________________________________________________________
Double_t TFitDiscreteFunc::Interpolate(const Point_t &point) const
{
  // Interpolate
  
  // empty point
  if (point.empty()) return kUndefined;
  
  if (fX.empty())
  {
    // In case of a constant function -> return the constant
    if (!fValue.empty()) return fValue.front();
    else return kUndefined;
  }
  
  // if var number in the point less then in the function -> return kUndefined
  if (point.size() < fX.size()) return kUndefined;
  
  //
  // *----*----*---*
  // |    | +  |   |
  // *----0----*---*
  // |    |    |   |
  // *----*----*---*
  //
  // Find the point "0", corresponding to the point "+"
  Point_t point0;
  set<Double_t>::const_iterator ix;
  for (UInt_t i=0; i<fX.size(); ++i)
  {
    ix = fX.at(i).upper_bound(point.at(i));
    if (ix == fX.at(i).begin()) return kUndefined;
    point0.push_back(*(--ix));
    cout << "x_"<<i<<": " << point.at(i) << "\tx0_"<<i<<": " << *(point0.end()-1) << endl;
  }
  
  // if the point0 doesn't exist -> return kUndefined // TODO сделать возможность учитывать пропуски (дырки) в сетке
  map<Point_t,Int_t>::const_iterator iPoint0 = fCoord.find(point0);
  if (iPoint0 == fCoord.end()) return kUndefined;
  Int_t indexPoint0 = iPoint0->second;
  cout << "indexPoint0: " << indexPoint0 << endl;
  
  // Make interpolation
  // if at least 1 derivative doesn't exist -> return kUndefined
  Double_t valueInterpolated = fValue.at(indexPoint0);
  Double_t deriv;
  for (UInt_t i=0; i<fX.size(); ++i)
  {
    deriv = fDeriv.at(indexPoint0).at(i);
    cout << "deriv"<<i<<": " << deriv << endl;
    if (deriv == kUndefined) return kUndefined;
    valueInterpolated += deriv*(point.at(i)-iPoint0->first.at(i));
    cout << "valueInterpolated"<<i<<": " << valueInterpolated << endl;
  }
  
  //
  return valueInterpolated;
}

//___________________________________________________________________
void TFitDiscreteFunc::SetData2(const Data_t &data)
{
  // NOTE: Temporary solution for Sigma_2x
  if (data.empty()) return;
  
  {
    // Fill fX container with an empty sets<Double_t> (assume that all rows from "data" have the same # of columns)
    set<Double_t> setX;
    fX.assign(data.at(0).size()-3,setX);
    // Fill fX, fCoord and fValue
    Point_t coord;
    fNpoints = 0;
    for (Data_t::const_iterator iData = data.begin(); iData != data.end(); ++iData)
    {
      Double_t sigmap = *(iData->end()-2);
      Double_t sigmam = *(iData->end()-1);
//       cout << sigmap + sigmam << endl;
      Double_t sigma2x = (sigmap - sigmam) / (sigmap + sigmam);
      fValue.push_back(sigma2x);
      coord.clear();
      for (UInt_t i=0; i<fX.size(); ++i)
      {
	fX.at(i).insert(iData->at(i));
	coord.push_back(iData->at(i));
      }
      fCoord.insert( make_pair(coord,fNpoints++) );
    }
  }
  
//   // fDeriv
//   map<Point_t,Int_t>::const_iterator iCoord[2];
//   Double_t x[2], value[2], dVdX;
//   vector<Double_t> deriv;
//   fDeriv.assign(fCoord.size(),deriv); // make array with fCoord.size() elements
//   for (iCoord[0] = fCoord.begin(); iCoord[0] != fCoord.end(); ++iCoord[0])
//   {
//     deriv.clear();
//     value[0] = fValue.at(iCoord[0]->second);
//     for (UInt_t i=0; i<fX.size(); ++i)
//     {
//       x[0] = iCoord[0]->first.at(i);
//       iCoord[1] = FindNextCoord(iCoord[0]->first,i);
//       dVdX = kUndefined;
//       if (iCoord[1] != fCoord.end())
//       {
// 	x[1] = iCoord[1]->first.at(i);
// 	value[1] = fValue.at(iCoord[1]->second);
// 	dVdX = (value[1]-value[0])/(x[1]-x[0]);
// 	// TEST
// // 	cout << iCoord[0]->second << "\t" << i << "\t" << value[1] << "-" << value[0] << "/" << x[1] << "-" << x[0] << " = " << dVdX << endl;	  
//       }
//       deriv.push_back(dVdX);
//     }
//     fDeriv[iCoord[0]->second] = deriv;
//   }
//   cout << "!!!!!" << endl;
  
}

//____________________________________________________________________________________
void TFitDiscreteFunc::SetData(const Data_t &data)
{
  //
  if (data.empty()) return;
  
  {
    // Fill fX container with an empty sets<Double_t> (assume that all rows from "data" have the same # of columns)
    set<Double_t> setX;
    fX.assign(data.at(0).size()-1,setX);
    // Fill fX, fCoord and fValue
    Point_t coord;
    fNpoints = 0;
    for (Data_t::const_iterator iData = data.begin(); iData != data.end(); ++iData)
    {
      fValue.push_back(*(iData->end()-1));
      coord.clear();
      for (UInt_t i=0; i<fX.size(); ++i)
      {
	fX.at(i).insert(iData->at(i));
	coord.push_back(iData->at(i));
      }
      fCoord.insert( make_pair(coord,fNpoints++) );
    }
  }
  
//   // fDeriv
//   map<Point_t,Int_t>::const_iterator iCoord[2];
//   Double_t x[2], value[2], dVdX;
//   vector<Double_t> deriv;
//   fDeriv.assign(fCoord.size(),deriv); // make array with fCoord.size() elements
//   for (iCoord[0] = fCoord.begin(); iCoord[0] != fCoord.end(); ++iCoord[0])
//   {
//     deriv.clear();
//     value[0] = fValue.at(iCoord[0]->second);
//     for (UInt_t i=0; i<fX.size(); ++i)
//     {
//       x[0] = iCoord[0]->first.at(i);
//       iCoord[1] = FindNextCoord(iCoord[0]->first,i);
//       dVdX = kUndefined;
//       if (iCoord[1] != fCoord.end())
//       {
// 	x[1] = iCoord[1]->first.at(i);
// 	value[1] = fValue.at(iCoord[1]->second);
// 	dVdX = (value[1]-value[0])/(x[1]-x[0]);
// 	// TEST
// // 	cout << iCoord[0]->second << "\t" << i << "\t" << value[1] << "-" << value[0] << "/" << x[1] << "-" << x[0] << " = " << dVdX << endl;	  
//       }
//       deriv.push_back(dVdX);
//     }
//     fDeriv[iCoord[0]->second] = deriv;
//   }
}
