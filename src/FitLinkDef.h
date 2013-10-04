// SVN: $Id: FitLinkDef.h 194 2013-07-16 08:17:34Z mushkar $
#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ namespace Fit;
#pragma link C++ class Fit::TFit+;
#pragma link C++ class Fit::TFitFunc+;
#pragma link C++ class Fit::TFitFuncLinApprox+;
#pragma link C++ class Fit::TFitDiscreteFunc+;
#pragma link C++ class Fit::DataExp_t+;
#pragma link C++ class Fit::TFitResult+;
#pragma link C++ class Fit::Contour_t+;
#pragma link C++ class vector<Fit::Contour_t>+;
#pragma link C++ class Fit::Par_t+;
#pragma link C++ class vector<Fit::Par_t>+;
#pragma link C++ class vector<TGraphErrors*>+;

#endif
