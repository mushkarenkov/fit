// SVN: $Id: FitMain.cxx 234 2013-08-22 11:10:07Z mushkar $

// STD
#include <getopt.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
using namespace std;
// OpenMP
// #include <omp.h>
// ROOT
#include <TString.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TFrame.h>
#include <TGraphErrors.h>
#include <TLatex.h>
#include <TApplication.h>
// Fit
#include "TFit.h"
#include "FitTypes.h"
#include "TFitFunc.h"
#include "FunComptonTheory.h"

void Usage(const char*);
void ParseCommandLineArguments(int, char**,TString&);
void PrintPartialDerivatives(Fit::TFit&,const TString&,const TString&,const TString&,Fit::Point_t&,Double_t*);

void test(Fit::TFit&);
void test2(Fit::TFit&);
void test3(Fit::TFit&,TString,TString,Double_t,TString,Double_t);
void test33(Fit::TFit&);
void test4(Fit::TFit&);
void test5(Fit::TFit&);
void test6(Fit::TFit&);
void test7(Fit::TFit&);
void test8();
void test9(Fit::TFit&);

static const Double_t kInfiniti = 10e+10;

//_____________________________________________
int main(int argc, char** argv)
{
  // Parse the command line arguments
  TString nameFileParam;
  ParseCommandLineArguments(argc,argv,nameFileParam);
  
  // Create the fitting functions and constraints from the compton theory
  Fit::InitComptonTheory();
  
  // Create a fitter
  Fit::TFit fit(nameFileParam);
  
  // Do minimization
  fit.Minimize();
  
  //
//   test9(fit);
  
//   Double_t p1[6], p2[6];
//   p1[0] = 13.8;
//   p1[1] = 7.6;
//   p1[2] = -3.0;
//   p1[3] = 3.5;
//   p1[4] = -1.01;
//   p1[5] = 11.6;
//   Fit::ConvertPolariz1(p1,p2);
//   for (int i=0;i<6;++i)
//   {
//     cout << p2[i] << endl;
//   }
  
  //
  return 0;
}

//_____________________________________________
void Usage(const char *nameProgram)
{
  // Print help
  cout << endl;
  cout << "Usage: " << nameProgram << " --param=file [--help]" << endl;
  cout << endl;
  cout << "Options: " << endl;
  cout << "\t-h  --help  \t print this help and exit" << endl;
  cout << "\t-p  --param \t a file with parameters" << endl;
  cout << endl;
  
}

//_____________________________________________
void ParseCommandLineArguments(int argc, char** argv, TString &nameFileParam)
{
  // Parse the command line arguments
  
  // When no option is provided => exit
  if (argc == 1)
  {
    cout << ">>> ERROR! No option provided!" << endl;
    Usage(argv[0]);
    exit(-1);
  }
  
  // Define options
  const char *optsShort = "hp:";
  const struct option optsLong[] = {
    {"help",  no_argument      ,NULL,'h'},
    {"param", required_argument,NULL,'p'},
    {NULL,   0                 ,NULL, 0 }
  };
  
  // Parse options
  Int_t rez, iOpt = -1;
  while ( (rez=getopt_long(argc,argv,optsShort,optsLong,&iOpt)) != -1 )
  {
    switch (rez)
    {
      case 'h':
	Usage(argv[0]);
	exit(0);
	break;
      case 'p':
	nameFileParam = optarg;
	break;
      case '?':
      default:
	exit(-1);
    }
  }
  
}

//_____________________________________________
void test(Fit::TFit &fit)
{
  //
  TCanvas *canv = new TCanvas();
//   TH1 *h = canv->DrawFrame(0.,0.,500.,400.);
  TH1 *h = canv->DrawFrame(200.,-.5,350.,.8);
  canv->GetFrame()->SetFillColor(0);
  canv->GetFrame()->SetBorderSize(0);
  h->SetXTitle("titleX");
  h->SetYTitle("titleY");
  //
  vector<pair<Double_t,Double_t> > lim(4,make_pair(-kInfiniti,kInfiniti));
  lim[0] = make_pair(100.,180.);
  vector<TGraphErrors*> gr = fit.GetGraphsExp("Sigma3",1,-1,2,3,lim);
  for (UInt_t i=0; i<gr.size(); ++i)
  {
    gr.at(i)->Draw("P");
    gr.at(i)->SetMarkerStyle(20);
  }
  canv->WaitPrimitive();
  //
//   lim[0] = make_pair(80.,90.);
//   TGraphErrors *gr2 = fit.PlotExp("sigma","wol01",1,-1,2,3,lim);
//   gr2->Draw("P");
//   gr2->SetMarkerStyle(20);
//   gr2->SetMarkerColor(kRed);
//   canv->WaitPrimitive();
  
//   test2(fit,canv);
  
//   delete gr;
  delete canv;
}

//_____________________________________________
void test2(Fit::TFit &fit)
{
  //
  TCanvas *canv = new TCanvas();
//   TH1 *h = canv->DrawFrame(190.,0.,410.,600.);
//   TH1 *h = canv->DrawFrame(11.,0.,17.,600.);
  TH1 *h = canv->DrawFrame(0.,-1.,180.,1.);
//   TH1 *h = canv->DrawFrame(200.,-1.,350.,1.);
  canv->GetFrame()->SetFillColor(0);
  canv->GetFrame()->SetBorderSize(0);
  h->SetXTitle("titleX");
  h->SetYTitle("titleY");
  //
  Fit::Point_t point(8);
//   point[0] = 120.;
//   point[1] = 200.;//12.;//13.;		// 12 - 16, 2
//   point[2] = 10.;//11.;//9.;		// 8 - 12, 2
//   point[3] = 2.;//-3.;//-3.3;		// -6.3 - -2.3, 2
//   point[4] = -4.3;//-1.;//0.99;	// -2.01 - 1.99, 2
//   point[5] = -0.01;//2.;//-0.9;		// 0.1 - 4.1, 2
//   point[6] = 2.1;//3.;//-0.1;		// 0.9 - 4.9, 2
//   point[7] = 0.9;			// 5
  
  point[0] = 90.;
  point[1] = 110.;
  point[2]= 11.2461;
  point[3]= 3.87054;
  point[4]= -3.34375;
  point[5]= 2.95934;
  point[6]= 0.185504;
  point[7]= 1.14328;
  
  
//   Fit::Point_t point(4);
//   point[0] = 134.;
//   point[1] = 4.;
//   point[2] = 2.;
//   point[3] = 103.;
  
//   TGraph* gr = fit.PlotTh("Sigma2x",0,20,point);
  TGraph* gr = fit.GetGraphFitFunc("PascalutsaSigma3",0,25,5.,175.,point);
  if (gr)
  {
    gr->Draw("L");
    gr->SetLineColor(kRed);
    gr->SetLineWidth(3);
    canv->WaitPrimitive();
  }
  else
  {
    cout << "Nothing to plot!" << endl;
  }
  //
  TGraph* gr2 = fit.GetGraphFitFunc("PasquiniSigma3",0,5,0.,180.,point);
  if (gr2)
  {
    gr2->Draw("Lsame");
    gr2->SetLineColor(kGreen);
    gr2->SetLineWidth(3);
    canv->WaitPrimitive();
  }
  else
  {
    cout << "Nothing to plot!" << endl;
  }
  //
  delete gr;
  delete gr2;
  delete canv;
}

//_____________________________________________
void test33(Fit::TFit &fit)
{
  test3(fit,"40","50",45.,"45",375.);
  test3(fit,"55","65",60.,"60",350.);
  test3(fit,"85","95",90.,"90",250.);
  test3(fit,"130","140",135.,"135",200.);
  test3(fit,"150","160",155.,"155",50.);
  test3(fit,"170","180",175.,"175",40.);
  
//   test3(fit,"50","60",55.,"55",375.);
//   test3(fit,"70","80",75.,"75",350.);
//   test3(fit,"90","100",95.,"95",250.);
//   test3(fit,"110","120",115.,"115",200.);
//   test3(fit,"150","160",155.,"155",50.);
//   test3(fit,"170","180",175.,"175",40.);
}

//_____________________________________________
void test3(Fit::TFit &fit,TString min, TString max, Double_t avg, TString savg,Double_t s)
{
  //
//   TFile *fileResult = new TFile("test.root");
//   Fit::TFitResult *result = (Fit::TFitResult*)fileResult->Get("Result");
//   vector<Double_t> par = result->fPar;
//   cout << par.at(0) << endl;
  
  //
  TCanvas *canv = new TCanvas();
  TH1 *h = canv->DrawFrame(0.,0.,500.,s);
//   TH1 *h = canv->DrawFrame(200.,-.5,350.,.8);
  canv->GetFrame()->SetFillColor(0);
  canv->GetFrame()->SetBorderSize(0);
  h->SetXTitle("E_{#gamma} (MeV)");
  h->SetYTitle("d#sigma/d#Omega (nb/sr)");
//   h->SetYTitle("#Sigma_{3}");
  TLatex text;
  text.SetNDC();
//   text.DrawLatex(0.4,0.93,"#theta_{#gamma} = 40#circ - 50#circ");
  text.DrawLatex(0.4,0.93,"#theta_{#gamma} = "+min+"#circ - "+max+"#circ");
//   text.DrawLatex(0.4,0.93,"#theta_{#gamma} = 85#circ - 95#circ");
//   text.DrawLatex(0.4,0.93,"#theta_{#gamma} = 130#circ - 140#circ");
//   text.DrawLatex(0.4,0.93,"#theta_{#gamma} = 150#circ - 160#circ");
//   text.DrawLatex(0.4,0.93,"#theta_{#gamma} = 150#circ - 160#circ");
//   text.DrawLatex(0.4,0.93,"#theta_{#gamma} = 170#circ - 180#circ");
  // exp
  vector<pair<Double_t,Double_t> > lim(4,make_pair(-kInfiniti,kInfiniti));
  lim[0] = make_pair(avg-5.,avg+5.);
  vector<TGraphErrors*> gr = fit.GetGraphsExp("sigma",1,-1,2,3,lim);
  for (UInt_t i=0; i<gr.size(); ++i)
  {
    gr.at(i)->Draw("P");
    gr.at(i)->SetMarkerStyle(20);
  }
//   canv->WaitPrimitive();
  // theor
  Fit::Point_t point(4);
  point[0] = 134.;
  point[1] = 11.0232;
  point[2] = 3.83113;
  point[3] = avg;
  TGraph* grTh = fit.GetGraphTh("sigma",0,10,0.,100.,point);
  point[1] = 12.;
  point[2] = 1.9;
  TGraph* grThPDG = fit.GetGraphTh("sigma",0,10,0.,100.,point);
  if (grTh)
  {
    grTh->Draw("L");
    grThPDG->Draw("L");
    grTh->SetLineColor(kRed);
    grThPDG->SetLineColor(kBlue);
    grTh->SetLineWidth(3);
    grThPDG->SetLineWidth(3);
    canv->WaitPrimitive();
    canv->SaveAs("Lvov_UnpolCS_LEGS_Baldin"+savg+".pdf");
  }
  else
  {
    cout << "Nothing to plot!" << endl;
  }
    
  delete canv;
}

//_____________________________________________
void test4(Fit::TFit &fit)
{
  Fit::Point_t point(8);
  point[0] = 100.;
  point[1] = 240.;
  point[2] = 11.;
  point[3] = 3.5;
  point[4] = -4.3;
  point[5] = -1.01;
  point[6] = 1.1;
  point[7] = 1.9;
  
  cout << fit.GetDataTh()->find("Sigma2x")->second->GetValue(point) << endl;
}

//_____________________________________________
void test5(Fit::TFit &fit)
{
  //
  
  TCanvas *canv = new TCanvas();
  TH1 *h = canv->DrawFrame(0.,-0.7,180.,0.7);
//   TH1 *h = canv->DrawFrame(200.,-.5,350.,.8);
  canv->GetFrame()->SetFillColor(0);
  canv->GetFrame()->SetBorderSize(0);
  h->SetXTitle("#theta (degrees)");
  h->SetYTitle("#Sigma_{2x}");
//   h->SetYTitle("#Sigma_{3}");
//   TLatex text;
//   text.SetNDC();
//   text.DrawLatex(0.4,0.93,"#theta_{#gamma} = 40#circ - 50#circ");
//   text.DrawLatex(0.4,0.93,"#theta_{#gamma} = "+min+"#circ - "+max+"#circ");
//   text.DrawLatex(0.4,0.93,"#theta_{#gamma} = 85#circ - 95#circ");
//   text.DrawLatex(0.4,0.93,"#theta_{#gamma} = 130#circ - 140#circ");
//   text.DrawLatex(0.4,0.93,"#theta_{#gamma} = 150#circ - 160#circ");
//   text.DrawLatex(0.4,0.93,"#theta_{#gamma} = 150#circ - 160#circ");
//   text.DrawLatex(0.4,0.93,"#theta_{#gamma} = 170#circ - 180#circ");
  // exp
  vector<pair<Double_t,Double_t> > lim(4,make_pair(-kInfiniti,kInfiniti));
//   Double_t avg = 286.5;
//   lim[1] = make_pair(avg-1.,avg+1.);
  vector<TGraphErrors*> gr = fit.GetGraphsExp("Sigma2x",0,-1,2,3,lim);
  for (UInt_t i=0; i<gr.size(); ++i)
  {
    gr.at(i)->Draw("P");
    gr.at(i)->SetMarkerStyle(20);
  }
//   canv->WaitPrimitive();
  // theor
  Fit::Point_t point(8);
  // Pasquini
//   point[0] = 134.;
//   point[1] = 285.;
//   point[2] = 11.5544+3.55136;
//   point[3] = 11.5544-3.55136;
//   point[4] = -4.31919;
//   point[5] = 0.849831;
//   point[6] = 3.57917;
//   point[7] = 0.900146;
  //
  //Pascalutsa
//   point[0] = 134.;
//   point[1] = 286.5;
//   point[2] = 12.1085;//8.81026;
//   point[3] = 5.;//4.9792;
//   point[4] = -1.8471;//-3.99551;
//   point[5] = 5.0429e-09;//0.144839;
//   point[6] = -2.;//1.99987;
//   point[7] = 4.86411;//2.849;
  
  point[0] = 134.;
  point[1] = 285.;
  point[2]=11.2461;
  point[3]=3.87054;
  point[4]=-3.34375;
  point[5]=2.95934;
  point[6]=0.185504;
  point[7]=1.14328;
  

//   TGraph* grThBaldin = fit.PlotTh("Sigma2x",0,25,point);
//   point[0] = 134.;
//   point[1] = 285.;
//   point[2] = 11.167;
//   point[3] = 3.86921;
//   point[4] = -4.;
//   point[5] = 0.0199647;
//   point[6] = 2.;
//   point[7] = 3.;
//   TGraph* grThChPT = fit.PlotTh("Sigma2x",0,25,point);
  TGraph* grThChPT = fit.GetGraphFitFunc("PascalutsaSigma2x",0,25,5.,175.,point);
  TGraph* grThDR = fit.GetGraphFitFunc("PasquiniSigma2x",0,10,5.,175.,point);
//   if (grTh)
//   {
    grThDR->Draw("L");
    grThChPT->Draw("L");
    grThDR->SetLineColor(kRed);
    grThChPT->SetLineColor(kBlue);
    grThDR->SetLineWidth(3);
    grThChPT->SetLineWidth(3);
    canv->WaitPrimitive();
    canv->SaveAs("Sigma2x_Phil_Pascalutsa-Pasquini.pdf");
//   }
//   else
//   {
//     cout << "Nothing to plot!" << endl;
//   }
    
//   delete canv;
}

//_____________________________________________
void test6(Fit::TFit &fit)
{
  //
    TCanvas *canv = new TCanvas();
  TH1 *h = canv->DrawFrame(0.,-0.7,180.,0.7);
//   TH1 *h = canv->DrawFrame(200.,-.5,350.,.8);
  canv->GetFrame()->SetFillColor(0);
  canv->GetFrame()->SetBorderSize(0);
  h->SetXTitle("#theta (degrees)");
//   h->SetYTitle("#Sigma_{2x}");
  h->SetYTitle("#Sigma_{3}");
  // exp
  vector<pair<Double_t,Double_t> > lim(4,make_pair(-kInfiniti,kInfiniti));
//   lim[1] = make_pair(333.8,333.8);
  lim[1] = make_pair(286.5,286.5);
//   vector<TGraphErrors*> gr = fit.PlotExp("Sigma2x",0,-1,2,3,lim);
  vector<TGraphErrors*> gr = fit.GetGraphsExp("Sigma3",0,-1,2,3,lim);
  for (UInt_t i=0; i<gr.size(); ++i)
  {
    gr.at(i)->Draw("P");
    gr.at(i)->SetMarkerStyle(20);
  }
  canv->WaitPrimitive();
  // theor
  Fit::Point_t point(8);
  // Pasquini
//   point[0] = 134.;
//   point[1] = 285.;
//   point[2] = 11.5544+3.55136;
//   point[3] = 11.5544-3.55136;
//   point[4] = -4.31919;
//   point[5] = 0.849831;
//   point[6] = 3.57917;
//   point[7] = 0.900146;
  //
  //Pascalutsa
  point[0] = 134.;
  point[1] = 286.5;//333.8;
  point[2] = 9.91563e+00;
  point[3] = 7.01934e+00;
  point[4] = 1.12228e+01;
  point[5] = 3.86969e+00;
  point[6] = -4.84450e-01;
  point[7] = 1.40502e+00;
  TGraph* grThBaldin = fit.GetGraphFitFunc("PascalutsaSigma3",0,25,5.,175.,point);
//   point[0] = 134.;
//   point[1] = 285.;
//   point[2] = 11.167;
//   point[3] = 3.86921;
//   point[4] = -4.;
//   point[5] = 0.0199647;
//   point[6] = 2.;
//   point[7] = 3.;
//   TGraph* grThChPT = fit.PlotTh("PascalutsaSigma3",0,25,0.,180.,point);
// //   if (grTh)
// //   {
    grThBaldin->Draw("L");
//     grThChPT->Draw("L");
    grThBaldin->SetLineColor(kRed);
//     grThChPT->SetLineColor(kBlue);
    grThBaldin->SetLineWidth(3);
//     grThChPT->SetLineWidth(3);
    canv->WaitPrimitive();
    canv->SaveAs("sigma3.pdf");
//   }
//   else
//   {
//     cout << "Nothing to plot!" << endl;
//   }
    
//   delete canv;
}

//_____________________________________________
void test7(Fit::TFit &fit)
{
  //
    TCanvas *canv = new TCanvas();
  TH1 *h = canv->DrawFrame(0.,-0.7,180.,0.7);
//   TH1 *h = canv->DrawFrame(200.,-.5,350.,.8);
  canv->GetFrame()->SetFillColor(0);
  canv->GetFrame()->SetBorderSize(0);
  h->SetXTitle("#theta (degrees)");
  h->SetYTitle("#Sigma_{2x}");
//   h->SetYTitle("#Sigma_{3}");
  // exp
  vector<pair<Double_t,Double_t> > lim(4,make_pair(-kInfiniti,kInfiniti));
//   lim[1] = make_pair(333.8,333.8);
  vector<TGraphErrors*> gr = fit.GetGraphsExp("Sigma2x",0,-1,2,3,lim);
//   vector<TGraphErrors*> gr = fit.PlotExp("Sigma3",0,-1,2,3,lim);
  for (UInt_t i=0; i<gr.size(); ++i)
  {
    gr.at(i)->Draw("P");
    gr.at(i)->SetMarkerStyle(20);
  }
  canv->WaitPrimitive();
  // theor
  Fit::Point_t point(8);
  // Pasquini
//   point[0] = 134.;
//   point[1] = 285.;
//   point[2] = 11.5544+3.55136;
//   point[3] = 11.5544-3.55136;
//   point[4] = -4.31919;
//   point[5] = 0.849831;
//   point[6] = 3.57917;
//   point[7] = 0.900146;
  //
  //Pascalutsa
  point[0] = 134.;
  point[1] = 285;
  point[2] = 9.91563e+00;
  point[3] = 7.01934e+00;
  point[4] = 1.12228e+01;
  point[5] = 3.86969e+00;
  point[6] = -4.84450e-01;
  point[7] = 1.40502e+00;
  TGraph* grThBaldin = fit.GetGraphFitFunc("PascalutsaSigma2x",0,25,5.,175.,point);
//   point[0] = 134.;
//   point[1] = 285.;
//   point[2] = 11.167;
//   point[3] = 3.86921;
//   point[4] = -4.;
//   point[5] = 0.0199647;
//   point[6] = 2.;
//   point[7] = 3.;
//   TGraph* grThChPT = fit.PlotTh("PascalutsaSigma3",0,25,0.,180.,point);
// //   if (grTh)
// //   {
   grThBaldin->Draw("L");
//     grThChPT->Draw("L");
   grThBaldin->SetLineColor(kRed);
//     grThChPT->SetLineColor(kBlue);
   grThBaldin->SetLineWidth(3);
//     grThChPT->SetLineWidth(3);
   canv->WaitPrimitive();
//    canv->SaveAs("sigma2x.pdf");
//   }
//   else
//   {
//     cout << "Nothing to plot!" << endl;
//   }
    
//   delete canv;
}

//_____________________________________________
void test8()
{
//   Double_t p[8];
  Fit::Point_t p(8);
  p[0] = 90.;
  p[1] = 285.;
//   p[2] = 12.1952;
//   p[3] = 5.00092;
//   p[4] = -0.000866842;
//   p[5] = 1.91106;
//   p[6] = 0.516603;
//   p[7] = 0.322618;
  p[2] = 12.;
  p[3] = 1.9;
  p[4] = -4.5; // ee
  p[5] = 3.6; // mm
  p[6] = -0.6; // em
  p[7] = 2.4; // me
//   e: 285  theta: 90       a: 12.1952      b: 5.00092      g1: -0.000866842        g2: 1.91106     g3: 0.516603    g4: 0.322618
//   cout << constraintAlphaPascalutsa(p) << endl;
  cout << Fit::Sigma2xPascalutsa(p) << endl;
  cout << Fit::Sigma3Pascalutsa(p) << endl;
  cout << Fit::Sigma2xPasquini(p) << endl;
  cout << Fit::Sigma3Pasquini(p) << endl;
}

//_____________________________________________
void PrintPartialDerivatives(Fit::TFit &fit, const TString &obs, const TString &expset, const TString &suffix, Fit::Point_t &point, Double_t *err)
{
  //
  ofstream fileOut(obs+"_"+expset+"_"+suffix+".dat");
  //
  Fit::DataSetExp_t::const_iterator exp = fit.GetDataExp(obs,expset);
  if (exp == fit.GetDataExp()->end())
  {
    cout << ">>> ERROR: " << obs << " " << expset << " not found!" << endl;
    return;
  }
  const Fit::Data_t &dataExp = exp->fData;
  const Fit::TFitFunc *dataTh = exp->fTh;

  //
  fileOut << "Alpha = " << point[2] << " +- " << err[0] << " x10^-4 fm^3" << endl;
  fileOut << "Beta = " << point[3] << " +- " << err[1] << " x10^-4 fm^3" << endl;
  fileOut << "gE1E1 = " << point[4] << " x10^-4 fm^4" << endl;
  fileOut << "gM1M1 = " << point[5] << " x10^-4 fm^4" << endl;
  fileOut << "gE1M2 = " << point[6] << " x10^-4 fm^4" << endl;
  fileOut << "gM1E2 = " << point[7] << " x10^-4 fm^4" << endl;
  fileOut << endl;
//   fileOut << "sigma_g0 = " <<  "0.18 x10^-4 fm^4" << endl;
//   fileOut << "sigma_gpi = " <<  "1.8 x10^-4 fm^4" << endl;
  fileOut << endl;
  fileOut << "A = " << obs << endl;
  fileOut << endl;
  fileOut << "Eg\ttheta\texp\texp_err\t   A\t\tdA/dAlpha\tdA/dBeta\tdA/dgE1E1\tdA/dgM1M1\tdA/dgE1M2\tdA/dgM1E2" << endl;
  //
  Int_t iPoint = 1;
  for (Fit::Data_t::const_iterator iDataExp = dataExp.begin(); iDataExp != dataExp.end(); ++iDataExp)
  {
    Int_t i;
    for (i=0; i<2; ++i) point[i] = iDataExp->at(i);
    Double_t observ    = iDataExp->at(i);   // value
    Double_t observErr = iDataExp->at(i+1); // err_value
    cout << "exp point "<< iPoint << ":\t" << "e: "<< point[1] << "\ttheta: "<< point[0] << 
    "\ta: " << point[2] << 
    "\tb: "<< point[3] << 
    "\tee: "<< point[4] << 
    "\tmm: "<< point[5] << 
    "\tem: "<< point[6] << 
    "\tme: "<< point[7];
    //
//     point[1] = 50.; // MeV
    //
    Double_t observTh = dataTh->GetValue(point);
    // deriv
    Double_t deriv[6];
    for (Int_t i=0; i<6;++i)
    {
      Fit::Point_t point2 = point;
      Int_t ii = i+2;
      point2[ii] += 0.1;
      Double_t observTh2 = dataTh->GetValue(point2);
      deriv[i] = (observTh2-observTh)/(point2[ii]-point[ii]);
    }
    //
    cout << "\texp: "<< observ << "\tth: " << observTh << "\texp_err: " << observErr << endl;
    cout << "\td_alpha: "<< deriv[0] << "\td_beta: " << deriv[1] << endl;
    ++iPoint; // only for TEST
    //
    fileOut <<
    point[1] << "\t" << point[0] <<
    "\t" << observ << "\t" << observErr <<
    "\t" << observTh <<
    "\t" << deriv[0] <<
    "\t" << deriv[1] <<
    "\t" << deriv[2] <<
    "\t" << deriv[3] <<
    "\t" << deriv[4] <<
    "\t" << deriv[5] << endl;
  }
  
  //
  fileOut.close();
}

//_____________________________________________
void test9(Fit::TFit &fit)
{
  //
  Fit::Point_t point1(8), point2(8), point3(8);
  Double_t err[2];
  
  //
  point1[2] =  10.7;
  point1[3] =  3.1;
  point1[4] = -3.0;
  point1[5] =  3.5;
  point1[6] =  -2.295;
  point1[7] =  2.805;
  //
  point2[2] =  10.7;
  point2[3] =  3.1;
  point2[4] = -3.0;
  point2[5] =  3.0;
  point2[6] = -0.5;
  point2[7] =  1.51;
  //
  point3[2] =  10.7;
  point3[3] =  3.1;
  point3[4] = -3.0;
  point3[5] =  3.0;
  point3[6] = -1.4;
  point3[7] =  2.41;
  
  err[0] = err[1] = 0.;
  
//
// #pragma omp parallel
//    {
// #pragma omp sections
//      {
// #pragma omp section
//        {
	PrintPartialDerivatives(fit,"Sigma2x_Pascalutsa","mar","1",point1,err);
// 	PrintPartialDerivatives(fit,"Sigma2x_Pascalutsa","mar","2",point2,err);
// 	PrintPartialDerivatives(fit,"Sigma2x_Pascalutsa","mar","3",point3,err);
	PrintPartialDerivatives(fit,"Sigma3_Pascalutsa","bla01","1",point1,err);
// 	PrintPartialDerivatives(fit,"Sigma3_Pascalutsa","bla01","2",point2,err);
// 	PrintPartialDerivatives(fit,"Sigma3_Pascalutsa","bla01","3",point3,err);
	
//        }
// #pragma omp section
//        {
 	PrintPartialDerivatives(fit,"Sigma2x_Pasquini","mar","1",point1,err);
//        }
// #pragma omp section
//        {
	 PrintPartialDerivatives(fit,"Sigma3_Pasquini","bla01","1",point1,err);
//        }
// #pragma omp section
//        {
// 	PrintPartialDerivatives(fit,"Sigma2x_Pasquini","mar","2",point2,err);
//        }
// #pragma omp section
//        {
// 	PrintPartialDerivatives(fit,"Sigma3_Pasquini","bla01","2",point2,err);
//        }
//      }
//    }
	
//         PrintPartialDerivatives(fit,"Sigma2x_Pasquini","mar","3",point3,err);
// 	PrintPartialDerivatives(fit,"Sigma3_Pasquini","bla01","3",point3,err);

}
