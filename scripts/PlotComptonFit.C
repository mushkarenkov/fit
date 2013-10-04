#include "../src/FitTypes.h"
// #include "../src/TFitFuncLinApprox.h"

//__________________________________________________________________________________
void PlotComptonFit()
{
  // Main
  
  // Load libFit
  gROOT->ProcessLine(".L "+TString(getenv("Fit"))+"/lib/libFit.so");
    
  // Init compton theory
  Fit::InitComptonTheory();
  
  // fit object to plot exp and th data
//   Fit::TFit fit("plot.par");
  
  // Plotting style
  gROOT->ProcessLine(".x Style_Default_2.C");
  gROOT->SetStyle("Default_2");
  gStyle->SetOptStat(0);
  gStyle->SetHistFillColor(32);
  gStyle->SetHistLineColor(1);
  gStyle->SetHistFillStyle(1001);
  gStyle->SetHistLineStyle(1);
  gStyle->SetHistLineWidth(1);
  gROOT->ForceStyle();
  
  // Plot Sigma2x
//   PlotSigma2x_NoFit();
//      PlotSigma2x_FitSigma2x();
//      PlotSigma2x_FitSigma2xSigma3();
  
  // Plot Sigma3
//   PlotSigma3_NoFit();
  
  // Plot cross sections
//   PlotCS_NoFit();
  
  // TEST
  Plot_test1();
}

//__________________________________________________________________________________
void Plot_test1()
{
  // test1
  
  Fit::TFit fit("../plot.par");
  
  Fit::Point_t point1(8);
  point1[0] = 120.;
  point1[1] = 285.;
  point1[2] = 11.2461;
  point1[3] = 3.87054;
  point1[4] = -3.34375;
  point1[5] = 2.95934;
  point1[6] = 0.185504;
  point1[7] = 1.14328;
  
  //
  TGraph *gr1 = fit.GetGraphFitFunc("PascalutsaSigma3",2,40,7,17,point1);
  gr1->Draw("AL");
  gr1->SetLineColor(kRed);
  
  //
  Fit::TFitFuncLinApprox *func = (Fit::TFitFuncLinApprox*)fit.GetFitFunc("PascalutsaSigma3LinApprox");
  func->SetNvar(2);
  for (Int_t i=0; i<6; ++i) func->AddPar(point1[i+2]);
  func->AddApproxAt(point1);
  TGraph *gr2 = func->GetGraph(2,40,7,17,point1);
  gr2->Draw("L");
  gr2->SetLineColor(kGreen);
}

//__________________________________________________________________________________
void PlotCS_NoFit()
{
  // Pascalutsa
  Fit::Point_t point1(8);
  point1[2] = 11.2461;
  point1[3] = 3.87054;
  point1[4] = -3.34375;
  point1[5] = 2.95934;
  point1[6] = 0.185504;
  point1[7] = 1.14328;
  
  // Pasquini
  Fit::Point_t point2(8);
  point2[2] = 12.1;
  point2[3] = 1.6;
  point2[4] = -4.3;
  point2[5] = 2.9;
  point2[6] = 0.;
  point2[7] = 2.1;
  
  // Lvov
  Fit::Point_t point3(4);
  point2[2] = 12.;
  point2[3] = 1.9;
  
  //
  Fit::Point_t point_dummy;
  
  // Root file
//   TFile file("CS_canvas.root","RECREATE","CS canvas output");
  
  // Plot Sigma3
//   Double_t eg[12] = {333.8,323.8,310.1,298.2,286.5,275.6,265.1,254.1,244.8,234.7,224.3,213.1};
  stringstream ss;
  TString egamma;
  for (Int_t i=0; i<60; ++i)
  {
    Double_t eg = 5. + 10.*i;
    stringstream ss;
    TString egamma;
    point1[1] = point2[1] = point3[1] = eg;
    ss << eg; ss >> egamma;
    TCanvas *canv = PlotCS(eg,
			     point1, "PascalutsaCS", 50, "Pascalutsa",
			     point2, "PasquiniCS", 20, "Pasquini",
			     point3, "LvovCS", 50, "Lvov",
			     point_dummy, "", 0, "",
			     point_dummy, "", 0, "",
			     point_dummy, "", 0, ""
			    );
    if (canv) canv->SaveAs("CS_NoFit_"+egamma+"MeV.pdf");
//     file.cd();
//     canv->SetName("CS"+egamma+"MeV")
//     canv->Write();
  }
  
}

//__________________________________________________________________________________
void PlotSigma2x_NoFit()
{
  // Pascalutsa
  Fit::Point_t point1(8);
  point1[2] = 11.2461;
  point1[3] = 3.87054;
  point1[4] = -3.34375;
  point1[5] = 2.95934;
  point1[6] = 0.185504;
  point1[7] = 1.14328;
  
  // Pasquini
  Fit::Point_t point2(8);
  point2[2] = 12.1;
  point2[3] = 1.6;
  point2[4] = -4.3;
  point2[5] = 2.9;
  point2[6] = 0.;
  point2[7] = 2.1;
  
  //
  Fit::Point_t point_dummy;
  
  //
  Double_t eg = 285.;
  stringstream ss;
  TString egamma;
  point1[1] = point2[1] = eg;
  ss << eg; ss >> egamma;
  TCanvas *canv = PlotSigma2x(eg,
			     point1, "PascalutsaSigma2x", 50, "Pascalutsa (ChPT)",
			     point2, "PasquiniSigma2x", 15, "Pasquini (DR)",
			     point_dummy, "", 0, "",
			     point_dummy, "", 0, "",
			     point_dummy, "", 0, "",
			     point_dummy, "", 0, ""
			    );
  canv->SaveAs("Sigma2x_NoFit_"+egamma+"MeV.pdf");
}

//__________________________________________________________________________________
void PlotSigma2x_FitSigma2x()
{
  // Pascalutsa-1 (2)
  Fit::Point_t point1(8);
  point1[2] = 11.1712;
  point1[3] = 3.8871;
  point1[4] = -11.0987;
  point1[5] = 3;
  point1[6] = 7.60487;
  point1[7] = 1.50385;
  
  // Pascalutsa-2 (3)
  Fit::Point_t point11(8);
  point11[2] = 8.78534;
  point11[3] = 4.99999;
  point11[4] = -10.2791;
  point11[5] = 3;
  point11[6] = 6.78394;
  point11[7] = 1.50502;
  
  // Pasquini-1 (2)
  Fit::Point_t point2(8);
  point2[2] = 11.174;
  point2[3] = 3.83527;
  point2[4] = -4.48416;
  point2[5] = 3.;
  point2[6] = 1.09065;
  point2[7] = 1.40619;
  
  // Pasquini-2 (3)
  Fit::Point_t point22(8);
  point22[2] = 10.4661;
  point22[3] = 3.36451;
  point22[4] = -4.30308;
  point22[5] = 3.;
  point22[6] = 0.753415;
  point22[7] = 1.56836;
  
  //
  Fit::Point_t point_dummy;
  
  //
  Double_t eg = 285.;
  stringstream ss;
  TString egamma;
  point1[1] = point2[1] = point11[1] = point22[1] = eg;
  ss << eg; ss >> egamma;
  TCanvas *canv = PlotSigma2x(eg,
			     point1, "PascalutsaSigma2x", 50, "Pascalutsa-1",
			     point2, "PasquiniSigma2x", 15, "Pasquini-1",
			     point11, "PascalutsaSigma2x", 50, "Pascalutsa-2",
			     point22, "PasquiniSigma2x", 15, "Pasquini-2",
			     point_dummy, "", 0, "",
			     point_dummy, "", 0, ""
			    );
  canv->SaveAs("Sigma2x_FitSigma2x_"+egamma+"MeV.pdf");
}

//__________________________________________________________________________________
void PlotSigma2x_FitSigma2xSigma3()
{
  // Pascalutsa-1 (1)
  Fit::Point_t point1(8);
  point1[2] = 13.5452;
  point1[3] = 5;
  point1[4] = -17.9999;
  point1[5] = 2;
  point1[6] = -2;
  point1[7] = 5;
  
  // Pascalutsa-2 (2)
  Fit::Point_t point11(8);
  point11[2] = 11.5921;
  point11[3] = 4.6093;
  point11[4] = -2.57298;
  point11[5] = -8.31962;
  point11[6] = -0.945504;
  point11[7] = 12.8488;
  
  // Pasquini-1 (1)
  Fit::Point_t point2(8);
//   point2[2] = 11.174;
//   point2[3] = 3.83527;
//   point2[4] = -4.48416;
//   point2[5] = 3.;
//   point2[6] = 1.09065;
//   point2[7] = 1.40619;
  
  // Pasquini-2 (2)
  Fit::Point_t point22(8);
//   point22[2] = 10.4661;
//   point22[3] = 3.36451;
//   point22[4] = -4.30308;
//   point22[5] = 3.;
//   point22[6] = 0.753415;
//   point22[7] = 1.56836;
  
  //
  Fit::Point_t point_dummy;
  
  // Sigma2x
  Double_t egam = 285.;
  stringstream ss;
  TString egamma;
  point1[1] = point2[1] = point11[1] = point22[1] = egam;
  ss << egam; ss >> egamma;
  TCanvas *canv = PlotSigma2x(egam,
			     point1, "PascalutsaSigma2x", 50, "Pascalutsa-1",
			     point2, "PasquiniSigma2x", 0, "Pasquini-1",
			     point11, "PascalutsaSigma2x", 0, "Pascalutsa-2",
			     point22, "PasquiniSigma2x", 0, "Pasquini-2",
			     point_dummy, "", 0, "",
			     point_dummy, "", 0, ""
			    );
  canv->SaveAs("Sigma2x_FitSigma2xSigma3_"+egamma+"MeV.pdf");
  
  // Sigma3
  Double_t eg[12] = {333.8,323.8,310.1,298.2,286.5,275.6,265.1,254.1,244.8,234.7,224.3,213.1};
  stringstream ss;
  TString egamma;
  for (Int_t i=0; i<12; ++i)
  {
    stringstream ss;
    TString egamma;
    point1[1] = point11[1] = eg[i];
    ss << eg[i]; ss >> egamma;
    TCanvas *canv = PlotSigma3(eg[i],
			     point1, "PascalutsaSigma3", 50, "Pascalutsa-1",
			     point2, "PasquiniSigma3", 0, "Pasquini",
			     point11, "PascalutsaSigma3", 0, "Pascalutsa-2",
			     point_dummy, "", 0, "",
			     point_dummy, "", 0, "",
			     point_dummy, "", 0, ""
			    );
    canv->SaveAs("Sigma3_FitSigma2xSigma3_"+egamma+"MeV.pdf");
  }
}


//__________________________________________________________________________________
void PlotSigma2x_1()
{
  // Sigma2x, pascalutsa_2, 20130702, Global Chi2: 2.42015/3 = 0.807
  Fit::Point_t point1(8);
  point1[2] = 11.1712;
  point1[3] = 3.8871;
  point1[4] = -11.0987;
  point1[5] = 3.;
  point1[6] = 7.60487;
  point1[7] = 1.50385;
  
  // Sigma2x, pascalutsa_3, 20130702,Global Chi2: 1.83835/2 = 0.919
  Fit::Point_t point2(8);
  point2[2] = 8.78534;
  point2[3] = 4.99999;
  point2[4] = -10.2791;
  point2[5] = 3.;
  point2[6] = 6.78394;
  point2[7] = 1.50502;
  
  //
  Fit::Point_t point3;
  
  //
  Double_t eg = 285.;
  point1[1] = point2[1] = eg;
  TCanvas *canv = PlotSigma2x(eg, point1, "Pascalutsa-1", point2, "Pascalutsa-2", point3, "", "PascalutsaSigma2x");
  
  //
  canv->SaveAs("Sigma2x_Pascalutsa_1.pdf");
}

//__________________________________________________________________________________
void PlotSigma3_1()
{
  // pascalutsa_1, 20130702, Global Chi2: 77.7096/77 = 1.01
  Fit::Point_t point1(8);
  point1[2] = 13.5452;
  point1[3] = 5;
  point1[4] = -17.9999;
  point1[5] = 2;
  point1[6] = -2;
  point1[7] = 5;
  
  // pascalutsa_2, 20130702,Global Chi2: 75.3823/81 = 0.931
  Fit::Point_t point2(8);
  point2[2] = 11.5921;
  point2[3] = 4.6093;
  point2[4] = -2.57298;
  point2[5] = -8.31962;
  point2[6] = -0.945504;
  point2[7] = 12.8488;
  
  //
  Fit::Point_t point3;
  
  // Sigma3
  Double_t eg;
  stringstream ss;
  TString egamma;
  point1[1] = eg = 333.8;
  point2[1] = point1[1];
  ss << eg; ss >> egamma;
  TCanvas *canv = PlotSigma3(eg, point1, "Pascalutsa-1", point2, "Pascalutsa-2", point3, "", "PascalutsaSigma3");
  canv->SaveAs("Sigma3_Pascalutsa_1_"+egamma+"MeV.pdf");
  point1[1] = eg = 323.8;
  point2[1] = point1[1];
  stringstream ss;
  TString egamma;
  ss << eg; ss >> egamma; ss.str("");
  TCanvas *canv = PlotSigma3(eg, point1, "Pascalutsa-1", point2, "Pascalutsa-2", point3, "", "PascalutsaSigma3");
  canv->SaveAs("Sigma3_Pascalutsa_1_"+egamma+"MeV.pdf");
  stringstream ss;
  TString egamma;
  point1[1] = eg = 310.1;
  point2[1] = point1[1];
  ss << eg; ss >> egamma;
  TCanvas *canv = PlotSigma3(eg, point1, "Pascalutsa-1", point2, "Pascalutsa-2", point3, "", "PascalutsaSigma3");
  canv->SaveAs("Sigma3_Pascalutsa_1_"+egamma+"MeV.pdf");
  stringstream ss;
  TString egamma;
  point1[1] = eg = 298.2;
  point2[1] = point1[1];
  ss << eg; ss >> egamma;
  TCanvas *canv = PlotSigma3(eg, point1, "Pascalutsa-1", point2, "Pascalutsa-2", point3, "", "PascalutsaSigma3");
  canv->SaveAs("Sigma3_Pascalutsa_1_"+egamma+"MeV.pdf");
   stringstream ss;
  TString egamma;
  point1[1] = eg = 286.5;
  point2[1] = point1[1];
  ss << eg; ss >> egamma;
  TCanvas *canv = PlotSigma3(eg, point1, "Pascalutsa-1", point2, "Pascalutsa-2", point3, "", "PascalutsaSigma3");
  canv->SaveAs("Sigma3_Pascalutsa_1_"+egamma+"MeV.pdf");
   stringstream ss;
  TString egamma;
  point1[1] = eg = 275.6;
  point2[1] = point1[1];
  ss << eg; ss >> egamma;
  TCanvas *canv = PlotSigma3(eg, point1, "Pascalutsa-1", point2, "Pascalutsa-2", point3, "", "PascalutsaSigma3");
  canv->SaveAs("Sigma3_Pascalutsa_1_"+egamma+"MeV.pdf");
  stringstream ss;
  TString egamma;
  point1[1] = eg = 265.1;
  point2[1] = point1[1];
  ss << eg; ss >> egamma;
  TCanvas *canv = PlotSigma3(eg, point1, "Pascalutsa-1", point2, "Pascalutsa-2", point3, "", "PascalutsaSigma3");
  canv->SaveAs("Sigma3_Pascalutsa_1_"+egamma+"MeV.pdf");
   stringstream ss;
  TString egamma;
  point1[1] = eg = 254.1;
  point2[1] = point1[1];
  ss << eg; ss >> egamma;
  TCanvas *canv = PlotSigma3(eg, point1, "Pascalutsa-1", point2, "Pascalutsa-2", point3, "", "PascalutsaSigma3");
  canv->SaveAs("Sigma3_Pascalutsa_1_"+egamma+"MeV.pdf");
  stringstream ss;
  TString egamma;
  point1[1] = eg = 244.8;
  point2[1] = point1[1];
  ss << eg; ss >> egamma;
  TCanvas *canv = PlotSigma3(eg, point1, "Pascalutsa-1", point2, "Pascalutsa-2", point3, "", "PascalutsaSigma3");
  canv->SaveAs("Sigma3_Pascalutsa_1_"+egamma+"MeV.pdf");
  stringstream ss;
  TString egamma;
  point1[1] = eg = 234.7;
  point2[1] = point1[1];
  ss << eg; ss >> egamma;
  TCanvas *canv = PlotSigma3(eg, point1, "Pascalutsa-1", point2, "Pascalutsa-2", point3, "", "PascalutsaSigma3");
  canv->SaveAs("Sigma3_Pascalutsa_1_"+egamma+"MeV.pdf");
  stringstream ss;
  TString egamma;
  point1[1] = eg = 224.3;
  point2[1] = point1[1];
  ss << eg; ss >> egamma;
  TCanvas *canv = PlotSigma3(eg, point1, "Pascalutsa-1", point2, "Pascalutsa-2", point3, "", "PascalutsaSigma3");
  canv->SaveAs("Sigma3_Pascalutsa_1_"+egamma+"MeV.pdf");
  stringstream ss;
  TString egamma;
  point1[1] = eg = 213.1;
  point2[1] = point1[1];
  ss << eg; ss >> egamma;
  TCanvas *canv = PlotSigma3(eg, point1, "Pascalutsa-1", point2, "Pascalutsa-2", point3, "", "PascalutsaSigma3");
  canv->SaveAs("Sigma3_Pascalutsa_1_"+egamma+"MeV.pdf");
  
  // Sigma2x
  Double_t eg = 285.;
  point1[1] = point2[1] = eg;
  TCanvas *canv = PlotSigma2x(eg, point1, "Pascalutsa-1", point2, "Pascalutsa-2", point3, "", "PascalutsaSigma2x");
  canv->SaveAs("Sigma2x_Pascalutsa_1_285MeV.pdf");
  
}

//__________________________________________________________________________________
void PlotSigma3_NoFit()
{
  // Pascalutsa
  Fit::Point_t point1(8);
  point1[2] = 11.2461;
  point1[3] = 3.87054;
  point1[4] = -3.34375;
  point1[5] = 2.95934;
  point1[6] = 0.185504;
  point1[7] = 1.14328;
  
  // Pasquini
  Fit::Point_t point2(8);
  point2[2] = 12.1;
  point2[3] = 1.6;
  point2[4] = -4.3;
  point2[5] = 2.9;
  point2[6] = 0.;
  point2[7] = 2.1;
  
  // Lvov
  Fit::Point_t point3(4);
  point2[2] = 12.;
  point2[3] = 1.9;
  
  //
  Fit::Point_t point_dummy;
  
  // Plot Sigma3
  Double_t eg[12] = {333.8,323.8,310.1,298.2,286.5,275.6,265.1,254.1,244.8,234.7,224.3,213.1};
  stringstream ss;
  TString egamma;
  for (Int_t i=0; i<12; ++i)
  {
    stringstream ss;
    TString egamma;
    point1[1] = point2[1] = point3[1] = eg[i];
    ss << eg[i]; ss >> egamma;
    TCanvas *canv = PlotSigma3(eg[i],
			     point1, "PascalutsaSigma3", 50, "Pascalutsa (ChPT)",
			     point2, "PasquiniSigma3", 15, "Pasquini (DR)",
			     point3, "LvovSigma3", 50, "Lvov (DR)",
			     point_dummy, "", 0, "",
			     point_dummy, "", 0, "",
			     point_dummy, "", 0, ""
			    );
    canv->SaveAs("Sigma3_NoFit_"+egamma+"MeV.pdf");
  }
  
}

//__________________________________________________________________________________
void PlotSigma3()
{
  // Sigma2x Sigma3
  Fit::Point_t pointPascalutsa(8);
  pointPascalutsa[0] = 0.;
  pointPascalutsa[1] = 244.8;
  pointPascalutsa[2] = 12.9892;
  pointPascalutsa[3] = 5.;
  pointPascalutsa[4] = -1.79461;
  pointPascalutsa[5] = 2.32828e-06;
  pointPascalutsa[6] = -1.99999;
  pointPascalutsa[7] = 4.80911;
  //
  Fit::Point_t pointPasquini = pointPascalutsa;
//   pointPasquini[0] = 
//   pointPasquini[0] = 
//   pointPasquini[0] = 
//   pointPasquini[0] = 
//   pointPasquini[0] = 
//   pointPasquini[0] = 
//   pointPasquini[0] = 
//   pointPasquini[0] = 
  //
  Fit::Point_t pointLvov(4);
  pointLvov[0] = 0.;
  pointLvov[1] = 244.8;
  pointLvov[2] = 12.9892;
  pointLvov[3] = 5.;
  //
//   pointPascalutsa[1] = pointPasquini[1] = pointLvov[1] = 333.8;
  pointPascalutsa.clear();
  pointPasquini.clear();
  pointLvov[1] = 333.8;
  PlotSigma3(pointLvov[1], pointPascalutsa, pointPasquini, pointLvov);
//   pointPascalutsa[1] = pointPasquini[1] = pointLvov[1] = 323.8;
//   PlotSigma3(pointPascalutsa[1], pointPascalutsa, pointPasquini, pointLvov);
//   pointPascalutsa[1] = pointPasquini[1] = pointLvov[1] = 310.1;
//   PlotSigma3(pointPascalutsa[1], pointPascalutsa, pointPasquini, pointLvov);
//   pointPascalutsa[1] = pointPasquini[1] = pointLvov[1] = 298.2;
//   PlotSigma3(pointPascalutsa[1], pointPascalutsa, pointPasquini, pointLvov);
//   pointPascalutsa[1] = pointPasquini[1] = pointLvov[1] = 286.5;
//   PlotSigma3(pointPascalutsa[1], pointPascalutsa, pointPasquini, pointLvov);
//   pointPascalutsa[1] = pointPasquini[1] = pointLvov[1] = 275.6;
//   PlotSigma3(pointPascalutsa[1], pointPascalutsa, pointPasquini, pointLvov);
//   pointPascalutsa[1] = pointPasquini[1] = pointLvov[1] = 265.1;
//   PlotSigma3(pointPascalutsa[1], pointPascalutsa, pointPasquini, pointLvov);
//   pointPascalutsa[1] = pointPasquini[1] = pointLvov[1] = 254.1;
//   PlotSigma3(pointPascalutsa[1], pointPascalutsa, pointPasquini, pointLvov);
//   pointPascalutsa[1] = pointPasquini[1] = pointLvov[1] = 244.8;
//   PlotSigma3(pointPascalutsa[1], pointPascalutsa, pointPasquini, pointLvov);
//   pointPascalutsa[1] = pointPasquini[1] = pointLvov[1] = 234.7;
//   PlotSigma3(pointPascalutsa[1], pointPascalutsa, pointPasquini, pointLvov);
//   pointPascalutsa[1] = pointPasquini[1] = pointLvov[1] = 224.3;
//   PlotSigma3(pointPascalutsa[1], pointPascalutsa, pointPasquini, pointLvov);
//   pointPascalutsa[1] = pointPasquini[1] = pointLvov[1] = 213.1;
//   PlotSigma3(pointPascalutsa[1], pointPascalutsa, pointPasquini, pointLvov);
}

//__________________________________________________________________________________
void PlotSigma2x(const TString &nameFileResultPascalutsa, const TString &nameFileResultPasquini)
{
  //
  Fit::Point_t pointPascalutsa;
  if (nameFileResultPascalutsa != "")
  {
    TFile file(nameFileResultPascalutsa);
    Fit::TFitResult *result = (Fit::TFitResult*)file->Get("Result");
    if (!result) return;
    vector<Fit::Par_t> params = result->fParams;
    pointPascalutsa.push_back(0.); // can be any
    pointPascalutsa.push_back(285.);
    for (UInt_t i=0;i<params.size();++i) pointPascalutsa.push_back(params.at(i));
  }
  //
  Fit::Point_t pointPasquini;
  if (nameFileResultPasquini != "")
  {
    TFile file(nameFileResultPasquini);
    Fit::TFitResult *result = (Fit::TFitResult*)file->Get("Result");
    if (!result) return;
    vector<Fit::Par_t> params = result->fParams;
    pointPasquini.push_back(0.); // can be any
    pointPasquini.push_back(285.);
    for (UInt_t i=0;i<params.size();++i) pointPasquini.push_back(params.at(i));
  }
  //
  PlotSigma2x(pointPascalutsa, pointPasquini);
}


//__________________________________________________________________________________
TCanvas *PlotSigma2x(const Fit::Point_t &pointPascalutsa, const Fit::Point_t &pointPasquini)
{
  // Plot Sigma2x
  
  //
  Fit::TFit fit("../plot.par");
  
  // Canvas
  TCanvas *canv = new TCanvas();
  TH1 *h = canv->DrawFrame(0.,-0.7,180.,0.7);
  canv->GetFrame()->SetFillColor(0);
  canv->GetFrame()->SetBorderSize(0);
  h->SetXTitle("#theta (degrees)");
  h->SetYTitle("#Sigma_{2x}");
    
  // Phil's data
  Double_t min[4], max[4];
  for (Int_t i=0; i<4; ++i)
  {
    min[i] = -1e+10;
    max[i] =  1e+10;
  }
  TGraphErrors *grExp = fit.GetGraphExp("Sigma2x","mar", 0,-1,2,3,min,max,4);
  if (grExp)
  {
    grExp->Draw("P");
    grExp->SetMarkerStyle(20);
  }
  
  // Sigma2x Pascalutsa
  if (!pointPascalutsa.empty())
  {
    TGraph* grThPascalutsa = fit.GetGraphFitFunc("PascalutsaSigma2x",0,25,5.,175.,pointPascalutsa);
    if (grThPascalutsa)
    {
      grThPascalutsa->Draw("L");
      grThPascalutsa->SetLineColor(kBlue);
      grThPascalutsa->SetLineWidth(3);
    }
  }
  
  // Sigma2x Pasquini
  if (!pointPasquini.empty())
  {
    TGraph* grThPasquini = fit.GetGraphFitFunc("Sigma2xPasquini",0,15,0.,180.,pointPasquini);
    if (grThPasquini)
    {
      grThPasquini->Draw("L");
      grThPasquini->SetLineColor(kRed);
      grThPasquini->SetLineWidth(3);
    }
  }
  
  //
  return canv;
  
}

//__________________________________________________________________________________
TCanvas *PlotSigma2x(const Double_t &eg,
		     const Fit::Point_t &point1, const TString &leg1,
		     const Fit::Point_t &point2, const TString &leg2,
		     const Fit::Point_t &point3, const TString &leg3,
		     const TString &nameFitFunc)
{
  // Plot Sigma2x
  
  //
  Fit::TFit fit("../plot.par");
  
  // Canvas
  TCanvas *canv = new TCanvas();
  TH1 *h = canv->DrawFrame(0.,-0.7,180.,0.7);
  canv->GetFrame()->SetFillColor(0);
  canv->GetFrame()->SetBorderSize(0);
  h->SetXTitle("#theta (degrees)");
  h->SetYTitle("#Sigma_{2x}");
  
  // Legend
  TLegend *lg = new TLegend(0.55,0.75,0.89,0.89,"","brNDC");
  lg->SetLineColor(0);
  lg->SetTextSize(0.03);
    
  // Phil's data
  Double_t min[4], max[4];
  for (Int_t i=0; i<4; ++i)
  {
    min[i] = -1e+10;
    max[i] =  1e+10;
  }
  TGraphErrors *grExp = fit.GetGraphExp("Sigma2x","mar", 0,-1,2,3,min,max,4);
  if (grExp)
  {
    grExp->Draw("P");
    grExp->SetMarkerStyle(20);
    lg->AddEntry(grExp,"CB@MAMI","p");
  }
  
  // point1
  if (!point1.empty())
  {
    TGraph* grTh1 = fit.GetGraphFitFunc(nameFitFunc,0,50,5.,175.,point1);
    if (grTh1)
    {
      grTh1->Draw("L");
      grTh1->SetLineColor(kBlue);
      grTh1->SetLineWidth(3);
      lg->AddEntry(grTh1,leg1,"l");
    }
  }
  
  // point2
  if (!point2.empty())
  {
    TGraph* grTh2 = fit.GetGraphFitFunc(nameFitFunc,0,50,5.,175.,point2);
    if (grTh2)
    {
      grTh2->Draw("L");
      grTh2->SetLineColor(kRed);
      grTh2->SetLineWidth(3);
      lg->AddEntry(grTh2,leg2,"l");
    }
  }
  
  // point3
  if (!point3.empty())
  {
    TGraph* grTh3 = fit.GetGraphFitFunc(nameFitFunc,0,50,5.,175.,point3);
    if (grTh3)
    {
      grTh3->Draw("L");
      grTh3->SetLineColor(kGreen);
      grTh3->SetLineWidth(3);
      lg->AddEntry(grTh3,leg3,"l");
    }
  }
  
  //
  lg->Draw();
  
  //
  return canv;
  
}

//__________________________________________________________________________________
void PlotSigma3(const TString &nameFileResultPascalutsa, const TString &nameFileResultPasquini)
{
  //
  Fit::Point_t pointPascalutsa;
  if (nameFileResultPascalutsa != "")
  {
    TFile file(nameFileResultPascalutsa);
    Fit::TFitResult *result = (Fit::TFitResult*)file->Get("Result");
    if (!result) return;
    vector<Fit::Par_t> params = result->fParams;
    pointPascalutsa.push_back(0.); // can be any
    for (UInt_t i=0;i<params.size();++i) pointPascalutsa.push_back(params.at(i));
  }
  //
  Fit::Point_t pointPasquini;
  if (nameFileResultPasquini != "")
  {
    TFile file(nameFileResultPasquini);
    Fit::TFitResult *result = (Fit::TFitResult*)file->Get("Result");
    if (!result) return;
    vector<Fit::Par_t> params = result->fParams;
    pointPasquini.push_back(0.); // can be any
    for (UInt_t i=0;i<params.size();++i) pointPasquini.push_back(params.at(i));
  }
  //
  PlotSigma3(285.,pointPascalutsa, pointPasquini);
}

//__________________________________________________________________________________
void PlotSigma3(const Double_t &eg, const Fit::Point_t &pointPascalutsa, const Fit::Point_t &pointPasquini, const Fit::Point_t &pointLvov)
{
  // Plot Sigma3
  
  //
  Fit::TFit fit("../plot.par");
  
  // Canvas
  TCanvas *canv = new TCanvas();
  TH1 *h = canv->DrawFrame(0.,-0.7,180.,0.7);
  canv->GetFrame()->SetFillColor(0);
  canv->GetFrame()->SetBorderSize(0);
  h->SetXTitle("#theta (degrees)");
  h->SetYTitle("#Sigma_{3}");
    
  // LEGS data
  Double_t min[4], max[4];
  for (Int_t i=0; i<4; ++i)
  {
    min[i] = -1e+10;
    max[i] =  1e+10;
  }
  min[1] = max[1] = eg;
  TGraphErrors *grExp = fit.GetGraphExp("Sigma3","bla01", 0,-1,2,3,min,max,4);
  if (grExp)
  {
    grExp->Draw("P");
    grExp->SetMarkerStyle(20);
  }
  
  // Pascalutsa
  if (!pointPascalutsa.empty())
  {
    TGraph* grThPascalutsa = fit.GetGraphFitFunc("PascalutsaSigma3",0,25,0.,180.,pointPascalutsa);
    if (grThPascalutsa)
    {
      grThPascalutsa->Draw("L");
      grThPascalutsa->SetLineColor(kBlue);
      grThPascalutsa->SetLineWidth(3);
    }
  }
  
  // Pasquini
  if (!pointPasquini.empty())
  {
    TGraph* grThPasquini = fit.GetGraphFitFunc("PasquiniSigma3",0,10,0.,180.,pointPasquini);
    if (grThPasquini)
    {
      grThPasquini->Draw("L");
      grThPasquini->SetLineColor(kRed);
      grThPasquini->SetLineWidth(3);
    }
  }
  
  // Lvov
  if (!pointLvov.empty())
  {
    TGraph* grThLvov = fit.GetGraphFitFunc("LvovSigma3",0,100,0.,180.,pointLvov);
    if (grThLvov)
    {
      grThLvov->Draw("L");
      grThLvov->SetLineColor(kGreen);
      grThLvov->SetLineWidth(3);
    }
  }
  
}

//__________________________________________________________________________________
TCanvas *PlotSigma3(const Double_t &eg,
		    const Fit::Point_t &point1, const TString &nameFitFunc1, const Int_t n1, const TString &leg1,
		    const Fit::Point_t &point2, const TString &nameFitFunc2, const Int_t n2, const TString &leg2,
		    const Fit::Point_t &point3, const TString &nameFitFunc3, const Int_t n3, const TString &leg3,
		    const Fit::Point_t &point4, const TString &nameFitFunc4, const Int_t n4, const TString &leg4,
		    const Fit::Point_t &point5, const TString &nameFitFunc5, const Int_t n5, const TString &leg5,
		    const Fit::Point_t &point6, const TString &nameFitFunc6, const Int_t n6, const TString &leg6
		   )
{
  // Plot Sigma3
  
  //
  Fit::TFit fit("../plot.par");
  
  // Canvas
  TCanvas *canv = new TCanvas();
  TH1 *h = canv->DrawFrame(0.,-0.7,180.,0.7);
  canv->GetFrame()->SetFillColor(0);
  canv->GetFrame()->SetBorderSize(0);
  h->SetXTitle("#theta (degrees)");
  h->SetYTitle("#Sigma_{3}");
  
  // Legend
  TLegend *lg = new TLegend(0.55,0.75,0.89,0.89,"","brNDC");
  lg->SetLineColor(0);
  lg->SetTextSize(0.03);
  
  // Text Eg
  stringstream ss;
  TString egamma;
  ss << eg; ss >> egamma;
  TText text;
  text.DrawTextNDC(0.43,0.93,egamma+" MeV");
    
  // Phil's data
  Double_t min[4], max[4];
  for (Int_t i=0; i<4; ++i)
  {
    min[i] = -1e+10;
    max[i] =  1e+10;
  }
  min[1] = max[1] = eg;
  TGraphErrors *grExp = fit.GetGraphExp("Sigma3","bla01", 0,-1,2,3,min,max,4);
  if (grExp)
  {
    grExp->Draw("P");
    grExp->SetMarkerStyle(20);
    lg->AddEntry(grExp,"LEGS","p");
  }
  
  // point1
  if (n1>1)
  {
    TGraph* grTh1 = fit.GetGraphFitFunc(nameFitFunc1,0,n1,5.,175.,point1);
    if (grTh1)
    {
      grTh1->Draw("L");
      grTh1->SetLineColor(kBlue);
      grTh1->SetLineWidth(3);
      lg->AddEntry(grTh1,leg1,"l");
    }
  }
  
  // point2
  if (n2>1)
  {
    TGraph* grTh2 = fit.GetGraphFitFunc(nameFitFunc2,0,n2,5.,175.,point2);
    if (grTh2)
    {
      grTh2->Draw("L");
      grTh2->SetLineColor(kRed);
      grTh2->SetLineWidth(3);
      lg->AddEntry(grTh2,leg2,"l");
    }
  }
  
  // point3
  if (n3>1)
  {
    TGraph* grTh3 = fit.GetGraphFitFunc(nameFitFunc3,0,n3,5.,175.,point3);
    if (grTh3)
    {
      grTh3->Draw("L");
      grTh3->SetLineColor(kGreen);
      grTh3->SetLineWidth(3);
      lg->AddEntry(grTh3,leg3,"l");
    }
  }
  
    // point4
  if (n4>1)
  {
    TGraph* grTh4 = fit.GetGraphFitFunc(nameFitFunc4,0,n4,5.,175.,point4);
    if (grTh4)
    {
      grTh4->Draw("L");
      grTh4->SetLineColor(kBlue);
      grTh4->SetLineWidth(3);
      grTh4->SetLineStyle(2);
      lg->AddEntry(grTh4,leg4,"l");
    }
  }
  
  // point5
  if (n5>1)
  {
    TGraph* grTh5 = fit.GetGraphFitFunc(nameFitFunc5,0,n5,5.,175.,point5);
    if (grTh5)
    {
      grTh5->Draw("L");
      grTh5->SetLineColor(kRed);
      grTh5->SetLineWidth(3);
      grTh5->SetLineStyle(2);
      lg->AddEntry(grTh5,leg5,"l");
    }
  }
  
  // point6
  if (n6>1)
  {
    TGraph* grTh6 = fit.GetGraphFitFunc(nameFitFunc6,0,n6,5.,175.,point6);
    if (grTh6)
    {
      grTh6->Draw("L");
      grTh6->SetLineColor(kGreen);
      grTh6->SetLineWidth(3);
      grTh6->SetLineStyle(2);
      lg->AddEntry(grTh6,leg6,"l");
    }
  }
  
  //
  lg->Draw();
  
  //
  return canv;
  
}

//__________________________________________________________________________________
TCanvas *PlotSigma2x(const Double_t &eg,
		    const Fit::Point_t &point1, const TString &nameFitFunc1, const Int_t n1, const TString &leg1,
		    const Fit::Point_t &point2, const TString &nameFitFunc2, const Int_t n2, const TString &leg2,
		    const Fit::Point_t &point3, const TString &nameFitFunc3, const Int_t n3, const TString &leg3,
		    const Fit::Point_t &point4, const TString &nameFitFunc4, const Int_t n4, const TString &leg4,
		    const Fit::Point_t &point5, const TString &nameFitFunc5, const Int_t n5, const TString &leg5,
		    const Fit::Point_t &point6, const TString &nameFitFunc6, const Int_t n6, const TString &leg6
		   )
{
  // Plot Sigma2x
  
  //
  Fit::TFit fit("../plot.par");
  
  // Canvas
  TCanvas *canv = new TCanvas();
  TH1 *h = canv->DrawFrame(0.,-0.7,180.,0.7);
  canv->GetFrame()->SetFillColor(0);
  canv->GetFrame()->SetBorderSize(0);
  h->SetXTitle("#theta (degrees)");
  h->SetYTitle("#Sigma_{2x}");
  
  // Legend
  TLegend *lg = new TLegend(0.55,0.75,0.89,0.89,"","brNDC");
  lg->SetLineColor(0);
  lg->SetTextSize(0.03);
  
  // Text Eg
  stringstream ss;
  TString egamma;
  ss << eg; ss >> egamma;
  TText text;
  text.DrawTextNDC(0.43,0.93,egamma+" MeV");
    
  // Phil's data
  Double_t min[4], max[4];
  for (Int_t i=0; i<4; ++i)
  {
    min[i] = -1e+10;
    max[i] =  1e+10;
  }
  min[1] = max[1] = eg;
  TGraphErrors *grExp = fit.GetGraphExp("Sigma2x","mar", 0,-1,2,3,min,max,4);
  if (grExp)
  {
    grExp->Draw("P");
    grExp->SetMarkerStyle(20);
    lg->AddEntry(grExp,"CB@MAMI","p");
  }
  
  // point1
  if (n1>1)
  {
    TGraph* grTh1 = fit.GetGraphFitFunc(nameFitFunc1,0,n1,5.,175.,point1);
    if (grTh1)
    {
      grTh1->Draw("L");
      grTh1->SetLineColor(kBlue);
      grTh1->SetLineWidth(3);
      lg->AddEntry(grTh1,leg1,"l");
    }
  }
  
  // point2
  if (n2>1)
  {
    TGraph* grTh2 = fit.GetGraphFitFunc(nameFitFunc2,0,n2,5.,175.,point2);
    if (grTh2)
    {
      grTh2->Draw("L");
      grTh2->SetLineColor(kRed);
      grTh2->SetLineWidth(3);
      lg->AddEntry(grTh2,leg2,"l");
    }
  }
  
  // point3
  if (n3>1)
  {
    TGraph* grTh3 = fit.GetGraphFitFunc(nameFitFunc3,0,n3,5.,175.,point3);
    if (grTh3)
    {
      grTh3->Draw("L");
      grTh3->SetLineColor(kGreen);
      grTh3->SetLineWidth(3);
      lg->AddEntry(grTh3,leg3,"l");
    }
  }
  
    // point4
  if (n4>1)
  {
    TGraph* grTh4 = fit.GetGraphFitFunc(nameFitFunc4,0,n4,5.,175.,point4);
    if (grTh4)
    {
      grTh4->Draw("L");
      grTh4->SetLineColor(kBlue);
      grTh4->SetLineWidth(3);
      grTh4->SetLineStyle(2);
      lg->AddEntry(grTh4,leg4,"l");
    }
  }
  
  // point5
  if (n5>1)
  {
    TGraph* grTh5 = fit.GetGraphFitFunc(nameFitFunc5,0,n5,5.,175.,point5);
    if (grTh5)
    {
      grTh5->Draw("L");
      grTh5->SetLineColor(kRed);
      grTh5->SetLineWidth(3);
      grTh5->SetLineStyle(2);
      lg->AddEntry(grTh5,leg5,"l");
    }
  }
  
  // point6
  if (n6>1)
  {
    TGraph* grTh6 = fit.GetGraphFitFunc(nameFitFunc6,0,n6,5.,175.,point6);
    if (grTh6)
    {
      grTh6->Draw("L");
      grTh6->SetLineColor(kGreen);
      grTh6->SetLineWidth(3);
      grTh6->SetLineStyle(2);
      lg->AddEntry(grTh6,leg6,"l");
    }
  }
  
  //
  lg->Draw();
  
  //
  return canv;
  
}

//__________________________________________________________________________________
TCanvas *PlotCS(const Double_t &eg,
		    const Fit::Point_t &point1, const TString &nameFitFunc1, Int_t n1, const TString &leg1,
		    const Fit::Point_t &point2, const TString &nameFitFunc2, Int_t n2, const TString &leg2,
		    const Fit::Point_t &point3, const TString &nameFitFunc3, Int_t n3, const TString &leg3,
		    const Fit::Point_t &point4, const TString &nameFitFunc4, Int_t n4, const TString &leg4,
		    const Fit::Point_t &point5, const TString &nameFitFunc5, Int_t n5, const TString &leg5,
		    const Fit::Point_t &point6, const TString &nameFitFunc6, Int_t n6, const TString &leg6
		   )
{
  // Plot SigmaCS
  
  //
  Fit::TFit fit("../plot.par");
  
  // Canvas
  TCanvas *canv = new TCanvas();
  canv->SetRightMargin(0.15);
  Int_t csmax = 400.;
  if (eg < 100.) csmax = 30.;
  else if (eg < 170.) csmax = 50.;
  else if (eg < 220.) csmax = 100.;
  else if (eg < 300.) csmax = 300.;
  TH1 *h = canv->DrawFrame(0.,0.,180.,csmax);
  canv->GetFrame()->SetFillColor(0);
  canv->GetFrame()->SetBorderSize(0);
  h->SetXTitle("#theta (degrees)");
  h->SetYTitle("d#sigma/d#Omega, nb/sr");
  
  // Legend
  TLegend *lg = new TLegend(0.86,0.5,1.,0.89,"","brNDC");
  lg->SetLineColor(0);
  lg->SetTextSize(0.03);
  
  // Text Eg
  Double_t de = 5;
  stringstream ss;
  TString egamma;
  ss << eg; ss >> egamma;
  stringstream ss;
  TString degamma;
  ss << de; ss >> degamma;
  TLatex text;
  text.SetNDC();
  text.DrawLatex(0.43,0.93,egamma+" #pm "+degamma+" MeV");
    
  // Phil's data
  Double_t min[4], max[4];
  for (Int_t i=0; i<4; ++i)
  {
    min[i] = -1e+10;
    max[i] =  1e+10;
  }
  min[1] = max[1] = eg;
  min[1] = eg - de;
  max[1] = eg + de;
  Int_t nPointsExp = 0;
  for (Int_t i=0; i<fit.GetDataExp()->size(); ++i)
  {
    Fit::DataExp_t &data = fit.GetDataExp()->at(i);
    if (data.fKeyObserv != "sigma") continue;
    TGraphErrors *grExp = fit.GetGraphExp("sigma",data.fKeyData, 0,-1,2,3,min,max,4);
    if (grExp)
    {
      nPointsExp += grExp->GetN();
      grExp->Draw("P");
      grExp->SetMarkerStyle(20);
      Int_t icolor = (i==0||i==10||i==18||i==19) ? i+31 : i;
      grExp->SetMarkerColor(icolor);
      lg->AddEntry(grExp,data.fKeyData,"p");
    }
  }
  
  // Don't plot theory when no exp
  if (nPointsExp < 15)
  {
    delete lg;
    delete canv;
    return NULL;
  }
  
  // point1
  if (n1>1)
  {
    TGraph* grTh1 = fit.GetGraphFitFunc(nameFitFunc1,0,n1,5.,175.,point1);
    if (grTh1)
    {
      grTh1->Draw("L");
      grTh1->SetLineColor(kBlue);
      grTh1->SetLineWidth(3);
      lg->AddEntry(grTh1,leg1,"l");
    }
  }
  
  // point2
  if (n2>1)
  {
    TGraph* grTh2 = fit.GetGraphFitFunc(nameFitFunc2,0,n2,5.,175.,point2);
    if (grTh2)
    {
      grTh2->Draw("L");
      grTh2->SetLineColor(kRed);
      grTh2->SetLineWidth(3);
      lg->AddEntry(grTh2,leg2,"l");
    }
  }
  
  // point3
  if (n3>1)
  {
    TGraph* grTh3 = fit.GetGraphFitFunc(nameFitFunc3,0,n3,5.,175.,point3);
    if (grTh3)
    {
      grTh3->Draw("L");
      grTh3->SetLineColor(kGreen);
      grTh3->SetLineWidth(3);
      lg->AddEntry(grTh3,leg3,"l");
    }
  }
  
    // point4
  if (n4>1)
  {
    TGraph* grTh4 = fit.GetGraphFitFunc(nameFitFunc4,0,n4,5.,175.,point4);
    if (grTh4)
    {
      grTh4->Draw("L");
      grTh4->SetLineColor(kBlue);
      grTh4->SetLineWidth(3);
      grTh4->SetLineStyle(2);
      lg->AddEntry(grTh4,leg4,"l");
    }
  }
  
  // point5
  if (n5>1)
  {
    TGraph* grTh5 = fit.GetGraphFitFunc(nameFitFunc5,0,n5,5.,175.,point5);
    if (grTh5)
    {
      grTh5->Draw("L");
      grTh5->SetLineColor(kRed);
      grTh5->SetLineWidth(3);
      grTh5->SetLineStyle(2);
      lg->AddEntry(grTh5,leg5,"l");
    }
  }
  
  // point6
  if (n6>1)
  {
    TGraph* grTh6 = fit.GetGraphFitFunc(nameFitFunc6,0,n6,5.,175.,point6);
    if (grTh6)
    {
      grTh6->Draw("L");
      grTh6->SetLineColor(kGreen);
      grTh6->SetLineWidth(3);
      grTh6->SetLineStyle(2);
      lg->AddEntry(grTh6,leg6,"l");
    }
  }
  
  //
  lg->Draw();
  
  //
  return canv;
  
}

// //__________________________________________________________________________________
// TCanvas *PlotCS_E(const Double_t &th,
// 		    const Fit::Point_t &point1, const TString &nameFitFunc1, Int_t n1, const TString &leg1,
// 		    const Fit::Point_t &point2, const TString &nameFitFunc2, Int_t n2, const TString &leg2,
// 		    const Fit::Point_t &point3, const TString &nameFitFunc3, Int_t n3, const TString &leg3,
// 		    const Fit::Point_t &point4, const TString &nameFitFunc4, Int_t n4, const TString &leg4,
// 		    const Fit::Point_t &point5, const TString &nameFitFunc5, Int_t n5, const TString &leg5,
// 		    const Fit::Point_t &point6, const TString &nameFitFunc6, Int_t n6, const TString &leg6
// 		   )
// {
//   // Plot SigmaCS (Eg)
//   
//   //
//   Fit::TFit fit("../plot.par");
//   
//   // Canvas
//   TCanvas *canv = new TCanvas();
//   canv->SetRightMargin(0.15);
//   Int_t csmax = 200.;
//   if (eg < 100.) csmax = 30.;
//   else if (eg < 170.) csmax = 50;
//   else if (eg < 250.) csmax 100.;
//   TH1 *h = canv->DrawFrame(0.,0.,180.,csmax);
//   canv->GetFrame()->SetFillColor(0);
//   canv->GetFrame()->SetBorderSize(0);
//   h->SetXTitle("#theta (degrees)");
//   h->SetYTitle("d#sigma/d#Omega, nb/sr");
//   
//   // Legend
//   TLegend *lg = new TLegend(0.86,0.5,1.,0.89,"","brNDC");
//   lg->SetLineColor(0);
//   lg->SetTextSize(0.03);
//   
//   // Text Eg
//   Double_t de = 5;
//   stringstream ss;
//   TString egamma;
//   ss << eg; ss >> egamma;
//   stringstream ss;
//   TString degamma;
//   ss << de; ss >> degamma;
//   TLatex text;
//   text.SetNDC();
//   text.DrawLatex(0.43,0.93,egamma+" #pm "+degamma+" MeV");
//     
//   // Phil's data
//   Double_t min[4], max[4];
//   for (Int_t i=0; i<4; ++i)
//   {
//     min[i] = -1e+10;
//     max[i] =  1e+10;
//   }
//   min[1] = max[1] = eg;
//   min[1] = eg - de;
//   max[1] = eg + de;
//   Int_t nPointsExp = 0;
//   for (Int_t i=0; i<fit.GetDataExp()->size(); ++i)
//   {
//     Fit::DataExp_t &data = fit.GetDataExp()->at(i);
//     if (data.fKeyObserv != "sigma") continue;
//     TGraphErrors *grExp = fit.GetGraphExp("sigma",data.fKeyData, 0,-1,2,3,min,max,4);
//     if (grExp)
//     {
//       nPointsExp += grExp->GetN();
//       grExp->Draw("P");
//       grExp->SetMarkerStyle(20);
//       grExp->SetMarkerColor(19+i);
//       lg->AddEntry(grExp,data.fKeyData,"p");
//     }
//   }
//   
//   // Don't plot theory when no exp
//   if (nPointsExp < 15)
//   {
//     delete lg;
//     delete canv;
//     return NULL;
//   }
//   
//   // point1
//   if (n1>1)
//   {
//     TGraph* grTh1 = fit.GetGraphFitFunc(nameFitFunc1,0,n1,5.,175.,point1);
//     if (grTh1)
//     {
//       grTh1->Draw("L");
//       grTh1->SetLineColor(kBlue);
//       grTh1->SetLineWidth(3);
//       lg->AddEntry(grTh1,leg1,"l");
//     }
//   }
//   
//   // point2
//   if (n2>1)
//   {
//     TGraph* grTh2 = fit.GetGraphFitFunc(nameFitFunc2,0,n2,5.,175.,point2);
//     if (grTh2)
//     {
//       grTh2->Draw("L");
//       grTh2->SetLineColor(kRed);
//       grTh2->SetLineWidth(3);
//       lg->AddEntry(grTh2,leg2,"l");
//     }
//   }
//   
//   // point3
//   if (n3>1)
//   {
//     TGraph* grTh3 = fit.GetGraphFitFunc(nameFitFunc3,0,n3,5.,175.,point3);
//     if (grTh3)
//     {
//       grTh3->Draw("L");
//       grTh3->SetLineColor(kGreen);
//       grTh3->SetLineWidth(3);
//       lg->AddEntry(grTh3,leg3,"l");
//     }
//   }
//   
//     // point4
//   if (n4>1)
//   {
//     TGraph* grTh4 = fit.GetGraphFitFunc(nameFitFunc4,0,n4,5.,175.,point4);
//     if (grTh4)
//     {
//       grTh4->Draw("L");
//       grTh4->SetLineColor(kBlue);
//       grTh4->SetLineWidth(3);
//       grTh4->SetLineStyle(2);
//       lg->AddEntry(grTh4,leg4,"l");
//     }
//   }
//   
//   // point5
//   if (n5>1)
//   {
//     TGraph* grTh5 = fit.GetGraphFitFunc(nameFitFunc5,0,n5,5.,175.,point5);
//     if (grTh5)
//     {
//       grTh5->Draw("L");
//       grTh5->SetLineColor(kRed);
//       grTh5->SetLineWidth(3);
//       grTh5->SetLineStyle(2);
//       lg->AddEntry(grTh5,leg5,"l");
//     }
//   }
//   
//   // point6
//   if (n6>1)
//   {
//     TGraph* grTh6 = fit.GetGraphFitFunc(nameFitFunc6,0,n6,5.,175.,point6);
//     if (grTh6)
//     {
//       grTh6->Draw("L");
//       grTh6->SetLineColor(kGreen);
//       grTh6->SetLineWidth(3);
//       grTh6->SetLineStyle(2);
//       lg->AddEntry(grTh6,leg6,"l");
//     }
//   }
//   
//   //
//   lg->Draw();
//   
//   //
//   return canv;
//   
// }
// 
