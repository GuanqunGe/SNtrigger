#include <algorithm>
#include <iostream>
#include <vector>
#include <fstream>

#include "TFile.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TH1D.h"
#include "THStack.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TStyle.h"
#include "TText.h"
#include "TROOT.h"

#include "colors.h"

using namespace std;
std::vector<int> SelectedConfig = {0,1,2,3,4,5};

int main()
{
  TString s_Filename = "GH_SNMC";
  ifstream inFile;
  inFile.open("Analyse_SNBurst_"+s_Filename+".txt");

  std::map<double,std::pair<double,double>> map_ConfigToEffAndBkgd;
  double Config;
  int Timewindow;
  double Eff, Bkgd, dummy;
  std::vector<double> vec_Config;
  
  while(inFile >> Config >> Timewindow >> Eff >> Bkgd >> dummy){ // time window???? config = cut ??
    std::cout << "Config " << Config
      << ", Time window " << Timewindow
      << ", Eff " << Eff
      << ", Bkgd " << Bkgd
      << std::endl;
   std::cout << "works fine" << std::endl;
    if(std::find(vec_Config.begin(), vec_Config.end(), Config) == vec_Config.end())
      vec_Config.push_back(Config);
    
    map_ConfigToEffAndBkgd[Config] = {Eff,Bkgd};
  }
  
  int nConfig = (int)vec_Config.size();
//  std::cout << "There are " << nConfig << " configs." << std::endl;
 // std::cout << "Size of map is " << map_ConfigToEffAndBkgd.size() << std::endl; 
 
  gROOT->ForceStyle();
  std::vector<int> vec_Colors = getColors(2);// initialize vec_colors with 36 elements
  TFile *f_SNTheoryDistributions = new TFile("SNTheoryDistributions.root","READ");
  TH1D  *h_SNProbabilityVDistance = (TH1D*)f_SNTheoryDistributions->Get("h_SNProbabilityVDistance_LMC");
  h_SNProbabilityVDistance->SetLineWidth(3);
  h_SNProbabilityVDistance->SetLineColor(46);


  TFile *f_Input  = new TFile("Analyse_SNBurst_GH_SNMC.root", "READ");
  std::map<double,TH1D*>   map_h_FakeRateVNClusters;
  std::map<double,TH1D*>   map_h_FakeRateVNClustersLow;
  std::map<double,TH1D*>   map_h_EfficiencyVEvents;
  std::map<double,TH1D*>   map_h_EfficiencyVDistance;
  std::map<double,TH1D*>   map_h_EffGalaxy;
  std::map<double,TGraph*> map_g_ROC;

  std::map<double,std::pair<double,double>>::iterator it_Config;
  int globalIt=0;
  for(auto const& it_Config : map_ConfigToEffAndBkgd){
    int color = vec_Colors.at(globalIt % vec_Colors.size());

    double config     = it_Config.first;
    
    std::cout << "Dealing with Config " << config << std::endl;
    TString s_FakeRateVNClusters     = Form("h_FakeRateVNClusters_Cut%0.5f_TW%i",    config, Timewindow);
    TString s_FakeRateVNClustersLow  = Form("h_FakeRateVNClustersLow_Config%0.5f_TW%i", config, Timewindow);
    TString s_EfficiencyVEvents      = Form("h_EfficiencyVEvents_Cut%0.5f_TW%i",     config, Timewindow);
    TString s_EfficiencyVDistance    = Form("h_EfficiencyVDistance_Cut%0.5f_TW%i",   config, Timewindow);
    TString s_EffGalaxy              = Form("h_NeighbourhoodEffiency_Cut%0.5f_TW%i", config, Timewindow);
    TString s_ROC                    = Form("g_ROC_Cut%0.5f_TW%i",                   config, Timewindow);
    
    TH1D *h_FakeRateVNClusters    = (TH1D*)f_Input->Get(s_FakeRateVNClusters); 
    TH1D *h_FakeRateVNClustersLow = (TH1D*)f_Input->Get(s_FakeRateVNClustersLow);
    
    TH1D *h_EfficiencyVEvents = (TH1D*)f_Input->Get(s_EfficiencyVEvents);
    if(!h_EfficiencyVEvents){
      std::cout << "Erasing config " << config << " as there is no EfficiencyVEvents plot" << std::endl;
      map_ConfigToEffAndBkgd.erase(config);
      continue;
    }

    int   maxBin_Events = h_EfficiencyVEvents->GetMaximumBin();
    for(int i = maxBin_Events; i < h_EfficiencyVEvents->GetSize()-1 ; i++)
      h_EfficiencyVEvents->SetBinContent(i,1);



    TH1D *h_EfficiencyVDistance = (TH1D*)f_Input->Get(s_EfficiencyVDistance);
    int   maxBin_Distance = h_EfficiencyVDistance->GetMaximumBin();
    for(int i = maxBin_Distance; i > 0; i--)
      h_EfficiencyVDistance->SetBinContent(i,1);


    TH1D   *h_EffGalaxy           = (TH1D*)f_Input->Get(s_EffGalaxy);
    TGraph *g_ROC                 = (TGraph*)f_Input->Get(s_ROC);
    h_FakeRateVNClusters   ->SetLineColor  (color);
    h_FakeRateVNClustersLow->SetLineColor  (color);
    h_FakeRateVNClusters   ->SetMarkerColor(color);
    h_FakeRateVNClustersLow->SetMarkerColor(color);
    h_EfficiencyVEvents    ->SetLineColor  (color);
    h_EfficiencyVDistance  ->SetLineColor  (color);
    h_EffGalaxy            ->SetLineColor  (color);
    g_ROC                  ->SetMarkerColor(color);
    g_ROC                  ->SetLineColor  (color);

    g_ROC->SetMarkerStyle(3);

    map_h_FakeRateVNClusters   [config] = h_FakeRateVNClusters;
    map_h_FakeRateVNClustersLow[config] = h_FakeRateVNClustersLow;
    map_h_EfficiencyVEvents    [config] = h_EfficiencyVEvents;
    map_h_EfficiencyVDistance  [config] = h_EfficiencyVDistance;
    map_h_EffGalaxy            [config] = h_EffGalaxy;
    map_g_ROC                  [config] = g_ROC;
    globalIt++;    
  }
  std::cout << "number of config is "<<  globalIt << std::endl;
  TCanvas *c_Global = new TCanvas();
  c_Global->Print("alexResults.pdf[");
  std::string legHeader = "Individual Marley Eff & 10kt Bkgd Rate";
  std::string legEntryFormat = "Config: %.5f - Eff: %.3f & Bkgd rate: %.2f Hz";
  
  THStack *stk_FakeRateVNClusters = new THStack("stk_FakeRateVNClusters", "Number of Clusters in Time Window Required to Trigger vs. Trigger Rate");
  TLegend *leg_FakeRateVNClusters = new TLegend(0.05, 0.05, 0.95, 0.95);
  //leg_FakeRateVNClusters->SetTextSize(0.023);
  leg_FakeRateVNClusters->SetHeader(legHeader.c_str());

  for(auto const& it_Config : map_ConfigToEffAndBkgd){
    stk_FakeRateVNClusters->Add(map_h_FakeRateVNClusters     [it_Config.first]);
    stk_FakeRateVNClusters->Add(map_h_FakeRateVNClustersLow  [it_Config.first]);
    leg_FakeRateVNClusters->AddEntry(map_h_FakeRateVNClusters[it_Config.first],
                                     Form(legEntryFormat.c_str(), it_Config.first, it_Config.second.first, it_Config.second.second), "L");
  }
  //  std::cout << "size of the config map is" << map_ConfigToEffAndBkgd.size() << std::endl;
  leg_FakeRateVNClusters->Draw();
  c_Global->Print("alexResults.pdf");
  c_Global->SetLogy();
  stk_FakeRateVNClusters->SetMinimum(1e-9);
  stk_FakeRateVNClusters->Draw("NOSTACK C");  
  stk_FakeRateVNClusters->GetXaxis()->SetTitle("Number of Clusters/Time Window");
  stk_FakeRateVNClusters->GetYaxis()->SetTitle("Trigger Rate, (Hz)");
  stk_FakeRateVNClusters->GetXaxis()->SetLimits(0,map_h_FakeRateVNClusters.begin()->second->GetXaxis()->GetXmax());
  std::cout << "this quantity is "<< map_h_FakeRateVNClusters.begin()->second->GetXaxis()->GetXmax() << std::endl;
  double range = map_h_FakeRateVNClusters.begin()->second->GetXaxis()->GetXmax() - map_h_FakeRateVNClusters.begin()->second->GetXaxis()->GetXmin();
  TText *t_perMonth  = new TText(map_h_FakeRateVNClusters.begin()->second->GetXaxis()->GetXmax()-0.15*range, 2e-7,    "1/Month");
  TText *t_perWeek   = new TText(map_h_FakeRateVNClusters.begin()->second->GetXaxis()->GetXmax()-0.15*range, 1.65e-6, "1/Week");
  TText *t_perDay    = new TText(map_h_FakeRateVNClusters.begin()->second->GetXaxis()->GetXmax()-0.15*range, 1.16e-5, "1/Day");
  TLine *l_perMonth  = new TLine(0, 4.13e-7, map_h_FakeRateVNClusters.begin()->second->GetXaxis()->GetXmax(), 4.13e-7);
  TLine *l_perWeek   = new TLine(0, 1.65e-6, map_h_FakeRateVNClusters.begin()->second->GetXaxis()->GetXmax(), 1.65e-6);
  TLine *l_perDay    = new TLine(0, 1.16e-5, map_h_FakeRateVNClusters.begin()->second->GetXaxis()->GetXmax(), 1.16e-5);
//  std::cout << "check point 3-4-1" << std::endl;
  l_perMonth->SetLineColor(4);
  l_perWeek ->SetLineColor(4);
  l_perDay  ->SetLineColor(4);
  l_perMonth->SetLineWidth(3);
  l_perWeek ->SetLineWidth(3);
  l_perDay  ->SetLineWidth(3);
  t_perMonth->Draw();
  t_perWeek ->Draw();
  t_perDay  ->Draw();
  l_perMonth->Draw();
  l_perWeek ->Draw();
  l_perDay  ->Draw();
  gPad->RedrawAxis(); 
  c_Global->Print("alexResults.pdf");
  THStack *stk_EfficiencyVEvents = new THStack("stk_EfficiencyVEvents", "Burst Efficiency vs. Number of Events in SN Burst, Fake Trigger Rate: 1/Month");
  TLegend *leg_EfficiencyVEvents = new TLegend(0.53, 0.75, 0.85, 0.85);
  leg_EfficiencyVEvents->SetTextSize(0.023);
  
  leg_EfficiencyVEvents->SetHeader(legHeader.c_str());
  for(auto const& it_Config : map_ConfigToEffAndBkgd){
    stk_EfficiencyVEvents->Add(map_h_EfficiencyVEvents     [it_Config.first]);
    leg_EfficiencyVEvents->AddEntry(map_h_EfficiencyVEvents[it_Config.first], Form(legEntryFormat.c_str(), it_Config.first, it_Config.second.first, it_Config.second.second), "L");
  }
  c_Global->Clear();
  c_Global->SetLogy(false);
  c_Global->Draw();
  c_Global->SetLogx();
  //leg_EfficiencyVEvents->Draw(); // if you draw it here, legend box won't show up.
  stk_EfficiencyVEvents->Draw("NOSTACK");
  stk_EfficiencyVEvents->GetXaxis()->SetTitle("Number of Events in SN Burst");
  stk_EfficiencyVEvents->GetYaxis()->SetTitle("Burst Efficiency");
  gPad->RedrawAxis();
 // leg_EfficiencyVEvents->Draw();
  c_Global->Print("alexResults.pdf");

  THStack *stk_EfficiencyVDistance = new THStack("stk_EfficiencyVDistance", "Burst Efficiency vs. Distance to SN, Fake Trigger Rate: 1/Month");
  TLegend *leg_EfficiencyVDistance = new TLegend(0.15, 0.15, 0.48, 0.35);
  leg_EfficiencyVDistance->SetTextSize(0.023);

  leg_EfficiencyVDistance->SetHeader(legHeader.c_str());
  for(auto const& it_Config : map_ConfigToEffAndBkgd){
    stk_EfficiencyVDistance->Add(map_h_EfficiencyVDistance[it_Config.first]);
    leg_EfficiencyVDistance->AddEntry(map_h_EfficiencyVDistance[it_Config.first], 
                                      Form(legEntryFormat.c_str(), it_Config.first, it_Config.second.first, it_Config.second.second), "L");
  }
  c_Global->Clear();
  c_Global->Draw();
  c_Global->SetLogx();
  leg_EfficiencyVDistance->Draw();
  stk_EfficiencyVDistance->Draw("NOSTACK");
  stk_EfficiencyVDistance->GetXaxis()->SetTitle("SN Distance, (kpc)");
  stk_EfficiencyVDistance->GetYaxis()->SetTitle("Burst Efficiency");
  gPad->RedrawAxis();
  c_Global->Print("alexResults.pdf");


 // gStyle->SetOptStat(0);
  THStack *stk_EffGalaxy = new THStack("stk_EffGalaxy", "Galactic Neighbourhood Coverage, Fake Trigger Rate 1/Month");
  TLegend *leg_EffGalaxy = new TLegend(0.15, 0.15, 0.48, 0.35);
  leg_EffGalaxy->SetTextSize(0.023);

  leg_EffGalaxy->SetHeader(legHeader.c_str());
  stk_EffGalaxy->Add(h_SNProbabilityVDistance);
  leg_EffGalaxy->AddEntry(h_SNProbabilityVDistance, "SN Probability", "L");

  for(auto const& it_Config : map_ConfigToEffAndBkgd){
    stk_EffGalaxy->Add(map_h_EffGalaxy[it_Config.first]);
    leg_EffGalaxy->AddEntry(map_h_EffGalaxy[it_Config.first],
                            Form(legEntryFormat.c_str(), it_Config.first, it_Config.second.first, it_Config.second.second), "L");
  }
  c_Global->Clear();
  c_Global->Draw();
  c_Global->SetLogx(false);
  c_Global->SetLogy();
  leg_EffGalaxy->Draw();
  stk_EffGalaxy->Draw("NOSTACK");
  stk_EffGalaxy->GetXaxis()->SetTitle("SN Distance, (kpc)");
  stk_EffGalaxy->GetYaxis()->SetTitle("Burst Efficiency x SN Probability");
  stk_EffGalaxy->Draw("NOSTACK");
  gPad->RedrawAxis();
  //leg_EffGalaxy->Draw();
  c_Global->Print("alexResults.pdf");

  //TMultiGraph *mulgraph = new TMultiGraph;
//  mulgraph->SetTitle("Fake Trigger Rate vs. Galactic Neighbourhood Coverage;Galactic Neighbourhood Coverage;Fake Trigger Rate, (Hz)");
  TLegend *leg_ROC    = new TLegend(0.15, 0.68, 0.48, 0.88);
  leg_ROC->SetTextSize(0.023);
  leg_ROC->SetHeader(legHeader.c_str());
  c_Global->Clear();
  c_Global->Draw();
  c_Global->SetLogx(false);
  c_Global->SetLogy();
  double minX=0.6, maxX=1.0;
  double minY=10e-15, maxY=10e5;

  map_g_ROC.begin()->second->GetXaxis()->SetLimits(minX, maxX);
  map_g_ROC.begin()->second->SetMaximum(maxY);
  map_g_ROC.begin()->second->SetMinimum(minY);
  map_g_ROC.begin()->second->SetTitle("Fake Trigger Rate vs. Galactic Neighbourhood Coverage");
  map_g_ROC.begin()->second->GetXaxis()->SetTitle("Galactic Neighbourhood Coverage");
  map_g_ROC.begin()->second->GetYaxis()->SetTitle("Fake Trigger Rate, (Hz)");
  //map_g_ROC.begin()->second->SetMaximum(10);
  map_g_ROC.begin()->second->Draw("AP");

  for(auto const& it_Config : map_ConfigToEffAndBkgd){
    leg_ROC->AddEntry(map_g_ROC[it_Config.first], Form(legEntryFormat.c_str(), it_Config.first, it_Config.second.first, it_Config.second.second), "P");
   // mulgraph->Add(map_g_ROC[it_Config.first]);
    map_g_ROC[it_Config.first]->Draw("P");
  }
 // mulgraph->SetMinimum(minY);
//  mulgraph->SetMaximum(maxY);
//  mulgraph->GetXaxis()->SetLimits(minX, maxX);
//  mulgraph->Draw("AP");
  TLine *l_perMonth_2 = new TLine(minX, 4.13e-7, maxX, 4.13e-7);
  l_perMonth_2->SetLineColor(1);
  l_perMonth_2->SetLineWidth(3);
  TText *t_perMonth_2 = new TText(minX+0.015, 8e-7, "1/Month");
  gPad->RedrawAxis();
  //leg_ROC->Draw();
  l_perMonth_2->Draw();
  t_perMonth_2->Draw();
  c_Global->Print("alexResults.pdf");
  c_Global->Print("alexResults.pdf]");

  return 0;
}
