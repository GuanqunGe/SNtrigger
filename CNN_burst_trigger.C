#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include <TCanvas.h>
#include <TString.h>
#include <TFile.h>
#include <TF1.h>
#include <TH1D.h>
#include <TMath.h>
#include <TGraph.h>
#include <TRandom.h>

using namespace std;

TFile *f_out;
TFile *f1;
TFile *f_in_theory;
TFile *f_in_sim;

  //Your goal is to:
  //Compare multiplicity for F and S for different N, A
  //Find the combination of (multiplicity cut, N, A) that maximizes burst trigger efficiency and keeps burst trigger fake rate to below 1/month


int main(int argc, const char** argv){

  double signaleff = 0.8927; //1.00=100%
  double bkgdeff = 0.0056; 
  bool killframeswithmorethanoneint = false; //setting this to true assumes that frames with more than one SN neutrino interaction (we didn't train on these) are identified as background instead of signal
  double L = 10;//in kpc; distance to SN 

  int N = 20; //number of successive frames to integrate trigger primitives over
  //  int A = 200; //number of APA to integrate trigger primitives over; not yet implemented!!! Only 200 is possible.

  int nAPAs = 200;
 
 
  std::string seq(argv[1]); //char** is pointer to  pointer of char, so argv[1] is a pointer to char/string.  

  
  std::string filename=Form("/a/data/westside/guanqun/CNNburst/pdf/CNN_RESULT_%i_%0.2f_",N,L)  + seq +".pdf"; //file will be put into corresponding directory.
  std::string rootname = Form("/a/data/westside/guanqun/CNNburst/root/CNN_burst_%i_%0.2f_",N,L)  + seq + ".root";
  std::string root2 = Form("/a/data/westside/guanqun/CNNburst/root_Max/CNN_burst_%i_%0.2f_",N,L)  + seq + ".root";
  std::string filestart = filename + "["; //open the pdf file
  std::string fileend = filename + "]";   //close this pdf file

 
 
  //For a period of T=>10 seconds, generate a distribution of neutrinos vs. time from a SN burst at some baseline L. Baseline L defines overall normalization of time distribution to draw from.
  f_in_theory = new TFile("/a/share/amsterdam/gegq/bursttrigger_v1/SNTheoryDistributions.root","READ");
  f_in_sim = new TFile("/a/share/amsterdam/gegq/bursttrigger_v1/TimeProfile.root","READ");
  TF1 *f_EventsVSNDistance_10kt_100kpc = (TF1*)f_in_theory->Get("f_EventsVSNDistance_10kt_100kpc");
  TH1D *h_MarlTime = (TH1D*)f_in_sim->Get("h_MarlTime");//up to 10 sec, in bins of 100[ms]
  TH1D *h_MarlTimeZoom = (TH1D*)f_in_sim->Get("h_MarlTimeZoom");//up to 50 [msec], in bins of 0.25[ms] 
  
  double T = 10.00125;//in seconds; corresponds to 4445 frames of 2.25ms
  double tframe = 0.00225;//in seconds
  int nframes = int(T/tframe); //4445
  int zoombinsperframe = int(tframe*1000/h_MarlTimeZoom->GetBinWidth(1)); //2.25ms/0.25ms=9 bins
  double framesperunzoombin = h_MarlTime->GetBinWidth(1)/(tframe*1000.); //100ms/2.25ms=44.44frames per bin

  TH1D *h_SNEvents_at_L = new TH1D("h_SNEvents_at_L","SN Burst Events vs. Time for 10kton Module",nframes,0.,T);
  //first fill above histogram by merging bins from the zoomed histogram and by subdividing bins from the umerged one; filled histogram then got normalized.
  int inew,iold;
  for (int i=0; i<h_MarlTimeZoom->GetNbinsX(); i++){ //start from i = 0 in order to get the sum of h_MarlTimeZoom bins from bin 1 to bin 9.
    inew = double(i)/double(zoombinsperframe);
    h_SNEvents_at_L->SetBinContent(inew+1,h_SNEvents_at_L->GetBinContent(inew+1)+h_MarlTimeZoom->GetBinContent(i+1));
  }
  for (int i=inew; i<h_SNEvents_at_L->GetNbinsX(); i++){
    iold = (1000.*h_SNEvents_at_L->GetBinCenter(i+1))/h_MarlTime->GetBinWidth(1);
    if (iold!=0) 
      h_SNEvents_at_L->SetBinContent(i+1,h_MarlTime->GetBinContent(iold+1)/framesperunzoombin);
    else      
      h_SNEvents_at_L->SetBinContent(i+1,(h_MarlTime->GetBinContent(iold+1)/framesperunzoombin)/2); 
  }

  f_out = new TFile( rootname.c_str(), "RECREATE");
  f_out->cd();
  TCanvas* c = new TCanvas();
  c->Print(filestart.c_str());
  h_SNEvents_at_L->Write();
  h_SNEvents_at_L->Draw();
  c->Print(filename.c_str());
  


  //Next, scale h_SN_Events_at_L histogram according to L for SN (overall normalization
  //Total expected number of SN nu interactions should be:
  double n_SNEvents_tot = f_EventsVSNDistance_10kt_100kpc->Eval(L); //add randomization (poisson fluctuations) to this later, per bin
  double scalefactor = n_SNEvents_tot/h_SNEvents_at_L->Integral()/double(nAPAs); //events per APA(scale the total probability to 1 at the same time)
  h_SNEvents_at_L->Scale(scalefactor); //This is the number of events per APA
  //check:  
  //cout << "Expected number of SN interactions (1 APA): " << h_SNEvents_at_L->Integral() << endl;

  std::map<double,int>map_MaxtoFrequency_signal;//create map to store Maximum of hist and it's frequency for further analysis of M_cut
  std::map<double,int>map_MaxtoFrequency_bkgd;

  double Min_Max_signal=0.0; //keep the minimum of Maximum of signal histogram
  double Max_bkgd=0.0; // always keep the highest value of bkgd histogram
 
  int numT = 1040; //times of variation;
 
  for(int num = 0 ; num < numT  ; num++ ){   // number of total 10s-histograms in this code.
    TH1D *h_SNEvents_at_L_var[nAPAs];
    char hname[50];
    TRandom *rand = new TRandom(0);
    double poissonmean; 
   
    //Next, vary each bin according to poisson distribution; this is the fluctuated number of events per APA:
    for (int n=0; n<nAPAs; n++){
      sprintf(hname,"h_SNEvents_at_L_var_%i",n);
      h_SNEvents_at_L_var[n] = new TH1D(hname,"Random Throw of SN Burst Events at Given L",nframes,0.,T);
      for (int i=0; i<nframes; i++){
         poissonmean = h_SNEvents_at_L->GetBinContent(i+1);
         h_SNEvents_at_L_var[n]->SetBinContent(i+1,(double)rand->Poisson(poissonmean));
      }
    //check:
    //cout << "Poisson-fluctuated (bin-by-bin) number of SN interactions (APA " << n <<"): " << h_SNEvents_at_L_var[n]->Integral() << endl;
    }
  
   //Check and make sure that when the T is split into frames of 2.25ms, no frame contains more than one event; if so, output a warning
    int nframesmorethanoneint = 0;
    for (int n=0; n<nAPAs; n++){
      nframesmorethanoneint = 0;
      for (int i=0; i<h_SNEvents_at_L_var[n]->GetNbinsX(); i++){
        if (h_SNEvents_at_L_var[n]->GetBinContent(i+1)>1){
        	nframesmorethanoneint++;
        //	cout << "WARNING: More than one neutrino interactions in this frame (frame:" << i+1 << ", nint:" << h_SNEvents_at_L_var[n]->GetBinContent(i+1) << ")! Continuing anyway..." << endl;
      }
    }
      //if (nframesmorethanoneint>0) cout << "-------> Fraction of frames with more than one SN neutrino interaction: " << double(nframesmorethanoneint)*100./double(nframes) << "%" << endl;
    }
  
 
   //Now, apply the signal selection efficiency independently on each frame, to determine which frames are tagged as signal (S)
    TH1D *signal_tagged_frames[nAPAs];
    for (int n=0; n<nAPAs; n++){
      sprintf(hname,"signal_tagged_frames_%i",n);
      signal_tagged_frames[n] = new TH1D(hname,"SN burst frames selected based on efficiency",nframes,0.,double(nframes));
      for (int i=0; i<nframes; i++){
        if (h_SNEvents_at_L_var[n]->GetBinContent(i+1)>0){//if there's a SN interaction in this frame
	  if ((double)rand->Uniform(0,1)<=signaleff) 
	    signal_tagged_frames[n]->SetBinContent(i+1,1);//efficiency is applied as a probability of tagging particular frame as signal
	  if (h_SNEvents_at_L_var[n]->GetBinContent(i+1)>1){//special case: if there's more than one SN interaction in this frame
	    if (killframeswithmorethanoneint){//treat as background
	      if ((double)rand->Uniform(0,1)<=bkgdeff)
	        signal_tagged_frames[n]->SetBinContent(i+1,1);//very small probability of accepting
	      else
	        signal_tagged_frames[n]->SetBinContent(i+1,0);//otherwise signal frame is missed entirely
	    }
	  }
        }
        else { //empty frame in SN burst behaves as "bkgd"
	  if ((double)rand->Uniform(0,1)<=bkgdeff)
	    signal_tagged_frames[n]->SetBinContent(i+1,1);//efficiency is applied as a probability of tagging particular frame as signal, when empty
        }
      }
    }
  
    //For the same period T, generate a distribution of fake trigger frames, F, by applying the fake trigger effiency per frame
    TH1D *bkgd_tagged_frames[nAPAs];
    for (int n=0; n<nAPAs; n++){
      sprintf(hname,"bkgd_tagged_frames_%i",n);
      bkgd_tagged_frames[n] = new TH1D(hname,"Background frames selected based on efficiency",nframes,0.,double(nframes));
      for (int i=0; i<nframes; i++){
        if ((double)rand->Uniform(0,1)<=bkgdeff)
	  bkgd_tagged_frames[n]->SetBinContent(i+1,1);//background frames mis-tagged as signal
      }
    }

    //Sum the signal S and background B frame multiplicity over N successive frames and A APAs in the 10kton module
    //For now, only A=all 200 APAs is possible!!!
    TH1D S_multi(Form("S_multi_%i",num),"Signal Multiplicity",nframes-N+1,0,nframes-N+1);
    TH1D B_multi(Form("B_multi_%i",num),"Background Multiplicity",nframes-N+1,0,nframes-N+1);
    int i=0;
    while (i<nframes-N+1){
      for (int n=0; n<nAPAs; n++){
       for (int f=0; f<N; f++){
       	 S_multi.SetBinContent(i+1,S_multi.GetBinContent(i+1)+signal_tagged_frames[n]->GetBinContent(i+1+f));
	 B_multi.SetBinContent(i+1,B_multi.GetBinContent(i+1)+bkgd_tagged_frames[n]->GetBinContent(i+1+f));
       }
      }
      i++;
    }

   //Save the highest/loewst value of signal and bkgd 200APAs-summed multiplicity histogram
    double S_max = S_multi.GetBinContent(S_multi.GetMaximumBin());
    if( Min_Max_signal == 0.0) Min_Max_signal = S_max;
      else if( S_max < Min_Max_signal)  Min_Max_signal = S_max;
        
    if(map_MaxtoFrequency_signal.find(S_max) != map_MaxtoFrequency_signal.end())
       {map_MaxtoFrequency_signal[S_max]++ ;}
      else map_MaxtoFrequency_signal[S_max]=1;

    //now for bkgd
    double B_max = B_multi.GetBinContent(B_multi.GetMaximumBin());
    if (B_max > Max_bkgd) Max_bkgd = B_max;
    if(map_MaxtoFrequency_bkgd.find(B_max) != map_MaxtoFrequency_bkgd.end())
       {map_MaxtoFrequency_bkgd[B_max]++ ;}
      else map_MaxtoFrequency_bkgd[B_max] = 1;

    //Save some valuable information for checks
    TH1D h_SNEvents_at_L_var_total(Form("h_SNEvents_at_L_var_total_%i",num),"SN events for 10kton module vs time",nframes,0.,T);  
    for (int n=0; n<nAPAs; n++){
      for (int i=0; i<nframes; i++){
        h_SNEvents_at_L_var_total.SetBinContent(i+1,h_SNEvents_at_L_var_total.GetBinContent(i+1)+h_SNEvents_at_L_var[n]->GetBinContent(i+1));
      }
    }

    //write out totals
    f_out->cd();
    h_SNEvents_at_L_var_total.Write();
    S_multi.Write();
    S_multi.Draw();
    c->Print(filename.c_str());
    S_multi.GetXaxis()->SetRange(0,100);
    S_multi.Draw();
    c->Print(filename.c_str());
    B_multi.Write();
    B_multi.Draw();
    c->Print(filename.c_str());
   // h_SNEvents_at_L_var_total.Draw();
   // c->Print(filename.c_str());
   
    
    //delete the pointers, to avoid memory leak;
    delete rand;
    for(int n = 0; n < nAPAs; n ++ ){
       delete h_SNEvents_at_L_var[n];
       delete signal_tagged_frames[n];
       delete bkgd_tagged_frames[n];
    }
  }
    c->Print(fileend.c_str());

    f_out->Close();
   
/*
   std::cout << "********now start to print out maximum value of S_multi histogram and frequency******************" << std::endl;
  for(map<double,int>::iterator iter = map_MaxtoFrequency_signal.begin(); iter != map_MaxtoFrequency_signal.end(); iter++)     {std::cout << "Maximum of his is : " << iter->first  << ", frequency is " << iter->second << std::endl;}

   std::cout << "********now start to print out maximum value of B_multi histogram and frequency******************" << std
::endl;
  for(map<double,int>::iterator iter = map_MaxtoFrequency_bkgd.begin(); iter != map_MaxtoFrequency_bkgd.end(); iter++)     {std::cout << "Maximum of his is : " << iter->first  << ", frequency is " << iter->second << std::endl;}

  std::cout << "*******now start to print out extreme value***************"<<std::endl;
  std::cout << "the lowest of the maximum value of " << numT << "s signal histogram is " << Min_Max_signal <<std::endl;
  std::cout << "the highest of the maximum value of "<< numT << "s bkgd histogram is " << Max_bkgd << std::endl;
*/

  //directly print out highest value is not convenient, use histogram and add them.
  f1 = new TFile(root2.c_str(),"RECREATE");
  
  // draw histogram with maximum and frequency.
  TH1D* hMaxSignal = new TH1D(Form("hMaxSignal_%i_%i", N, (int)L),"hMaxSignal", 150 , 0.0, 150.0);
  TH1D* hMaxBkgd = new TH1D(Form("hMaxBkgd_%i_%i", N, (int)L), "hMaxBkgd", 100,0,100);
  for(map<double,int>::iterator iter = map_MaxtoFrequency_signal.begin(); iter != map_MaxtoFrequency_signal.end(); iter++) 
    {hMaxSignal->Fill(iter->first, (double)iter->second);}
  for(map<double,int>::iterator iter = map_MaxtoFrequency_bkgd.begin(); iter != map_MaxtoFrequency_bkgd.end(); iter++)
    {hMaxBkgd->Fill(iter->first,(double)iter->second);}

  f1->cd();
  hMaxSignal->Write();
  hMaxBkgd->Write();
  //f_out->Close();

  delete hMaxSignal;
  delete hMaxBkgd;

  f1->Close(); // close root file after delete pointers, otherwise it will crush.

  return 0;

}

