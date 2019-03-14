



void CountBackgroundEvents(){
  TFile* f = new TFile("/uboone/data/users/gge/copy/SNAna_EvNumber.root","READ");
  TTree* t = (TTree*)f->Get("snanagaushit/SNSimTree");
  TCanvas c= new TCanvas("c1","c1", 600,900);
  TH1D* histo = new TH1D("n_gen_background_evt",";n generated background events;n larsoft events", 10000,0,10000);
  
  t->Project("n_gen_background_evt","TotGen_APA+TotGen_CPA+TotGen_Ar39+TotGen_Neut+TotGen_Kryp+TotGen_Plon+TotGen_Rdon+TotGen_Ar42");
  histo->Draw();
  std::cout << "Number of background events: " << histo->GetMean() << " +/- " << histo->GetStdDev() << std::endl;
  histo->Reset();
  t->Project("n_gen_background_evt","TotGen_Ar39");
  std::cout << "Number of Ar39 background events: " << histo->GetMean() << " +/- " << histo->GetStdDev() << std::endl;
  
}
