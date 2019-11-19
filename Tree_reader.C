

{
  gROOT->ProcessLine(".L ./Loader.C+");
  TFile *f = new TFile("skim4_5201_Tree1.root");
  TTree *t1 = (TTree*)f->Get("skim4_5201_Tree");

  vector<TLorentzVector> *v_p4;
  //vector<TLorentzVector> *v_vertex;
  vector<double> *beta;
  Double_t start_time;
  vector<double> *energy;
  vector<double> *P;
  vector<double> *charge;
  vector<double> *chi2PID;
  Int_t readchargetracks;
  Int_t readprotonno;
  Int_t readkaonpno;
  Int_t readpipno;
  Int_t readpimno;

  t1->SetBranchAddress("p4",&v_p4);
  t1->SetBranchAddress("beta",&beta);
  t1->SetBranchAddress("start_time",&start_time);
  t1->SetBranchAddress("energy",&energy);
  t1->SetBranchAddress("P",&P);
  t1->SetBranchAddress("charge",&charge);
  t1->SetBranchAddress("chi2PID",&chi2PID);
  t1->SetBranchAddress("chargetracks",&readchargetracks);
  t1->SetBranchAddress("protonno",&readprotonno);
  t1->SetBranchAddress("kaonpno",&readkaonpno);
  t1->SetBranchAddress("pipno",&readpipno);
  t1->SetBranchAddress("pimno",&readpimno);


  auto* hpimno=new TH1F("hpimno","pimno",200,0,15);
  auto* hprotonno=new TH1F("hprotonno","protonno",200,0,15);


  Long64_t nentries = t1->GetEntries();

  for(Long64_t i=0; i<nentries;i++){
    t1->GetEvent(i);
    hpimno->Fill(readpimno);
    hprotonno->Fill(readprotonno);

  }

TCanvas *can1=new TCanvas("can1","My Plot", 600, 600);
can1->Divide(2,1);
can1->cd(1);
hpimno->Draw();
can1->cd(2);
hprotonno->Draw();

}
