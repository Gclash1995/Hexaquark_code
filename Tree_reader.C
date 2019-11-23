

{
  gROOT->ProcessLine(".L ./Loader.C+");
  TFile *f = new TFile("skim4_Tree1.root");
  TTree *t1 = (TTree*)f->Get("skim4_Tree");

  vector<TLorentzVector> *v_p4;
  vector<TLorentzVector> *v_vertex;
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
  Int_t readelno;
  TLorentzVector readbeam;
  TLorentzVector readtarget;

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
  t1->SetBranchAddress("target",&readtarget);
  t1->SetBranchAddress("beam",&readbeam);

  TLorentzVector miss;
  TLorentzVector miss2;
  TLorentzVector miss3;

  vector<TLorentzVector> v_miss;

  auto* hbeta_pip=new TH2F("hbeta_pip","hbeta_pip",100,0,3,200,-0.02,0.02);
  auto* hbeta_pim=new TH2F("hbeta_pim","hbeta_pim",100,0,3,200,-0.02,0.02);
  auto* hbeta_pr=new TH2F("hbeta_pr","hbeta_pr",100,0,3,200,-0.02,0.02);

  auto* hhist=new TH1F("hhist","hhist",100,0,1);



  Long64_t nentries = t1->GetEntries();

  for(Long64_t i=0; i<nentries;i++){
    t1->GetEvent(i);
    if(&readpipno==1 && &readpimno==1 && &readprotonno==1 && &readelno==1){
      for(Int_t j=0; j<v_p4->size();j++){
        miss = beam + target
        if(v_p4->at(j).M()>0.2){
          beta_calc_pr=P->at(j)/(sqrt((pow(P->at(j),2))+(pow((v_p4->at(j).M()),2))));
          beta_tof_pr=beta->at(j);

        }
    }

  }
  }

TCanvas *can1=new TCanvas("can1","My Plot", 600, 600);
can1->Divide(2,1);
can1->cd(1);
hhist->Draw();


}
