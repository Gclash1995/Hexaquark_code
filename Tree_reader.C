#include <cstdlib>
#include <iostream>
#include <chrono>
#include <TFile.h>
#include <TTree.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TLorentzVector.h>
#include <TH1.h>
#include <TH2.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TBenchmark.h>
#include <vector>
{
  gROOT->ProcessLine(".L ./Loader.C+");
  TFile *f = new TFile("skim4_Tree1.root");
  TTree *t1 = (TTree*)f->Get("skim4_Tree");

  vector<TLorentzVector> *v_p4;


  TLorentzVector *readbeam;
  TLorentzVector *readtarget;


  vector<TLorentzVector> *v_vertex;

  vector<double> *v_beta;

  Double_t start_time;
  vector<double> *energy;
  vector<double> *charge;
  vector<double> *PID;
  vector<double> *chi2PID;

  Int_t readchargetracks;
  Int_t readprotonno;
  Int_t readpipno;
  Int_t readpimno;
  Int_t readelno;


  t1->SetBranchAddress("p4",&v_p4);

  t1->SetBranchAddress("vertex",&v_vertex);

  t1->SetBranchAddress("beta",&v_beta);

  t1->SetBranchAddress("beam",&readbeam);
  t1->SetBranchAddress("target",&readtarget);

  t1->SetBranchAddress("start_time",&start_time);
  t1->SetBranchAddress("energy",&energy);
  t1->SetBranchAddress("charge",&charge);
  t1->SetBranchAddress("PID",&PID);
  t1->SetBranchAddress("chi2PID",&chi2PID);
  t1->SetBranchAddress("chargetracks",&readchargetracks);
  t1->SetBranchAddress("protonno",&readprotonno);
  t1->SetBranchAddress("pipno",&readpipno);
  t1->SetBranchAddress("pimno",&readpimno);
  t1->SetBranchAddress("elno",&readelno);

  auto* hmiss=new TH1F("hmiss","Missing Mass",200,0,3);
  hmiss->GetXaxis()->SetTitle("MM(e #pi^{+} #pi^{-}) [GeV/c^{2}]");
  hmiss->GetYaxis()->SetTitle("Counts");
  hmiss->GetXaxis()->SetLabelSize(0.05);
  hmiss->GetYaxis()->SetLabelSize(0.05);

  auto* hmK0=new TH1F("hmK0","Invariant Mass",200,0,2);
  hmK0->GetXaxis()->SetTitle("M(#pi^{+} #pi^{-}) [GeV/c^{2}]");
  hmK0->GetYaxis()->SetTitle("Counts");
  hmK0->GetXaxis()->SetLabelSize(0.05);
  hmK0->GetYaxis()->SetLabelSize(0.05);

  auto* hbeta_pip=new TH2F("hbeta_pip","#Delta#beta for #pi^{+} (post PID)",300,0,10,200,-0.02,0.02);
  hbeta_pip->GetXaxis()->SetTitle("P [GeV/c]");
  hbeta_pip->GetYaxis()->SetTitle("#Delta#beta");
  hbeta_pip->GetXaxis()->SetLabelSize(0.05);
  hbeta_pip->GetYaxis()->SetLabelSize(0.05);

  auto* hbeta_pipfit=new TH2F("hbeta_pipfit","#Delta#beta for #pi^{+} (post PID)",300,0,10,200,-0.02,0.02);
  hbeta_pipfit->GetXaxis()->SetTitle("P [GeV/c]");
  hbeta_pipfit->GetYaxis()->SetTitle("#Delta#beta");
  hbeta_pipfit->GetXaxis()->SetLabelSize(0.05);
  hbeta_pipfit->GetYaxis()->SetLabelSize(0.05);

  auto* hbeta_pim=new TH2F("hbeta_pim","#Delta#beta for #pi^{-} (post PID)",300,0,10,200,-0.02,0.02);
  hbeta_pim->GetXaxis()->SetTitle("P [GeV/c]");
  hbeta_pim->GetYaxis()->SetTitle("#Delta#beta");
  hbeta_pim->GetXaxis()->SetLabelSize(0.05);
  hbeta_pim->GetYaxis()->SetLabelSize(0.05);

  auto* hbeta_pimfit=new TH2F("hbeta_pimfit","#Delta#beta for #pi^{-} (post PID)",300,0,10,200,-0.02,0.02);
  hbeta_pimfit->GetXaxis()->SetTitle("P [GeV/c]");
  hbeta_pimfit->GetYaxis()->SetTitle("#Delta#beta");
  hbeta_pimfit->GetXaxis()->SetLabelSize(0.05);
  hbeta_pimfit->GetYaxis()->SetLabelSize(0.05);

  auto* hbeta_pr=new TH2F("hbeta_pr","#Delta#beta for p (post PID)",300,0,10,200,-0.02,0.02);
  hbeta_pr->GetXaxis()->SetTitle("P [GeV/c]");
  hbeta_pr->GetYaxis()->SetTitle("#Delta#beta");
  hbeta_pr->GetXaxis()->SetLabelSize(0.05);
  hbeta_pr->GetYaxis()->SetLabelSize(0.05);

  auto* hbeta_prfit=new TH2F("hbeta_prfit","#Delta#beta for p (post PID)",300,0,10,200,-0.02,0.02);
  hbeta_prfit->GetXaxis()->SetTitle("P [GeV/c]");
  hbeta_prfit->GetYaxis()->SetTitle("#Delta#beta");
  hbeta_prfit->GetXaxis()->SetLabelSize(0.05);
  hbeta_prfit->GetYaxis()->SetLabelSize(0.05);

//cout<<readtarget->M()<<endl;
  TLorentzVector miss;
  TLorentzVector K0;

  Double_t beta_tof_pip;
  Double_t P_pip;
  Double_t beta_calc_pip;
  Double_t delta_beta_pip;

  Double_t beta_tof_pim;
  Double_t P_pim;
  Double_t beta_calc_pim;
  Double_t delta_beta_pim;

  Double_t beta_tof_pr;
  Double_t P_pr;
  Double_t beta_calc_pr;
  Double_t delta_beta_pr;

  TLorentzVector el;
  TLorentzVector pip;
  TLorentzVector pim;
  TLorentzVector pr;

  TFile fileOutput1("Cleaner_pions_and_proton.root","recreate");

 //
  Long64_t nentries = t1->GetEntries();
 // Long64_t nentries = 1000000;
  for(Long64_t i=0; i<nentries;i++){
    t1->GetEvent(i);

      Int_t Np = v_p4->size();
      for(Int_t j=0; j<Np; j++){
        if(v_p4->at(j).M()<0.1){
          el.SetXYZM(v_p4->at(j).Px(),v_p4->at(j).Py(),v_p4->at(j).Pz(),v_p4->at(j).M());
        }
        else if(v_p4->at(j).M()<0.5 && charge->at(j)>0){
          pip.SetXYZM(v_p4->at(j).Px(),v_p4->at(j).Py(),v_p4->at(j).Pz(),v_p4->at(j).M());
          beta_tof_pip = v_beta->at(j);
        }
        else if(v_p4->at(j).M()<0.5 && charge->at(j)<0){
          pim.SetXYZM(v_p4->at(j).Px(),v_p4->at(j).Py(),v_p4->at(j).Pz(),v_p4->at(j).M());
          beta_tof_pim = v_beta->at(j);
        }
        else if(v_p4->at(j).M()>0.9){
          pr.SetXYZM(v_p4->at(j).Px(),v_p4->at(j).Py(),v_p4->at(j).Pz(),v_p4->at(j).M());
          beta_tof_pr = v_beta->at(j);
        }

      }
      miss = (TLorentzVector)*readbeam + (TLorentzVector)*readtarget - el - pim - pip;
      K0 = pim + pip;

      // cout<<miss.M()<<"   "<<abs(delta_beta_pip)<<"   "<<abs(delta_beta_pim)<<"   "<<abs(delta_beta_pim)<<"   "<<P_pip<<"   "<<P_pim<<"   "<<P_pr<<endl;
       //cout<<abs(delta_beta_pip)<<endl;
       //cout<<abs(delta_beta_pip)<<endl;
            // cout<<pip.Py()<<endl;
      // cout<<pim.Pz()<<endl;
      // cout<<pr.Pz()<<endl;

//cout<<abs(delta_beta_pr)<<endl;
      P_pip = sqrt((pow(pip.Px(),2))+(pow(pip.Py(),2))+(pow(pip.Pz(),2)));
      beta_calc_pip = P_pip/(sqrt((pow(P_pip,2))+(pow(pip.M(),2))));
      delta_beta_pip = beta_calc_pip-beta_tof_pip;
      //cout<<beta_tof_pip<<" "<<delta_beta_pip<<endl;

      P_pim = sqrt((pow(pim.Px(),2))+(pow(pim.Py(),2))+(pow(pim.Pz(),2)));
      beta_calc_pim = P_pim/(sqrt((pow(P_pim,2))+(pow(pim.M(),2))));
      delta_beta_pim = beta_calc_pim-beta_tof_pim;

      P_pr = sqrt((pow(pr.Px(),2))+(pow(pr.Py(),2))+(pow(pr.Pz(),2)));
      beta_calc_pr = P_pr/(sqrt((pow(P_pr,2))+(pow(pr.M(),2))));
      delta_beta_pr = beta_calc_pr-beta_tof_pr;

      hmiss->Fill(miss.M());

      if(miss.M()>0.8 && miss.M()<1.1 && abs(delta_beta_pip) < 0.02 && abs(delta_beta_pr) <0.02 &&
         abs(delta_beta_pim) < 0.02 && P_pip>0.2 && P_pim>0.2 && P_pr>0.2){
            //      cout<<miss.M()<<"   "<<beamA.M()<<"   "<<targetA.M()<<"   "<<el.M()<<endl;

           hmK0->Fill(K0.M());
           hbeta_pip->Fill(P_pip,delta_beta_pip);
           hbeta_pim->Fill(P_pim,delta_beta_pim);
           hbeta_pr->Fill(P_pr,delta_beta_pr);

           hbeta_pipfit->Fill(P_pip,delta_beta_pip);
           hbeta_pimfit->Fill(P_pim,delta_beta_pim);
           hbeta_prfit->Fill(P_pr,delta_beta_pr);
         }

  }

hbeta_pip->FitSlicesY();
// the mean and sigma
TH1D *hbeta_pip_1 = (TH1D*)gDirectory->Get("hbeta_pip_1");
TH1D *hbeta_pip_2 = (TH1D*)gDirectory->Get("hbeta_pip_2");

hbeta_pim->FitSlicesY();
// the mean and sigma
TH1D *hbeta_pim_1 = (TH1D*)gDirectory->Get("hbeta_pim_1");
TH1D *hbeta_pim_2 = (TH1D*)gDirectory->Get("hbeta_pim_2");

hbeta_pr->FitSlicesY();
// the mean and sigma
TH1D *hbeta_pr_1 = (TH1D*)gDirectory->Get("hbeta_pr_1");
TH1D *hbeta_pr_2 = (TH1D*)gDirectory->Get("hbeta_pr_2");

// TCanvas *can1=new TCanvas("can1","My Plot", 600, 600);
// can1->Divide(1,1);
// can1->cd(1);
// hmK0->Draw();
//
// TCanvas *can2=new TCanvas("can2","My Plot", 600, 600);
// can2->Divide(1,1);
// can2->cd(1);
// hbeta_pip->Draw("colz");
//
// TCanvas *can3=new TCanvas("can3","My Plot", 600, 600);
// can3->Divide(1,1);
// can3->cd(1);
// hbeta_pim->Draw("colz");
//
// TCanvas *can4=new TCanvas("can4","My Plot", 600, 600);
// can4->Divide(1,1);
// can4->cd(1);
// hbeta_pr->Draw("colz");
//
// TCanvas *can5=new TCanvas("can5","My Plot", 600, 600);
// can5->Divide(1,1);
// can5->cd(1);
// hmiss->Draw();

fileOutput1.Write();

}
