
#include <iostream>
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

  vector<TLorentzVector> *v_p4_el;
  vector<TLorentzVector> *v_p4_pip;
  vector<TLorentzVector> *v_p4_pim;
  vector<TLorentzVector> *v_p4_pr;

  TLorentzVector *readbeam=NULL;
  TLorentzVector *readtarget=NULL;

  vector<TLorentzVector> *v_vertex;

  vector<double> *v_beta_el;
  vector<double> *v_beta_pip;
  vector<double> *v_beta_pim;
  vector<double> *v_beta_pr;

  Double_t start_time;
  vector<double> *energy;
  vector<double> *P;
  vector<double> *charge;
  vector<double> *PID;
  vector<double> *chi2PID;
  Int_t readchargetracks;
  Int_t readprotonno;
  Int_t readpipno;
  Int_t readpimno;
  Int_t readelno;


  t1->SetBranchAddress("p4_el",&v_p4_el);
  t1->SetBranchAddress("p4_pip",&v_p4_pip);
  t1->SetBranchAddress("p4_pim",&v_p4_pim);
  t1->SetBranchAddress("p4_pr",&v_p4_pr);

  t1->SetBranchAddress("vertex",&v_vertex);

  t1->SetBranchAddress("beta_el",&v_beta_el);
  t1->SetBranchAddress("beta_pip",&v_beta_pip);
  t1->SetBranchAddress("beta_pim",&v_beta_pim);
  t1->SetBranchAddress("beta_pr",&v_beta_pr);

  t1->SetBranchAddress("beam",&readbeam);
  t1->SetBranchAddress("target",&readtarget);

  t1->SetBranchAddress("start_time",&start_time);
  t1->SetBranchAddress("energy",&energy);
  t1->SetBranchAddress("P",&P);
  t1->SetBranchAddress("charge",&charge);
  t1->SetBranchAddress("PID",&PID);
  t1->SetBranchAddress("chi2PID",&chi2PID);
  t1->SetBranchAddress("chargetracks",&readchargetracks);
  t1->SetBranchAddress("protonno",&readprotonno);
  t1->SetBranchAddress("pipno",&readpipno);
  t1->SetBranchAddress("pimno",&readpimno);
  t1->SetBranchAddress("pimno",&readelno);

  auto* hmiss=new TH1F("hmiss","Missing Mass",200,0,3);
  hmiss->GetXaxis()->SetTitle("MM(e #pi^{+} #pi^{-}) [GeV/c^{2}]");
  hmiss->GetYaxis()->SetTitle("Counts");
  hmiss->GetXaxis()->SetLabelSize(0.05);
  hmiss->GetYaxis()->SetLabelSize(0.05);

  auto* hmK0=new TH1F("hmK0","Invariant Mass",200,0,3);
  hmK0->GetXaxis()->SetTitle("M(#pi^{+} #pi^{-}) [GeV/c^{2}]");
  hmK0->GetYaxis()->SetTitle("Counts");
  hmK0->GetXaxis()->SetLabelSize(0.05);
  hmK0->GetYaxis()->SetLabelSize(0.05);

  auto* hbeta_pip=new TH2F("hbeta_pip","#Delta#beta for #pi^{+} (post PID)",300,0,10,200,-0.02,0.02);
  hbeta_pip->GetXaxis()->SetTitle("P [GeV/c]");
  hbeta_pip->GetYaxis()->SetTitle("#Delta#beta");
  hbeta_pip->GetXaxis()->SetLabelSize(0.05);
  hbeta_pip->GetYaxis()->SetLabelSize(0.05);

  auto* hbeta_pim=new TH2F("hbeta_pim","#Delta#beta for #pi^{-} (post PID)",300,0,10,200,-0.02,0.02);
  hbeta_pim->GetXaxis()->SetTitle("P [GeV/c]");
  hbeta_pim->GetYaxis()->SetTitle("#Delta#beta");
  hbeta_pim->GetXaxis()->SetLabelSize(0.05);
  hbeta_pim->GetYaxis()->SetLabelSize(0.05);

  auto* hbeta_pr=new TH2F("hbeta_pr","#Delta#beta for p (post PID)",300,0,10,200,-0.02,0.02);
  hbeta_pr->GetXaxis()->SetTitle("P [GeV/c]");
  hbeta_pr->GetYaxis()->SetTitle("#Delta#beta");
  hbeta_pr->GetXaxis()->SetLabelSize(0.05);
  hbeta_pr->GetYaxis()->SetLabelSize(0.05);




  Long64_t nentries = t1->GetEntries();
  for(Long64_t i=0; i<nentries;i++){
    t1->GetEvent(i);
    if(readprotonno==1 && readpipno==1 && readpimno==1 && readelno==1){

      TLorentzVector miss = beam + target - v_p4_el->at(0) - v_p4_pip->at(0) - v_p4_pim->at(0) - v_p4_pr->at(0);
      TLorentzVector K0 = v_p4_pim->at(0) + v_p4_pr->at(0);

      Double_t beta_tof_pip = v_beta_pip->at(0);
      Double_t P_pip = sqrt((pow(v_p4_pip->at(0).Px(),2))+(pow(v_p4_pip->at(0).Py(),2))+(pow(v_p4_pip->at(0).Px(),2)));
      Double_t beta_calc_pip = P_pip/(sqrt((pow(P_pip,2))+(pow(v_p4_pip->at(0).M(),2))));
      Double_t delta_beta_pip = beta_calc_pip-beta_tof_pip;

      Double_t beta_tof_pim = v_beta_pim->at(0);
      Double_t P_pim = sqrt((pow(v_p4_pim->at(0).Px(),2))+(pow(v_p4_pim->at(0).Py(),2))+(pow(v_p4_pim->at(0).Px(),2)));
      Double_t beta_calc_pim = P_pim/(sqrt((pow(P_pim,2))+(pow(v_p4_pim->at(0).M(),2))));
      Double_t delta_beta_pim = beta_calc_pim-beta_tof_pim;

      Double_t beta_tof_pr = v_beta_pr->at(0);
      Double_t P_pr = sqrt((pow(v_p4_pr->at(0).Px(),2))+(pow(v_p4_pr->at(0).Py(),2))+(pow(v_p4_pr->at(0).Px(),2)));
      Double_t beta_calc_pr = P_pr/(sqrt((pow(P_pr,2))+(pow(v_p4_pr->at(0).M(),2))));
      Double_t delta_beta_pr = beta_calc_pr-beta_tof_pr;

      if(miss.M()>0.8 && miss.M()<1.1 && abs(beta_tof_pip) < 0.02 && abs(delta_beta_pim) <0.02 &&
         abs(delta_beta_pim) < 0.02 && P_pip>0.2 && P_pim>0.2 && P_pr>0.2){
           hmK0->Fill(K0.M());
           hbeta_pip->Fill(P_pip,delta_beta_pip);
           hbeta_pim->Fill(P_pim,delta_beta_pim);
           hbeta_pr->Fill(P_pr,delta_beta_pr);
         }

    }

  }

TCanvas *can1=new TCanvas("can1","My Plot", 600, 600);
can1->Divide(1,1);
can1->cd(1);
hmK0->Draw();

TCanvas *can2=new TCanvas("can2","My Plot", 600, 600);
can2->Divide(1,1);
can2->cd(1);
hbeta_pip->Draw();

TCanvas *can3=new TCanvas("can3","My Plot", 600, 600);
can3->Divide(1,1);
can3->cd(1);
hbeta_pim->Draw();

TCanvas *can4=new TCanvas("can4","My Plot", 600, 600);
can4->Divide(1,1);
can4->cd(1);
hbeta_pr->Draw();

}
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
