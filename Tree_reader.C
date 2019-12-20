#include <cstdlib>
#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <TH1.h>
#include <TF1.h>
#include <TH2.h>
#include <TChain.h>
#include <TCanvas.h>
#include <vector>
#include <TLorentzVector.h>


void Tree_reader(){



  gROOT->ProcessLine(".L ./Loader.C+");
  TFile *f = new TFile("skim4_Tree1.root");
  TTree *t1 = (TTree*)f->Get("skim4_Tree");

  vector<TLorentzVector> *v_p4=0;


  TLorentzVector *readbeam=NULL;
  TLorentzVector *readtarget=NULL;


  vector<TLorentzVector> *v_vertex=0;

  vector<double> *v_beta=0;

  Double_t start_time;
  vector<double> *energy=0;
  vector<double> *charge=0;
  vector<double> *PID=0;
  vector<double> *chi2PID=0;

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

  TFile fileOutput1("Cleaner_pions_and_proton.root","recreate");

  auto* hmiss=new TH1F("hmiss","Missing Mass;MM(e #pi^{+} #pi^{-}) [GeV/c^{2}];Counts",250,0,4);
  hmiss->GetXaxis()->SetLabelSize(0.05);
  hmiss->GetYaxis()->SetLabelSize(0.05);

  auto* hmiss_sig=new TH1F("hmiss_sig","Missing Mass;MM(e #pi^{+} #pi^{-}) [GeV/c^{2}];Counts",250,0,4);
  hmiss_sig->GetXaxis()->SetLabelSize(0.05);
  hmiss_sig->GetYaxis()->SetLabelSize(0.05);

  auto* hmiss_cut=new TH1F("hmiss_cut","Missing Mass;MM(e #pi^{+} #pi^{-}) [GeV/c^{2}];Counts",250,0,4);
  hmiss_cut->GetXaxis()->SetLabelSize(0.05);
  hmiss_cut->GetYaxis()->SetLabelSize(0.05);

  auto* hmiss_all=new TH1F("hmiss_all","Missing Mass;MM(e #pi^{+} #pi^{-} p) [GeV/c^{2}];Counts",400,-3,3);
  hmiss_all->GetXaxis()->SetLabelSize(0.05);
  hmiss_all->GetYaxis()->SetLabelSize(0.05);

  auto* hThetaVsPhi_pip=new TH2F("hThetaVsPhi_pip","Theta Versus Phi;Phi [rad];Theta [rad]",200,-3.4,3.4,200,0,3.4);
  hThetaVsPhi_pip->GetXaxis()->SetLabelSize(0.05);
  hThetaVsPhi_pip->GetYaxis()->SetLabelSize(0.05);

  auto* hThetaVsPhi=new TH2F("hThetaVsPhi","Theta Versus Phi;Phi [rad];Theta [rad]",200,-3.4,3.4,200,0,3.4);
  hThetaVsPhi->GetXaxis()->SetLabelSize(0.05);
  hThetaVsPhi->GetYaxis()->SetLabelSize(0.05);

  auto* hThetaVsPhi_pim=new TH2F("hThetaVsPhi_pim","Theta Versus Phi;Phi [rad];Theta [rad]",200,-3.4,3.4,200,0,3.4);
  hThetaVsPhi_pim->GetXaxis()->SetLabelSize(0.05);
  hThetaVsPhi_pim->GetYaxis()->SetLabelSize(0.05);

  auto* hThetaVsPhi_pr=new TH2F("hThetaVsPhi_pr","Theta Versus Phi;Phi [rad];Theta [rad]",200,-3.4,3.4,200,0,3.4);
  hThetaVsPhi_pr->GetXaxis()->SetLabelSize(0.05);
  hThetaVsPhi_pr->GetYaxis()->SetLabelSize(0.05);

  auto* hThetaVsPhi_el=new TH2F("hThetaVsPhi_el","Theta Versus Phi",200,-3.4,3.4,200,0,3.4);
  hThetaVsPhi_el->GetXaxis()->SetTitle("Phi [rad]");
  hThetaVsPhi_el->GetXaxis()->SetTitle("Theta [rad]");
  hThetaVsPhi_el->GetXaxis()->SetLabelSize(0.05);
  hThetaVsPhi_el->GetYaxis()->SetLabelSize(0.05);

  auto* hmK0=new TH1F("hmK0","Invariant Mass;M(#pi^{+} #pi^{-}) [GeV/c^{2}];Counts",200,0,2);
  hmK0->GetXaxis()->SetLabelSize(0.05);
  hmK0->GetYaxis()->SetLabelSize(0.05);

  auto* hbeta_pip=new TH2F("hbeta_pip","#Delta#beta for #pi^{+} (post PID);P [GeV/c];#Delta#beta",1000,0,10,200,-0.02,0.02);
  hbeta_pip->GetXaxis()->SetLabelSize(0.05);
  hbeta_pip->GetYaxis()->SetLabelSize(0.05);

  auto* hbeta_pim=new TH2F("hbeta_pim","#Delta#beta for #pi^{-} (post PID);P [GeV/c];#Delta#beta",1000,0,10,200,-0.02,0.02);
  hbeta_pim->GetXaxis()->SetLabelSize(0.05);
  hbeta_pim->GetYaxis()->SetLabelSize(0.05);

  auto* hbeta_pr=new TH2F("hbeta_pr","#Delta#beta for p (post PID);P [GeV/c];#Delta#beta",1000,0,10,200,-0.02,0.02);
  hbeta_pr->GetXaxis()->SetLabelSize(0.05);
  hbeta_pr->GetYaxis()->SetLabelSize(0.05);


// cout<<readtarget->M()<<endl;
  TLorentzVector miss;
  TLorentzVector K0;
  TLorentzVector miss_all;

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


  Long64_t nentries = t1->GetEntries();
  cout<<nentries<<endl;

  for(Long64_t i=0; i<nentries;i++){
    t1->GetEntry(i);
    if (i % 1000 == 0){
      fprintf (stderr, "%lld\r", i/1000);
      fflush (stderr);
        }
      Int_t Np = v_p4->size();
      for(Int_t j=0; j<Np; j++){
        hThetaVsPhi->Fill(v_p4->at(j).Phi(),v_p4->at(j).Theta());
        if(PID->at(j)==11){
          el.SetXYZM(v_p4->at(j).Px(),v_p4->at(j).Py(),v_p4->at(j).Pz(),v_p4->at(j).M());
        }
        else if(PID->at(j)==211){
          pip.SetXYZM(v_p4->at(j).Px(),v_p4->at(j).Py(),v_p4->at(j).Pz(),v_p4->at(j).M());
          beta_tof_pip = v_beta->at(j);
        }
        else if(PID->at(j)==-211){
          pim.SetXYZM(v_p4->at(j).Px(),v_p4->at(j).Py(),v_p4->at(j).Pz(),v_p4->at(j).M());
          beta_tof_pim = v_beta->at(j);
        }
        else if(PID->at(j)==2212){
          pr.SetXYZM(v_p4->at(j).Px(),v_p4->at(j).Py(),v_p4->at(j).Pz(),v_p4->at(j).M());
          beta_tof_pr = v_beta->at(j);
        }

      }
      miss = (TLorentzVector)*readbeam + (TLorentzVector)*readtarget - el - pim - pip;

      miss_all = (TLorentzVector)*readbeam + (TLorentzVector)*readtarget - el - pim - pip - pr;

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

      hmiss_all->Fill(miss_all.M2());

      hThetaVsPhi_pip->Fill(pip.Phi(),pip.Theta());
      hThetaVsPhi_pim->Fill(pim.Phi(),pim.Theta());
      hThetaVsPhi_pr->Fill(pr.Phi(),pr.Theta());

      if(miss.M()>0.8 && miss.M()<1.1 && abs(miss_all.M2())<0.05 && abs(delta_beta_pip) < 0.02 && abs(delta_beta_pr) <0.02 &&
         abs(delta_beta_pim) < 0.02 && P_pip>0.2 && P_pim>0.2 && P_pr>0.2){
            //      cout<<miss.M()<<"   "<<beamA.M()<<"   "<<targetA.M()<<"   "<<el.M()<<endl;



           hbeta_pip->Fill(P_pip,delta_beta_pip);
           hbeta_pim->Fill(P_pim,delta_beta_pim);
           hbeta_pr->Fill(P_pr,delta_beta_pr);
         }
         if(abs(miss_all.M2())<0.05 && abs(delta_beta_pip) < 0.02 && abs(delta_beta_pr) <0.02 &&
            abs(delta_beta_pim) < 0.02 && P_pip>0.2 && P_pim>0.2 && P_pr>0.2){
              hmiss_cut->Fill(miss.M());
            }
        if(K0.M()>0.358 && K0.M()<638 && abs(miss_all.M2())<0.05 && abs(delta_beta_pip) < 0.02 && abs(delta_beta_pr) <0.02 &&
           abs(delta_beta_pim) < 0.02 && P_pip>0.2 && P_pim>0.2 && P_pr>0.2){
             hmiss_sig->Fill(miss.M());
           }
        if(miss.M()>1.1 && miss.M()<1.3 && abs(miss_all.M2())<0.05 && abs(delta_beta_pip) < 0.02 && abs(delta_beta_pr) <0.02 &&
           abs(delta_beta_pim) < 0.02 && P_pip>0.2 && P_pim>0.2 && P_pr>0.2){
            hmK0->Fill(K0.M());
          }

  }

TH1D * prjection_hbeta_pip = hbeta_pip->ProjectionX();
TH1D * prjection_hbeta_pim = hbeta_pim->ProjectionX();
TH1D * prjection_hbeta_pr = hbeta_pr->ProjectionX();

TCanvas *c1=new TCanvas("c1","My plots", 600, 600);
c1->Divide(1,2);
c1->cd(1);
hbeta_pip->Draw("colz");
c1->cd(2);
prjection_hbeta_pip->Draw();

TCanvas *c2=new TCanvas("c2","My plots", 600, 600);
c2->Divide(1,2);
c2->cd(1);
hbeta_pim->Draw("colz");
c2->cd(2);
prjection_hbeta_pim->Draw();

TCanvas *c3=new TCanvas("c3","My plots", 600, 600);
c3->Divide(1,2);
c3->cd(1);
hbeta_pr->Draw("colz");
c3->cd(2);
prjection_hbeta_pr->Draw();

c1->SaveAs("/shared/storage/physhad/JLab/gc1108/work/projPip.jpg");
c2->SaveAs("/shared/storage/physhad/JLab/gc1108/work/projPim.jpg");
c3->SaveAs("/shared/storage/physhad/JLab/gc1108/work/projPr.jpg");


// Int_t a=0;
// Int_t b=0;
// Int_t c=0;
//
// Int_t a_pip=0;
// Int_t S_a_pip=0;
//
// Int_t a_pim=0;
// Int_t S_a_pim=0;
//
// Int_t a_pr=0;
// Int_t S_a_pr=0;
//
//   for(Int_t i;i<1001;i++){
//       a_pip=prjection_hbeta_pip->GetBinContnet(i);
//       S_a_pip=S_a_pip+a_pip;
//
//       if(S_a_pip>100){
//         x_pip[a]=
//       }
//
//       a_pim=rjection_hbeta_pim->GetBinContnet(i);
//       S_a_pim=S_a_pim+a_pim;
//
//       a_pr=rjection_hbeta_pr->GetBinContnet(i);
//       S_a_pr=S_a_pr+a_pr;
//
//     }
//   }


// hbeta_pip->FitSlicesY();
// // the mean and sigma
// TH1D *hbeta_pip_1 = (TH1D*)gDirectory->Get("hbeta_pip_1");
// TH1D *hbeta_pip_2 = (TH1D*)gDirectory->Get("hbeta_pip_2");
//
// hbeta_pim->FitSlicesY();
// // the mean and sigma
// TH1D *hbeta_pim_1 = (TH1D*)gDirectory->Get("hbeta_pim_1");
// TH1D *hbeta_pim_2 = (TH1D*)gDirectory->Get("hbeta_pim_2");
//
// hbeta_pr->FitSlicesY();
// // the mean and sigma
// TH1D *hbeta_pr_1 = (TH1D*)gDirectory->Get("hbeta_pr_1");
// TH1D *hbeta_pr_2 = (TH1D*)gDirectory->Get("hbeta_pr_2");


// TF1 *meanPsig_pip = new TF1("meanPsig_pip","[0]+[1]*x+exp([2]+[3]*x)+exp([4]+[5]*x)",0,10);
// meanPsig_pip->SetParameter(0,-0.000848);
// meanPsig_pip->SetParameter(1,0.0002449);
// meanPsig_pip->SetParameter(2,-1.731);
// meanPsig_pip->SetParameter(3,-5.332);
// meanPsig_pip->SetParameter(4,-5.382);
// meanPsig_pip->SetParameter(5,-0.07091);
//
// TF1 *meanP2sig_pip = new TF1("meanP2sig_pip","[0]+[1]*x+2*(exp([2]+[3]*x)+exp([4]+[5]*x))",0,10);
// meanP2sig_pip->SetParameter(0,-0.000848);
// meanP2sig_pip->SetParameter(1,0.0002449);
// meanP2sig_pip->SetParameter(2,-1.731);
// meanP2sig_pip->SetParameter(3,-5.332);
// meanP2sig_pip->SetParameter(4,-5.382);
// meanP2sig_pip->SetParameter(5,-0.07091);
//
// TF1 *meanP3sig_pip = new TF1("meanP3sig_pip","[0]+[1]*x+3*(exp([2]+[3]*x)+exp([4]+[5]*x))",0,10);
// meanP3sig_pip->SetParameter(0,-0.000848);
// meanP3sig_pip->SetParameter(1,0.0002449);
// meanP3sig_pip->SetParameter(2,-1.731);
// meanP3sig_pip->SetParameter(3,-5.332);
// meanP3sig_pip->SetParameter(4,-5.382);
// meanP3sig_pip->SetParameter(5,-0.07091);
//
// TF1 *meanMsig_pip = new TF1("meanMsig_pip","-meanPsig_pip",0,10);
// TF1 *meanM2sig_pip = new TF1("meanM2sig_pip","-meanP2sig_pip",0,10);
// TF1 *meanM3sig_pip = new TF1("meanM3sig_pip","-meanP3sig_pip",0,10);
//
// TF1 *meanPsig_pim = new TF1("meanPsig_pim","[0]+[1]*x+exp([2]+[3]*x)+exp([4]+[5]*x)",0,10);
// meanPsig_pim->SetParameter(0,-0.001139);
// meanPsig_pim->SetParameter(1,0.0003158);
// meanPsig_pim->SetParameter(2,0.2751);
// meanPsig_pim->SetParameter(3,-7.967);
// meanPsig_pim->SetParameter(4,-5.701);
// meanPsig_pim->SetParameter(5,0.0687);
//
// TF1 *meanP2sig_pim = new TF1("meanP2sig_pim","[0]+[1]*x+2*(exp([2]+[3]*x)+exp([4]+[5]*x))",0,10);
// meanP2sig_pim->SetParameter(0,-0.001139);
// meanP2sig_pim->SetParameter(1,0.0003158);
// meanP2sig_pim->SetParameter(2,0.2751);
// meanP2sig_pim->SetParameter(3,-7.967);
// meanP2sig_pim->SetParameter(4,-5.701);
// meanP2sig_pim->SetParameter(5,0.0687);
//
// TF1 *meanP3sig_pim = new TF1("meanP3sig_pim","[0]+[1]*x+3*(exp([2]+[3]*x)+exp([4]+[5]*x))",0,10);
// meanP3sig_pim->SetParameter(0,-0.001139);
// meanP3sig_pim->SetParameter(1,0.0003158);
// meanP3sig_pim->SetParameter(2,0.2751);
// meanP3sig_pim->SetParameter(3,-7.967);
// meanP3sig_pim->SetParameter(4,-5.701);
// meanP3sig_pim->SetParameter(5,0.0687);
//
// TF1 *meanMsig_pim = new TF1("meanMsig_pim","-meanPsig_pim",0,10);
// TF1 *meanM2sig_pim = new TF1("meanM2sig_pim","-meanP2sig_pim",0,10);
// TF1 *meanM3sig_pim = new TF1("meanM3sig_pim","-meanP3sig_pim",0,10);
//
// TF1 *meanPsig_pr = new TF1("meanPsig_pr","[0]+[1]*x+exp([2]+[3]*x)+exp([4]+[5]*x)",0,10);
// meanPsig_pr->SetParameter(0,0.0004742);
// meanPsig_pr->SetParameter(1,0.0003901);
// meanPsig_pr->SetParameter(2,-3.267);
// meanPsig_pr->SetParameter(3,-3.201);
// meanPsig_pr->SetParameter(4,-6.237);
// meanPsig_pr->SetParameter(5,0.6004);
//
// TF1 *meanP2sig_pr = new TF1("meanP2sig_pr","[0]+[1]*x+2*(exp([2]+[3]*x)+exp([4]+[5]*x))",0,10);
// meanP2sig_pr->SetParameter(0,0.0004742);
// meanP2sig_pr->SetParameter(1,0.0003901);
// meanP2sig_pr->SetParameter(2,-3.267);
// meanP2sig_pr->SetParameter(3,-3.201);
// meanP2sig_pr->SetParameter(4,-6.237);
// meanP2sig_pr->SetParameter(5,0.6004);
//
// TF1 *meanP3sig_pr = new TF1("meanP3sig_pr","[0]+[1]*x+3*(exp([2]+[3]*x)+exp([4]+[5]*x))",0,10);
// meanP3sig_pr->SetParameter(0,0.0004742);
// meanP3sig_pr->SetParameter(1,0.0003901);
// meanP3sig_pr->SetParameter(2,-3.267);
// meanP3sig_pr->SetParameter(3,-3.201);
// meanP3sig_pr->SetParameter(4,-6.237);
// meanP3sig_pr->SetParameter(5,0.6004);
//
// TF1 *meanMsig_pr = new TF1("meanMsig_pr","-meanPsig_pr",0,10);
// TF1 *meanM2sig_pr = new TF1("meanM2sig_pr","-meanP2sig_pr",0,10);
// TF1 *meanM3sig_pr = new TF1("meanM3sig_pr","-meanP3sig_pr",0,10);
//
// TCanvas *c1=new TCanvas("c1","My plots", 600, 600);
// c1->Divide(1,1);
// c1->cd(1);
// hbeta_pip->DrawCopy("colz");
//
// meanPsig_pip->SetLineColor(1);
// meanPsig_pip->Draw("same");
//
// meanMsig_pip->SetLineColor(1);
// meanMsig_pip->Draw("same");
//
// meanP2sig_pip->SetLineColor(2);
// meanP2sig_pip->Draw("same");
//
// meanM2sig_pip->SetLineColor(2);
// meanM2sig_pip->Draw("same");
//
// meanP3sig_pip->SetLineColor(3);
// meanP3sig_pip->Draw("same");
//
// meanM3sig_pip->SetLineColor(3);
// meanM3sig_pip->Draw("same");
//
// // c1->SaveAs("/shared/storage/physhad/JLab/gc1108/work/root_stuff/clas12root-master/RunRoot/c1_1.pdf");
//
//
// TCanvas *c2=new TCanvas("c2","My plots", 600, 600);
// c2->Divide(1,1);
// c2->cd(1);
// hbeta_pim->Draw("colz");
//
// meanPsig_pim->SetLineColor(1);
// meanPsig_pim->Draw("same");
//
// meanMsig_pim->SetLineColor(1);
// meanMsig_pim->Draw("same");
//
// meanP2sig_pim->SetLineColor(2);
// meanP2sig_pim->Draw("same");
//
// meanM2sig_pim->SetLineColor(2);
// meanM2sig_pim->Draw("same");
//
// meanP3sig_pim->SetLineColor(3);
// meanP3sig_pim->Draw("same");
//
// meanM3sig_pim->SetLineColor(3);
// meanM3sig_pim->Draw("same");
//
// TCanvas *c3=new TCanvas("c3","My plots", 600, 600);
// c3->Divide(1,1);
// c3->cd(1);
// hbeta_pr->Draw("colz");
//
// meanPsig_pr->SetLineColor(1);
// meanPsig_pr->Draw("same");
//
// meanMsig_pr->SetLineColor(1);
// meanMsig_pr->Draw("same");
//
// meanP2sig_pr->SetLineColor(2);
// meanP2sig_pr->Draw("same");
//
// meanM2sig_pr->SetLineColor(2);
// meanM2sig_pr->Draw("same");
//
// meanP3sig_pr->SetLineColor(3);
// meanP3sig_pr->Draw("same");
//
// meanM3sig_pr->SetLineColor(3);
// meanM3sig_pr->Draw("same");



fileOutput1.Write();


}
