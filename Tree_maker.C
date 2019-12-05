#include <cstdlib>
#include <iostream>
#include <chrono>
#include <TFile.h>
#include <TTree.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <TH1.h>
#include <TH2.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TBenchmark.h>
#include "clas12reader.h"
#include <vector>
using namespace clas12;


        void Tree_maker(){
          gROOT->ProcessLine(".L ./Loader.C+");

          //Telling which files to run over
          TString inputFile("/work/clas12/rg-a/trains/v16_v2/skim4_inclusive/skim4_*.hipo");

          //Creating a chain for the data from different files
          TChain fake("hipo");
          fake.Add(inputFile.Data());

          //get the hipo data
          auto files=fake.GetListOfFiles();

          //Name of Tree file
          TFile f("skim4_Tree1.root","recreate");

          //Creating TTree object
          TTree skim4_Tree("skim4_Tree","a tree mate");

          //Creating variables and branches

          //Particle TLorentzVectors
          vector<TLorentzVector> v_p4_el;
          TLorentzVector p4_el;
          vector<TLorentzVector> v_p4_pip;
          TLorentzVector p4_pip;
          vector<TLorentzVector> v_p4_pim;
          TLorentzVector p4_pim;
          vector<TLorentzVector> v_p4_pr;
          TLorentzVector p4_pr;

          //Vertex position and time
          vector<TLorentzVector> v_vertex;
          TLorentzVector vertex;

          //Particle beta_tof
          Double_t beta_el;
          vector<double> v_beta_el;
          Double_t beta_pip;
          vector<double> v_beta_pip;
          Double_t beta_pim;
          vector<double> v_beta_pim;
          Double_t beta_pr;
          vector<double> v_beta_pr;

          Double_t start_time;
          vector<double> energy;
          vector<double> P;
          vector<double> charge;
          vector<double> PID;
          vector<double> chi2PID;
          Int_t chargetracks;
          Int_t protonno;
          Int_t kaonpno;
          Int_t pipno;
          Int_t pimno;
          Int_t elno;
          auto db=TDatabasePDG::Instance();
          TLorentzVector beam(0,0,10.6,10.6);
          TLorentzVector target(0,0,0,db->GetParticle(2212)->Mass());


          skim4_Tree.Branch("p4_el",&v_p4_el);
          skim4_Tree.Branch("p4_pip",&v_p4_pip);
          skim4_Tree.Branch("p4_pim",&v_p4_pim);
          skim4_Tree.Branch("p4_pr",&v_p4_pr);
          skim4_Tree.Branch("start_time",&start_time);
          skim4_Tree.Branch("beta_el",&v_beta_el);
          skim4_Tree.Branch("beta_pip",&v_beta_pip);
          skim4_Tree.Branch("beta_pim",&v_beta_pim);
          skim4_Tree.Branch("beta_pr",&v_beta_pr);
          skim4_Tree.Branch("energy",&energy);
          skim4_Tree.Branch("P",&P);
          skim4_Tree.Branch("charge",&charge);
          skim4_Tree.Branch("vertex",&v_vertex);
          skim4_Tree.Branch("PID",&PID);
          skim4_Tree.Branch("chi2PID",&chi2PID);
          skim4_Tree.Branch("chargetracks",&chargetracks,"chargetracks/I");
          skim4_Tree.Branch("protonno",&protonno,"protonno/I");
          skim4_Tree.Branch("kaonpno",&kaonpno,"kaonpno/I");
          skim4_Tree.Branch("pipno",&pipno,"pipno/I");
          skim4_Tree.Branch("pimno",&pimno,"pimno/I");
          skim4_Tree.Branch("elno",&elno,"elno/I");
          skim4_Tree.Branch("beam",&beam);
          skim4_Tree.Branch("target",&target);

          //Going over all files listed above and all events inside those
          for(Int_t i=0;i<files->GetEntries();i++){
            //create the event reader
            clas12reader c12(files->At(i)->GetTitle());

            //While loop covering multiple decisions on what data to store
            while(c12.next()==true){
              v_p4_el.clear();
              v_p4_pip.clear();
              v_p4_pim.clear();
              v_p4_pr.clear();

              v_beta_el.clear();
              v_beta_pip.clear();
              v_beta_pim.clear();
              v_beta_pr.clear();

              energy.clear();
              P.clear();
              charge.clear();
              v_vertex.clear();
              PID.clear();
              chi2PID.clear();

              //Getting start time for each event
              start_time=c12.event()->getStartTime();

              chargetracks = 0;
              protonno = 0;
              kaonpno = 0;
              pipno = 0;
              pimno = 0;
              elno = 0;

              for(auto& p : c12.getDetParticles()){
                 //  get predefined selected information

                 p4_el.SetXYZM(p->par()->getPx(), p->par()->getPy(), p->par()->getPz(),0);
                 p4_pip.SetXYZM(p->par()->getPx(), p->par()->getPy(), p->par()->getPz(),0);
                 p4_pim.SetXYZM(p->par()->getPx(), p->par()->getPy(), p->par()->getPz(),0);
                 p4_pr.SetXYZM(p->par()->getPx(), p->par()->getPy(), p->par()->getPz(),0);

                 vertex.SetXYZT(p->par()->getVx(), p->par()->getVy(), p->par()->getVz(),0);


                 //Setting 4 vector of momentum if it passes the deltabeta cuts
                 if(p->par()->getPid()==2212){
                   p4_pr.SetXYZM(p4_pr.Px(),p4_pr.Py(),p4_pr.Pz(),db->GetParticle(2212)->Mass());

                   p4_el.SetXYZM(0,0,0,0);
                   p4_pip.SetXYZM(0,0,0,0);
                   p4_pim.SetXYZM(0,0,0,0);

                   beta_pr=p->par()->getBeta();

                   beta_el=0;
                   beta_pip=0;
                   beta_pim=0;
                   //count the number of protons in this event
                   protonno=protonno+1;

                 }


                else if(p->par()->getPid()==211){
                  p4_pip.SetXYZM(p4_pip.Px(),p4_pip.Py(),p4_pip.Pz(),db->GetParticle(211)->Mass());

                  p4_el.SetXYZM(0,0,0,0);
                  p4_pr.SetXYZM(0,0,0,0);
                  p4_pim.SetXYZM(0,0,0,0);

                  beta_pip=p->par()->getBeta();
                  //count the number of positve pions in this event
                  pipno=pipno+1;

                }


                else if(p->par()->getPid()==-211){
                  p4_pim.SetXYZM(p4_pim.Px(),p4_pim.Py(),p4_pim.Pz(),db->GetParticle(-211)->Mass());

                  p4_el.SetXYZM(0,0,0,0);
                  p4_pip.SetXYZM(0,0,0,0);
                  p4_pr.SetXYZM(0,0,0,0);

                  beta_pim=p->par()->getBeta();

                  beta_el=0;
                  beta_pip=0;
                  beta_pr=0;
                  //count the number of negative pions in this event
                  pimno=pimno+1;

                }
                else if(p->par()->getPid()==11){
                  p4_el.SetXYZM(p4_el.Px(),p4_el.Py(),p4_el.Pz(),db->GetParticle(11)->Mass());

                  p4_pr.SetXYZM(0,0,0,0);
                  p4_pip.SetXYZM(0,0,0,0);
                  p4_pim.SetXYZM(0,0,0,0);

                  beta_el=p->par()->getBeta();

                  beta_pr=0;
                  beta_pip=0;
                  beta_pim=0;
                  //count the number of negative pions in this event
                  elno=elno+1;
                }



               //counting number of charge particles in each event
               if(p->par()->getCharge()!=0){
                 chargetracks = chargetracks + 1;
               }

               if(p4_el.M()>0.0001 || p4_pip.M()>0.0001 || p4_pim.M()>0.0001 || p4_pr.M()>0.0001){
                 v_p4_el.push_back(p4_el);
                 v_p4_pip.push_back(p4_pip);
                 v_p4_pim.push_back(p4_pim);
                 v_p4_pr.push_back(p4_pr);
                 v_beta_el.push_back(beta_el);
                 v_beta_pip.push_back(beta_pip);
                 v_beta_pim.push_back(beta_pim);
                 v_beta_pr.push_back(beta_pr);
                 energy.push_back(p->getDetEnergy());
                 P.push_back(p->par()->getP());
                 charge.push_back(p->par()->getCharge());
                 v_vertex.push_back(vertex);
                 PID.push_back(p->par()->getPid());
                 chi2PID.push_back(p->par()->getChi2Pid());
               }

}
                   //if statement to decide which data to store based on number of charged particles in event
                   if(chargetracks>2){
                     //Fill skim4_Tree
                     skim4_Tree.Fill();

                   }
                 }
               }

               skim4_Tree.Write();
             }

