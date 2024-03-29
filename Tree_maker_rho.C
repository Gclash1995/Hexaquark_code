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


        void Tree_maker_rho(){
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
          vector<TLorentzVector> v_p4;
          TLorentzVector p4;

          //Vertex position and time
          vector<TLorentzVector> v_vertex;
          TLorentzVector vertex;

          //Particle beta_tof
          vector<double> v_beta;
          Double_t beta;

          Double_t start_time;
          vector<double> energy;
          vector<double> charge;
          vector<double> PID;
          vector<double> chi2PID;
          Int_t chargetracks;
          Int_t protonno;
          Int_t kaonpno;
          Int_t kaonmno;
          Int_t pipno;
          Int_t pimno;
          Int_t elno;
          auto db=TDatabasePDG::Instance();
          TLorentzVector beam(0,0,10.6,10.6);
          TLorentzVector target(0,0,0,db->GetParticle(2212)->Mass());


          skim4_Tree.Branch("p4",&v_p4);

          skim4_Tree.Branch("beta",&v_beta);

          skim4_Tree.Branch("start_time",&start_time);
          skim4_Tree.Branch("energy",&energy);
          skim4_Tree.Branch("charge",&charge);
          skim4_Tree.Branch("vertex",&v_vertex);
          skim4_Tree.Branch("PID",&PID);
          skim4_Tree.Branch("chi2PID",&chi2PID);
          skim4_Tree.Branch("chargetracks",&chargetracks,"chargetracks/I");

          skim4_Tree.Branch("protonno",&protonno,"protonno/I");
          skim4_Tree.Branch("kaonpno",&kaonpno,"kaonpno/I");
          skim4_Tree.Branch("kaonmno",&kaonmno,"kaonmno/I");
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
              v_p4.clear();

              v_beta.clear();

              energy.clear();
              charge.clear();
              v_vertex.clear();
              PID.clear();
              chi2PID.clear();

              //Getting start time for each event
              start_time=c12.event()->getStartTime();

              protonno = 0;
              kaonpno = 0;
              kaonmno = 0;
              pipno = 0;
              pimno = 0;
              elno = 0;

              for(auto& p : c12.getDetParticles()){

                 //  get predefined selected information

                   p4.SetXYZM(p->par()->getPx(), p->par()->getPy(), p->par()->getPz(),0);

                   vertex.SetXYZT(p->par()->getVx(), p->par()->getVy(), p->par()->getVz(),0);


                   //Setting 4 vector of momentum if it passes the deltabeta cuts
                   if(p->par()->getPid()==2212){
                     p4.SetXYZM(p4.Px(),p4.Py(),p4.Pz(),db->GetParticle(2212)->Mass());

                     beta=p->par()->getBeta();

                     //count the number of protons in this event
                     protonno=protonno+1;

                   }

                  else if(p->par()->getPid()==211){
                    p4.SetXYZM(p4.Px(),p4.Py(),p4.Pz(),db->GetParticle(211)->Mass());

                    beta=p->par()->getBeta();
                    //count the number of positve pions in this event
                    pipno=pipno+1;

                  }

                  else if(p->par()->getPid()==-211){
                    p4.SetXYZM(p4.Px(),p4.Py(),p4.Pz(),db->GetParticle(-211)->Mass());

                    beta=p->par()->getBeta();

                    //count the number of negative pions in this event
                    pimno=pimno+1;

                  }

                  else if(p->par()->getPid()==11){
                    p4.SetXYZM(p4.Px(),p4.Py(),p4.Pz(),db->GetParticle(11)->Mass());

                    beta=p->par()->getBeta();

                    //count the number of negative pions in this event
                    elno=elno+1;
                  }

                  else if(p->par()->getPid()==321){
                    p4.SetXYZM(p4.Px(),p4.Py(),p4.Pz(),db->GetParticle(321)->Mass());

                    beta=p->par()->getBeta();

                    //count the number of negative pions in this event
                    kaonpno=kaonpno+1;
                  }

                  else if(p->par()->getPid()==321){
                    p4.SetXYZM(p4.Px(),p4.Py(),p4.Pz(),db->GetParticle(321)->Mass());

                    beta=p->par()->getBeta();

                    //count the number of negative pions in this event
                    kaonpno=kaonpno+1;
                  }

                  else if(p->par()->getPid()==-321){
                    p4.SetXYZM(p4.Px(),p4.Py(),p4.Pz(),db->GetParticle(-321)->Mass());

                    beta=p->par()->getBeta();

                    //count the number of negative pions in this event
                    kaonmno=kaonmno+1;
                  }


                   v_p4.push_back(p4);

                   v_beta.push_back(beta);

                   energy.push_back(p->getDetEnergy());
                   charge.push_back(p->par()->getCharge());
                   v_vertex.push_back(vertex);
                   PID.push_back(p->par()->getPid());
                   chi2PID.push_back(p->par()->getChi2Pid());


                 }
                   //if statement to decide which data to store based on number of certain particles.
                   if(elno == 1 && pipno == 1 && pimno == 1 && protonno == 1){
                     //Fill skim4_Tree
                     skim4_Tree.Fill();

                   }
                 }
               }

               skim4_Tree.Write();
             }
