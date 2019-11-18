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
            TString inputFile("/work/clas12/rg-a/trains/v16_v2/skim4_inclusive/skim4_5201.hipo");

          //Creating a chain for the data from different files
          TChain fake("hipo");
          fake.Add(inputFile.Data());

          //get the hipo data
          auto files=fake.GetListOfFiles();

          //Name of Tree file
          TFile f("skim4_5201_Tree1.root","recreate");

          //Creating TTree object
          TTree skim4_5201_Tree("skim4_5201_Tree","a tree mate");




          //Creating variables and branches
          vector<TLorentzVector> v_p4;
          TLorentzVector p4;
          //vector<TLorentzVector> v_vertex;
          //TLorentzVector vertex;
          vector<double> beta;
          Double_t start_time;
          vector<double> energy, P;
          vector<double> charge;
          vector<double> vertex_position;
          vector<double> chi2PID;

          skim4_5201_Tree.Branch("p4",&v_p4);
          skim4_5201_Tree.Branch("start_time",&start_time);
          skim4_5201_Tree.Branch("beta",&beta);
          skim4_5201_Tree.Branch("energy",&energy);
          skim4_5201_Tree.Branch("P",&P);
          skim4_5201_Tree.Branch("charge",&charge);
          //skim4_5201_Tree.Branch("vertex",&v_vertex);
          skim4_5201_Tree.Branch("chi2PID",&chi2PID);



          //Going over all files listed above and all events inside those
          for(Int_t i=0;i<files->GetEntries();i++){
            //create the event reader
            clas12reader c12(files->At(i)->GetTitle());



            //While loop covering multiple decisions on what data to store
            while(c12.next()==true){
              v_p4.clear();
              beta.clear();
              energy.clear();
              P.clear();
              charge.clear();
              //v_vertex.clear();
              chi2PID.clear();

              //Getting start time for each event
              start_time=c12.event()->getStartTime();

               Int_t chargetracks = 0;
               Int_t protonno = 0;
               Int_t kaonpno = 0;
               Int_t pipno = 0;
               Int_t pimno = 0;

              for(auto& p : c12.getDetParticles()){
                 //  get predefined selected information

                 p4.SetXYZM(p->par()->getPx(), p->par()->getPy(), p->par()->getPz(), 0);
                 //vertex.SetXYZT(p->par()->getVx(), p->par()->getVy(), p->par()->getVz(), p->par()->getVt());


                 //Deltabetas for different particles
                 Double_t Beta_calcpr=p->par()->getP()/sqrt(pow(0.938,2)+pow(p->par()->getP(),2));
                 Double_t DeltaBetapr=Beta_calcpr - p->par()->getBeta();

                 Double_t Beta_calckp=p->par()->getP()/sqrt(pow(0.49368,2)+pow(p->par()->getP(),2));
                 Double_t DeltaBetakp=Beta_calckp-p->par()->getBeta();

                 Double_t Beta_calcpic=p->par()->getP()/sqrt(pow(0.13957,2)+pow(p->par()->getP(),2));
                 Double_t DeltaBetapic=Beta_calcpic-p->par()->getBeta();


                 //Setting 4 vector of momentum if it passes the deltabeta cuts
                 if(p->par()->getPid()==2212 && abs(DeltaBetapr)<0.02){
                   p4.SetXYZM(p4.Px(),p4.Py(), p4.Pz(),0.938);
                   //count the number of protons in this event
                   protonno=protonno+1;

                 }


                 else if(p->par()->getPid()==321 && abs(DeltaBetakp)<0.02){
                   p4.SetXYZM(p4.Px(),p4.Py(), p4.Pz(),0.49368);
                   //count the number of positve kaons in this event
                   kaonpno=kaonpno+1;

                 }

                else if(p->par()->getPid()==211 && abs(DeltaBetapic)<0.02){
                  p4.SetXYZM(p4.Px(),p4.Py(), p4.Pz(),0.13957);
                  //count the number of positve pions in this event
                  pipno=pipno+1;

                }


                else if(p->par()->getPid()==-211 && abs(DeltaBetapic)<0.02){
                  p4.SetXYZM(p4.Px(),p4.Py(), p4.Pz(),0.13957);
                  //count the number of negative pions in this event
                  pimno=pimno+1;

                }



               //counting number of charge particles in each event
               if(p->par()->getCharge()!=0){
                 chargetracks = chargetracks + 1;


               }

               if(p4.M()>0.01){
                 v_p4.push_back(p4);
                 beta.push_back(p->par()->getBeta());
                 energy.push_back(p->getDetEnergy());
                 P.push_back(p->par()->getP());
                 charge.push_back(p->par()->getCharge());
                 //v_vertex.push_back(vertex);
                 chi2PID.push_back(p->par()->getChi2Pid());
            }

}
                   //for loop to decide which data to store based on number of charged particles in event
                   if(chargetracks>1){
                     //Fill skim4_5201_Tree
                     skim4_5201_Tree.Fill();





                   }





                 }
               }

               skim4_5201_Tree.Write();
             }