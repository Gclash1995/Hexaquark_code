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

using namespace clas12;


void SetLorentzVector(TLorentzVector &p4,clas12::region_part_ptr rp){
  p4.SetXYZM(rp->par()->getPx(),rp->par()->getPy(),
	      rp->par()->getPz(),p4.M());

}

void Ex1_CLAS12Reader_edit3(){
  // Record start time
  auto start = std::chrono::high_resolution_clock::now();


  /////////////////////////////////////
  //ignore this just getting file name!
   TString inputFile;
   TString outputFile;

   for(Int_t i=1;i<gApplication->Argc();i++){
    TString opt=gApplication->Argv(i);
    if((opt.Contains(".hipo"))){
      inputFile=opt(5,opt.Sizeof());
    }
   }
   if(inputFile==TString())  {
      std::cout << " *** please provide a file name..." << std::endl;
     exit(0);
   }
   /////////////////////////////////////


   cout<<"Analysing hipo file "<<inputFile<<endl;

   TChain fake("hipo");
   fake.Add(inputFile.Data());
   //get the hipo data
   //   reader.open(inputFile.Data());
   auto files=fake.GetListOfFiles();

   //some particles
   auto db=TDatabasePDG::Instance();
   TLorentzVector beam(0,0,10.6,10.6);
   TLorentzVector target(0,0,0,db->GetParticle(2212)->Mass());
   TLorentzVector el(0,0,0,db->GetParticle(11)->Mass());

   TLorentzVector pip(0,0,0,db->GetParticle(211)->Mass());
   TLorentzVector pim(0,0,0,db->GetParticle(-211)->Mass());

   auto* hmiss=new TH1F("missM","missM",200,0,3);
   //auto* hmiss2=new TH1F("missM2","missM2",200,0,2);
   //auto* hmiss3=new TH1F("missM3","missM3",200,-1,1);
   auto* hmK0=new TH1F("mK0","mK0",100,0,1);
   auto* hBeta=new TH2F("mBeta","mBeta",100,0,3,100,-2,2);
   auto* hBeta1=new TH2F("hBeta1","hBeta1",100,0,3,100,-2,2);
   auto* hBeta2=new TH2F("hBeta2","hBeta2",100,0,3,100,-2,2);
   auto* hBeta3=new TH2F("hBeta3","hBeta3",100,0,3,100,-2,2);
   auto* hBeta4=new TH2F("hBeta4","hBeta4",100,0,3,100,-2,2);
   auto* hBeta5=new TH2F("hBeta5","hBeta5",100,0,3,100,-2,2);
   //auto* hm2gCut=new TH1F("m2gCut","m2g",200,0,1);

   gBenchmark->Start("timer");
   int counter=0;


for(Int_t i=0;i<files->GetEntries();i++){
     //create the event reader
 clas12reader c12(files->At(i)->GetTitle());
     //  clas12reader c12(files->At(i)->GetTitle(),{0});//add tags {tag1,tag2,tag3,...}

      //Add some event Pid based selections
      //////////c12.AddAtLeastPid(211,1); //at least 1 pi+
      //c12.addExactPid(11,1);    //exactly 1 electron
      //c12.addExactPid(211,1);    //exactly 1 pi+
      //c12.addExactPid(-211,1);    //exactly 1 pi-
      //c12.addExactPid(2212,1);    //exactly 1 proton
      //c12.addExactPid(22,2);    //exactly 2 gamma
      //////c12.addZeroOfRestPid();  //nothing else
      //////c12.useFTBased(); //and use the Pids from RECFT

 while(c12.next()==true){
       // c12.event()->getStartTime(); //hipo4
       // c12.head()->getStartTime(); //hipo3
        //Loop over all particles to see how to access detector info.

 	for(auto& p : c12.getDetParticles()){
  	 //  get predefined selected information
   p->getTime();
   p->getBeta();
	 p->getP();
   p->par()->getCharge();
	 p->getDetEnergy();
	 p->getDeltaEnergy();

	 Double_t Beta_calc=p->getP()/sqrt(pow(0.13957,2)+pow(p->getP(),2));
   Double_t DeltaBeta=Beta_calc-p->getBeta();

   if(p->par()->getCharge()>0)hBeta1->Fill(p->getP(),DeltaBeta);
   if(p->par()->getCharge()>0)hBeta2->Fill(p->getP(),Beta_calc);
   if(p->par()->getCharge()>0)hBeta3->Fill(p->getP(),p->getBeta());


/*
         p->getPx();
         p->getPy();
         p->getPz();
*/



//get beta (getBeta()) this is beta from TOF
//calc beta from momentum assuming pion mass
//For + or negative charges fill 2D histogram Delta Beta vs P
//getCharge()

	 // get any detector information (if exists for this particle)
	 // there should be a get function for any entry in the bank
	 switch(p->getRegion()) {
	 case FD :
	   p->cal(PCAL)->getEnergy();
	   p->cal(ECIN)->getEnergy();
	   p->cal(ECOUT)->getEnergy();
	   p->sci(FTOF1A)->getEnergy();
	   p->sci(FTOF1B)->getEnergy();
	   p->sci(FTOF2)->getEnergy();
	   p->trk(DC)->getSector();
	   p->che(HTCC)->getNphe();
	   p->che(LTCC)->getNphe();
	   //trajectories
	   p->traj(LTCC)->getX();
	   // p->traj(DC,DC1)->getCx();; //First layer of DC, hipo4
	   break;
	 case FT :
	   p->ft(FTCAL)->getEnergy();
	   p->ft(FTHODO)->getEnergy();
	   break;
	 case CD:
	   p->sci(CTOF)->getEnergy();
	   p->sci(CND)->getEnergy();
	   break;
	 }
	 //   covariance matrix (comment in to see!)
	 // p->covmat()->print();
	 p->cmat();
       }

       // get particles by type
       auto electrons=c12.getByID(11);
       auto protons=c12.getByID(2212);
       auto pips=c12.getByID(211);
       auto pims=c12.getByID(-211);

 if(electrons.size() == 1 && protons.size() == 0 && pips.size() == 1 && pims.size() == 1){

	 // set the particle momentum
	 SetLorentzVector(el,electrons[0]);

	 SetLorentzVector(pip,pips[0]);
	 SetLorentzVector(pim,pims[0]);

	 TLorentzVector miss=beam+target-el-pip-pim;
         //TLorentzVector miss2=beam+target-el-pr;
         //TLorentzVector miss3=beam+target-el-pr-pip-pim;
         //hmiss2->Fill(miss2.M2());
         //hmiss3->Fill(miss3.M2());
	 hmiss->Fill(miss.M2());
	 TLorentzVector K0 = pip+pim;
	 hmK0->Fill(K0.M());
	 //if(TMath::Abs(miss.M2())<0.5)hm2gCut->Fill(pi0.M());
   Double_t DeltaBeta=pip.Beta()-pips[0]->getBeta();
   hBeta->Fill(pip.P(),DeltaBeta);
   hBeta4->Fill(pip.P(),pip.Beta());
   hBeta5->Fill(pip.P(),pips[0]->getBeta());

	 //could also get particle time etc. here too
	 //Double_t eTime=electrons[0]->sci(FTOF1A)->getTime();
       }



       counter++;
     }
   }
   gBenchmark->Stop("timer");
   gBenchmark->Print("timer");
   TCanvas* can=new TCanvas();
   can->Divide(2,2);
   can->cd(1);
   hmiss->DrawCopy();
   can->cd(2);
   hmK0->DrawCopy();
   can->cd(3);
   hBeta->DrawCopy("colz");
   can->cd(4);
   hBeta1->DrawCopy("colz");

   TCanvas* can2=new TCanvas();
   can2->Divide(2,1);
   can2->cd(1);
   hBeta2->DrawCopy("colz");
   can2->cd(2);
   hBeta3->DrawCopy("colz");
   can2->cd(3);
   hBeta4->DrawCopy("colz");
   can2->cd(4);
   hBeta5->DrawCopy("colz");
   //hm2gCut->SetLineColor(2);
   //hm2gCut->DrawCopy("same");

   auto finish = std::chrono::high_resolution_clock::now();
   std::chrono::duration<double> elapsed = finish - start;
   std::cout << "Elapsed time: " << elapsed.count()<< " events = "<<counter<< " s\n";

}
