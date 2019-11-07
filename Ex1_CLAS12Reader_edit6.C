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

void Ex1_CLAS12Reader_edit5(){
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
  //TString inputFile("/work/clas12/rg-a/trains/v16_v2/skim4_inclusive/skim4_*.hipo");

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
   TLorentzVector pr(0,0,0,db->GetParticle(2212)->Mass());

   TFile fileOutput1("Testing_cuts.root","recreate");

   auto* hmiss=new TH1F("hmissM","Missing Mass #Sigma^{+}",200,0,3);
   auto* hmK0=new TH1F("hmK0","Invariant Mass",100,0,1);
   auto* hBeta_pip=new TH2F("hBeta_pip","#Delta#beta for #pi^{+} (post PID)",100,0,3,200,-0.1,0.1);
   auto* hBeta_pim=new TH2F("hBeta_pim","#Delta#beta for #pi^{-} (post PID)",100,0,3,200,-0.1,0.1);
   auto* hBeta_pr=new TH2F("hBeta_pr","#Delta#beta for P (post PID)",100,0,3,200,-2,2);


   hmiss->GetXaxis()->SetTitle("MM(e #pi^{+} #pi^{-}) [GeV/c^{2}]");
   hmiss->GetYaxis()->SetTitle("Counts");
   hmiss->GetXaxis()->SetLabelSize(0.05);
   hmiss->GetYaxis()->SetLabelSize(0.05);

   hmK0->GetXaxis()->SetTitle("M(#pi^{+} #pi^{-}) [GeV/c^{2}]");
   hmK0->GetYaxis()->SetTitle("Counts");
   hmK0->GetXaxis()->SetLabelSize(0.05);
   hmK0->GetYaxis()->SetLabelSize(0.05);

   hBeta_pip->GetXaxis()->SetTitle("P [GeV/c]");
   hBeta_pip->GetYaxis()->SetTitle("#Delta#beta");
   hBeta_pip->GetXaxis()->SetLabelSize(0.05);
   hBeta_pip->GetYaxis()->SetLabelSize(0.05);

   hBeta_pim->GetXaxis()->SetTitle("P [GeV/c]");
   hBeta_pim->GetYaxis()->SetTitle("#Delta#beta");
   hBeta_pim->GetXaxis()->SetLabelSize(0.05);
   hBeta_pim->GetYaxis()->SetLabelSize(0.05);

   hBeta_pr->GetXaxis()->SetTitle("P [GeV/c]");
   hBeta_pr->GetYaxis()->SetTitle("#Delta#beta");
   hBeta_pr->GetXaxis()->SetLabelSize(0.05);
   hBeta_pr->GetYaxis()->SetLabelSize(0.05);


/*
   hBeta1->GetXaxis()->SetTitle("P [GeV/c]");
   hBeta1->GetYaxis()->SetTitle("#Delta#beta");

   hBeta2->GetXaxis()->SetTitle("P [GeV/c]");
   hBeta2->GetYaxis()->SetTitle("#beta_{calc}");

   hBeta3->GetXaxis()->SetTitle("P [GeV/c]");
   hBeta3->GetYaxis()->SetTitle("#beta_{tof}");

   hBeta4->GetXaxis()->SetTitle("P [GeV/c]");
   hBeta4->GetYaxis()->SetTitle("#beta_{calc}");

   hBeta5->GetXaxis()->SetTitle("P [GeV/c]");
   hBeta5->GetYaxis()->SetTitle("#beta_{tof}");
*/

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
   p->par()->getBeta();
   p->par()->getFTBBeta();
	 p->par()->getP();
   p->par()->getCharge();
	 p->getDetEnergy();
	 p->getDeltaEnergy();


/*
Delta Beta pre PID

	 Double_t Beta_calc=p->par()->getP()/sqrt(pow(0.13957,2)+pow(p->par()->getP(),2));
   Double_t DeltaBeta=Beta_calc-p->getBeta();

   if(p->par()->getCharge()>0)hBeta1->Fill(p->par()->getP(),DeltaBeta);

   */
   //if(p->par()->getCharge()>0)hBeta2->Fill(p->par()->getP(),Beta_calc);
   //if(p->par()->getCharge()>0)hBeta3->Fill(p->par()->getP(),p->par()->getBeta());

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
       auto pims=c12.getByID(-211);
       auto pips=c12.getByID(211);
       auto protons=c12.getByID(2212);


 if(electrons.size() == 1 && pips.size() == 1 && pims.size() == 1 && protons.size()==1){

	 // set the particle momentum
	 SetLorentzVector(el,electrons[0]);
	 SetLorentzVector(pip,pips[0]);
	 SetLorentzVector(pim,pims[0]);
   SetLorentzVector(pr,protons[0]);

   TLorentzVector miss=beam+target-el-pip-pim;
   TLorentzVector K0 = pip+pim;

	 //if(TMath::Abs(miss.M2())<0.5)hm2gCut->Fill(pi0.M());
   Double_t DeltaBeta_pip=pip.Beta()-pips[0]->par()->getBeta();

   //if(DeltaBeta_pip>-0.06 && DeltaBeta_pip<0.06 && pip.P()>0.18)

   Double_t DeltaBeta_pim=pim.Beta()-pims[0]->par()->getBeta();

   Double_t DeltaBeta_pr=pr.Beta()-pr[0]->par()->getBeta();

   //if(DeltaBeta_pim>-0.06 && DeltaBeta_pim<0.06 && pim.P()>0.18)
   if(miss.M()>0.8 && miss.M()<1.1){
   hBeta_pip->Fill(pip.P(),DeltaBeta_pip);
   hBeta_pim->Fill(pim.P(),DeltaBeta_pim);
   hBeta_pr->Fill(pr.P(),DeltaBeta_pr);
 }



	 // if(DeltaBeta_pim>-0.04 && DeltaBeta_pim<0.04 && pim.P()>0.18 && DeltaBeta_pip>-0.04 && DeltaBeta_pip<0.04 &&
   //    miss.M()<1.1 && miss.M()>0.8){
   //   hmiss->Fill(miss.M());
   //   hmK0->Fill(K0.M());
   // }


	 //could also get particle time etc. here too
	 //Double_t eTime=electrons[0]->sci(FTOF1A)->getTime();
       }



       counter++;
     }
   }


//cut everything with less than 50 events
//    for(Int_t i=1;i<101;i++){
//     for(Int_t j=1;j<201;j++){
//       if(hBeta_pip_50cut->GetBinContent(i,j)<50)hBeta_pip_50cut->SetBinContent(i,j,0);
//      }
//      }
//
//
//    hBeta_pip_50cut->FitSlicesY();
// //the mean and sigma
//    TH1D *hBeta_pip_50cut_1 = (TH1D*)gDirectory->Get("hBeta_pip_50cut_1");
//    TH1D *hBeta_pip_50cut_2 = (TH1D*)gDirectory->Get("hBeta_pip_50cut_2");
// //to become mean plus and minus 3 sigma
//    TH1D *hBeta_pip_50cut_1_plus_3sig = (TH1D*)hBeta_pip_50cut_1->Clone("hBeta_pip_50cut_1_plus_3sig");
//    TH1D *hBeta_pip_50cut_1_minus_3sig = (TH1D*)hBeta_pip_50cut_1->Clone("hBeta_pip_50cut_1_minus_3sig");
// //To become 3 and -3 Sigma
//    TH1D *hBeta_pip_50cut_2_plus_3sig = (TH1D*)hBeta_pip_50cut_2->Clone("hBeta_pip_50cut_2_plus_3sig");
//    TH1D *hBeta_pip_50cut_2_minus_3sig = (TH1D*)hBeta_pip_50cut_2->Clone("hBeta_pip_50cut_2_minus_3sig");
//
//    hBeta_pip_50cut_2_plus_3sig->Scale(3);
//    hBeta_pip_50cut_2_minus_3sig->Scale(-3);
//
//    hBeta_pip_50cut_1_plus_3sig->Add(hBeta_pip_50cut_2_plus_3sig);
//    hBeta_pip_50cut_1_minus_3sig->Add(hBeta_pip_50cut_2_minus_3sig);
//
//    hBeta_pip_50cut_1->GetXaxis()->SetTitle("P [GeV/c]");
//    hBeta_pip_50cut_1->GetYaxis()->SetTitle("Mean of #Delta#beta");
//    hBeta_pip_50cut_1->SetTitle("Mean #Delta#beta for #pi^{+}");
//    hBeta_pip_50cut_1->GetYaxis()->SetLabelSize(0.05);
//    hBeta_pip_50cut_1->GetYaxis()->SetLabelSize(0.05);
//
//    hBeta_pip_50cut_1_plus_3sig->SetTitle("Mean #Delta#beta plus 3#sigma for #pi^{+}");
//    hBeta_pip_50cut_1_plus_3sig->GetXaxis()->SetTitle("P [GeV/c]");
//    hBeta_pip_50cut_1_plus_3sig->GetYaxis()->SetTitle("Mean plus #sigma of #Delta#beta");
//    hBeta_pip_50cut_1_plus_3sig->GetXaxis()->SetLabelSize(0.05);
//    hBeta_pip_50cut_1_plus_3sig->GetYaxis()->SetLabelSize(0.05);
//
//    hBeta_pip_50cut_1_minus_3sig->SetTitle("Mean #Delta#beta minus 3#sigma for #pi^{+}");
//    hBeta_pip_50cut_1_minus_3sig->GetXaxis()->SetTitle("P [GeV/c]");
//    hBeta_pip_50cut_1_minus_3sig->GetYaxis()->SetTitle("Mean mins #sigma of #Delta#beta");
//    hBeta_pip_50cut_1_minus_3sig->GetXaxis()->SetLabelSize(0.05);
//    hBeta_pip_50cut_1_minus_3sig->GetYaxis()->SetLabelSize(0.05);
//
//    //to become mean plus and minus 2 sigma
//    TH1D *hBeta_pip_50cut_1_plus_2sig = (TH1D*)hBeta_pip_50cut_1->Clone("hBeta_pip_50cut_1_plus_2sig");
//    TH1D *hBeta_pip_50cut_1_minus_2sig = (TH1D*)hBeta_pip_50cut_1->Clone("hBeta_pip_50cut_1_minus_2sig");
//    //To become 2 and -2 Sigma
//    TH1D *hBeta_pip_50cut_2_plus_2sig = (TH1D*)hBeta_pip_50cut_2->Clone("hBeta_pip_50cut_2_plus_2sig");
//    TH1D *hBeta_pip_50cut_2_minus_2sig = (TH1D*)hBeta_pip_50cut_2->Clone("hBeta_pip_50cut_2_minus_2sig");
//
//    hBeta_pip_50cut_2_plus_2sig->Scale(2);
//    hBeta_pip_50cut_2_minus_2sig->Scale(-2);
//
//    hBeta_pip_50cut_1_plus_2sig->Add(hBeta_pip_50cut_2_plus_2sig);
//    hBeta_pip_50cut_1_minus_2sig->Add(hBeta_pip_50cut_2_minus_2sig);
//
//    hBeta_pip_50cut_1_plus_2sig->SetTitle("Mean #Delta#beta plus 2#sigma for #pi^{+}");
//    hBeta_pip_50cut_1_plus_2sig->GetXaxis()->SetTitle("P [GeV/c]");
//    hBeta_pip_50cut_1_plus_2sig->GetYaxis()->SetTitle("Mean plus #sigma of #Delta#beta");
//    hBeta_pip_50cut_1_plus_2sig->GetXaxis()->SetLabelSize(0.05);
//    hBeta_pip_50cut_1_plus_2sig->GetYaxis()->SetLabelSize(0.05);
//
//    hBeta_pip_50cut_1_minus_2sig->SetTitle("Mean #Delta#beta minus 2#sigma for #pi^{+}");
//    hBeta_pip_50cut_1_minus_2sig->GetXaxis()->SetTitle("P [GeV/c]");
//    hBeta_pip_50cut_1_minus_2sig->GetYaxis()->SetTitle("Mean mins #sigma of #Delta#beta");
//    hBeta_pip_50cut_1_minus_2sig->GetXaxis()->SetLabelSize(0.05);
//    hBeta_pip_50cut_1_minus_2sig->GetYaxis()->SetLabelSize(0.05);
//
//    //to become mean plus and minus 1 sigma
//    TH1D *hBeta_pip_50cut_1_plus_sig = (TH1D*)hBeta_pip_50cut_1->Clone("hBeta_pip_50cut_1_plus_sig");
//    TH1D *hBeta_pip_50cut_1_minus_sig = (TH1D*)hBeta_pip_50cut_1->Clone("hBeta_pip_50cut_1_minus_sig");
//    //To become 1 and -1 Sigma
//    TH1D *hBeta_pip_50cut_2_plus_sig = (TH1D*)hBeta_pip_50cut_2->Clone("hBeta_pip_50cut_2_plus_sig");
//    TH1D *hBeta_pip_50cut_2_minus_sig = (TH1D*)hBeta_pip_50cut_2->Clone("hBeta_pip_50cut_2_minus_sig");
//
//    hBeta_pip_50cut_2_minus_sig->Scale(-1);
//
//    hBeta_pip_50cut_1_plus_sig->Add(hBeta_pip_50cut_2_plus_sig);
//    hBeta_pip_50cut_1_minus_sig->Add(hBeta_pip_50cut_2_minus_sig);
//
//    hBeta_pip_50cut_1_plus_sig->SetTitle("Mean #Delta#beta plus #sigma for #pi^{+}");
//    hBeta_pip_50cut_1_plus_sig->GetXaxis()->SetTitle("P [GeV/c]");
//    hBeta_pip_50cut_1_plus_sig->GetYaxis()->SetTitle("Mean plus #sigma of #Delta#beta");
//    hBeta_pip_50cut_1_plus_sig->GetXaxis()->SetLabelSize(0.05);
//    hBeta_pip_50cut_1_plus_sig->GetYaxis()->SetLabelSize(0.05);
//
//    hBeta_pip_50cut_1_minus_sig->SetTitle("Mean #Delta#beta minus #sigma for #pi^{+}");
//    hBeta_pip_50cut_1_minus_sig->GetXaxis()->SetTitle("P [GeV/c]");
//    hBeta_pip_50cut_1_minus_sig->GetYaxis()->SetTitle("Mean mins #sigma of #Delta#beta");
//    hBeta_pip_50cut_1_minus_sig->GetXaxis()->SetLabelSize(0.05);
//    hBeta_pip_50cut_1_minus_sig->GetYaxis()->SetLabelSize(0.05);
//
//    //cut everything with less than 50 events
//       for(Int_t i=1;i<101;i++){
//        for(Int_t j=1;j<201;j++){
//          if(hBeta_pim_100cut->GetBinContent(i,j)<100)hBeta_pim_100cut->SetBinContent(i,j,0);
//         }
//         }
//
//       hBeta_pim_100cut->FitSlicesY();
//    //the mean and sigma
//       TH1D *hBeta_pim_100cut_1 = (TH1D*)gDirectory->Get("hBeta_pim_100cut_1");
//       TH1D *hBeta_pim_100cut_2 = (TH1D*)gDirectory->Get("hBeta_pim_100cut_2");
//    //to become mean plus and minus 3 sigma
//       TH1D *hBeta_pim_100cut_1_plus_3sig = (TH1D*)hBeta_pim_100cut_1->Clone("hBeta_pim_100cut_1_plus_3sig");
//       TH1D *hBeta_pim_100cut_1_minus_3sig = (TH1D*)hBeta_pim_100cut_1->Clone("hBeta_pim_100cut_1_minus_3sig");
//    //To become 3 and -3 Sigma
//       TH1D *hBeta_pim_100cut_2_plus_3sig = (TH1D*)hBeta_pim_100cut_2->Clone("hBeta_pim_100cut_2_plus_3sig");
//       TH1D *hBeta_pim_100cut_2_minus_3sig = (TH1D*)hBeta_pim_100cut_2->Clone("hBeta_pim_100cut_2_minus_3sig");
//
//       hBeta_pim_100cut_2_plus_3sig->Scale(3);
//       hBeta_pim_100cut_2_minus_3sig->Scale(-3);
//
//       hBeta_pim_100cut_1_plus_3sig->Add(hBeta_pim_100cut_2_plus_3sig);
//       hBeta_pim_100cut_1_minus_3sig->Add(hBeta_pim_100cut_2_minus_3sig);
//
//       hBeta_pim_100cut_1->GetXaxis()->SetTitle("P [GeV/c]");
//       hBeta_pim_100cut_1->GetYaxis()->SetTitle("Mean of #Delta#beta");
//       hBeta_pim_100cut_1->SetTitle("Mean #Delta#beta for #pi^{+}");
//       hBeta_pim_100cut_1->GetYaxis()->SetLabelSize(0.05);
//       hBeta_pim_100cut_1->GetYaxis()->SetLabelSize(0.05);
//
//       hBeta_pim_100cut_1_plus_3sig->SetTitle("Mean #Delta#beta plus 3#sigma for #pi^{+}");
//       hBeta_pim_100cut_1_plus_3sig->GetXaxis()->SetTitle("P [GeV/c]");
//       hBeta_pim_100cut_1_plus_3sig->GetYaxis()->SetTitle("Mean plus #sigma of #Delta#beta");
//       hBeta_pim_100cut_1_plus_3sig->GetXaxis()->SetLabelSize(0.05);
//       hBeta_pim_100cut_1_plus_3sig->GetYaxis()->SetLabelSize(0.05);
//
//       hBeta_pim_100cut_1_minus_3sig->SetTitle("Mean #Delta#beta minus 3#sigma for #pi^{+}");
//       hBeta_pim_100cut_1_minus_3sig->GetXaxis()->SetTitle("P [GeV/c]");
//       hBeta_pim_100cut_1_minus_3sig->GetYaxis()->SetTitle("Mean mins #sigma of #Delta#beta");
//       hBeta_pim_100cut_1_minus_3sig->GetXaxis()->SetLabelSize(0.05);
//       hBeta_pim_100cut_1_minus_3sig->GetYaxis()->SetLabelSize(0.05);
//
//       //to become mean plus and minus 2 sigma
//       TH1D *hBeta_pim_100cut_1_plus_2sig = (TH1D*)hBeta_pim_100cut_1->Clone("hBeta_pim_100cut_1_plus_2sig");
//       TH1D *hBeta_pim_100cut_1_minus_2sig = (TH1D*)hBeta_pim_100cut_1->Clone("hBeta_pim_100cut_1_minus_2sig");
//       //To become 2 and -2 Sigma
//       TH1D *hBeta_pim_100cut_2_plus_2sig = (TH1D*)hBeta_pim_100cut_2->Clone("hBeta_pim_100cut_2_plus_2sig");
//       TH1D *hBeta_pim_100cut_2_minus_2sig = (TH1D*)hBeta_pim_100cut_2->Clone("hBeta_pim_100cut_2_minus_2sig");
//
//       hBeta_pim_100cut_2_plus_2sig->Scale(2);
//       hBeta_pim_100cut_2_minus_2sig->Scale(-2);
//
//       hBeta_pim_100cut_1_plus_2sig->Add(hBeta_pim_100cut_2_plus_2sig);
//       hBeta_pim_100cut_1_minus_2sig->Add(hBeta_pim_100cut_2_minus_2sig);
//
//       hBeta_pim_100cut_1_plus_2sig->SetTitle("Mean #Delta#beta plus 2#sigma for #pi^{+}");
//       hBeta_pim_100cut_1_plus_2sig->GetXaxis()->SetTitle("P [GeV/c]");
//       hBeta_pim_100cut_1_plus_2sig->GetYaxis()->SetTitle("Mean plus #sigma of #Delta#beta");
//       hBeta_pim_100cut_1_plus_2sig->GetXaxis()->SetLabelSize(0.05);
//       hBeta_pim_100cut_1_plus_2sig->GetYaxis()->SetLabelSize(0.05);
//
//       hBeta_pim_100cut_1_minus_2sig->SetTitle("Mean #Delta#beta minus 2#sigma for #pi^{+}");
//       hBeta_pim_100cut_1_minus_2sig->GetXaxis()->SetTitle("P [GeV/c]");
//       hBeta_pim_100cut_1_minus_2sig->GetYaxis()->SetTitle("Mean mins #sigma of #Delta#beta");
//       hBeta_pim_100cut_1_minus_2sig->GetXaxis()->SetLabelSize(0.05);
//       hBeta_pim_100cut_1_minus_2sig->GetYaxis()->SetLabelSize(0.05);
//
//       //to become mean plus and minus 1 sigma
//       TH1D *hBeta_pim_100cut_1_plus_sig = (TH1D*)hBeta_pim_100cut_1->Clone("hBeta_pim_100cut_1_plus_sig");
//       TH1D *hBeta_pim_100cut_1_minus_sig = (TH1D*)hBeta_pim_100cut_1->Clone("hBeta_pim_100cut_1_minus_sig");
//       //To become 1 and -1 Sigma
//       TH1D *hBeta_pim_100cut_2_plus_sig = (TH1D*)hBeta_pim_100cut_2->Clone("hBeta_pim_100cut_2_plus_sig");
//       TH1D *hBeta_pim_100cut_2_minus_sig = (TH1D*)hBeta_pim_100cut_2->Clone("hBeta_pim_100cut_2_minus_sig");
//
//       hBeta_pim_100cut_2_minus_sig->Scale(-1);
//
//       hBeta_pim_100cut_1_plus_sig->Add(hBeta_pim_100cut_2_plus_sig);
//       hBeta_pim_100cut_1_minus_sig->Add(hBeta_pim_100cut_2_minus_sig);
//
//       hBeta_pim_100cut_1_plus_sig->SetTitle("Mean #Delta#beta plus #sigma for #pi^{+}");
//       hBeta_pim_100cut_1_plus_sig->GetXaxis()->SetTitle("P [GeV/c]");
//       hBeta_pim_100cut_1_plus_sig->GetYaxis()->SetTitle("Mean plus #sigma of #Delta#beta");
//       hBeta_pim_100cut_1_plus_sig->GetXaxis()->SetLabelSize(0.05);
//       hBeta_pim_100cut_1_plus_sig->GetYaxis()->SetLabelSize(0.05);
//
//       hBeta_pim_100cut_1_minus_sig->SetTitle("Mean #Delta#beta minus #sigma for #pi^{+}");
//       hBeta_pim_100cut_1_minus_sig->GetXaxis()->SetTitle("P [GeV/c]");
//       hBeta_pim_100cut_1_minus_sig->GetYaxis()->SetTitle("Mean mins #sigma of #Delta#beta");
//       hBeta_pim_100cut_1_minus_sig->GetXaxis()->SetLabelSize(0.05);
//       hBeta_pim_100cut_1_minus_sig->GetYaxis()->SetLabelSize(0.05);

   gBenchmark->Stop("timer");
   gBenchmark->Print("timer");

   fileOutput1.Write();
/*
   hBeta_1->Write();
   hBeta_1_copy1->Write();
   hBeta_1_copy2->Write();
*/





/*
   TCanvas* can2=new TCanvas();
   can2->Divide(2,2);
   can2->cd(1);
   hBeta_1->DrawCopy();
   can2->cd(2);
   hBeta_2->DrawCopy();
   can2->cd(3);
   hBeta_1_copy1->DrawCopy();
   can2->cd(4);
   hBeta_1_copy2->DrawCopy();

*/
/*
   can2->Divide(2,2);
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
*/
   auto finish = std::chrono::high_resolution_clock::now();
   std::chrono::duration<double> elapsed = finish - start;
   std::cout << "Elapsed time: " << elapsed.count()<< " events = "<<counter<< " s\n";

}
