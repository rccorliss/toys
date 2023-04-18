//basic framework to read in a time ordered set of distortions and randomly draw from it
//to generate a randomly shuffled set of distortions
//also can randomly rotate individual events if desired.and compare events in a time ordered distortion file
//these files contain a TTree with xingnum, and three distortion map branches
#include "TTree.h" //this prevents a lazy binding issue and/or is a magic spell.
#include "TCanvas.h" //this prevents a lazy binding issue and/or is a magic spell.
#include "TH3F.h"
#include "TH2F.h"
#include "TH1F.h"

void shiftPhiWithGuardBins(TH3 *hist, int phiShift);


void shuffleDistortions(char *inputFileName="TimeOrderedDistortions.root", char *outputFileName="temp.hist.root", int nRequested=1, bool doRotation=false){


  printf("loading %s.  Drawing %d random events%s and outputing to %s\n", inputFileName,nRequested,doRotation?", rotating them,":"",outputFileName);
  TFile *file=TFile::Open(inputFileName,"READ");
  TTree *tree=(TTree*)(file->Get("TimeDists"));

  int nhists=6;
  TH3F *inhist[nhists];
  TH3F *outhist[nhists];

  
 std::string branchname[]={"hIntDistortionP_negz",
			    "hIntDistortionPhi_negz",
			    "hIntDistortionR_negz",
			    "hIntDistortionZ_negz",
			    "hIntDistortionP_posz",
			    "hIntDistortionPhi_posz",
			    "hIntDistortionR_posz",
			    "hIntDistortionZ_posz"};
  

  TFile *outfile=TFile::Open(outputFileName,"RECREATE");
  TTree *outtree=new TTree("TimeDists","TimeDists");

  int xingnum;
  outtree->Branch("xingnum",&xingnum);
 for (int i=0;i<nhists;i++){
    inhist[i]=new TH3F(Form("inhist%d",i),Form("inhist%d",i),10,0,10,20,0,20,30,0,30);
    outhist[i]=new TH3F(Form("outhist%d",i),Form("outhist%d",i),10,0,10,20,0,20,30,0,30);
    tree->SetBranchAddress(branchname[i].c_str(),&(inhist[i]));
    //outtree->Branch(branchname[i].c_str(),&(outhist[i]));
  }
  printf("histograms built and branches set.\n");

  TRandom *rand=new TRandom3();
  

  for (int i=0;i<nRequested;i++){
    //    if ((i)%(nRequested/10)==0) printf(".");
    //select a random event:
    int event=rand->Integer(tree->GetEntries());
    printf("i=%d getting event %d of %d\n",i,event,tree->GetEntries());

    tree->GetEntry(event);
    xingnum=i;
    for (int j=0;j<nhists;j++){
      outhist[j]=inhist[j];//probably unnecessary to clone.
    }
    if (doRotation){
      int phiShift=rand->Integer(outhist[0]->GetXaxis()->GetNbins()-2);
      for (int j=0;j<nhists;j++){
	shiftPhiWithGuardBins(outhist[j],phiShift);
      }
    }
    outtree->Fill();
  }

   outtree->Write();
  outfile->Close();
  return;
}


  
  
void shiftPhiWithGuardBins(TH3 *hist, int phiShift){

  //get dimensions:
  TAxis *ax[3]={nullptr,nullptr,nullptr};
  ax[0]=hist->GetXaxis();
  ax[1]=hist->GetYaxis();
  ax[2]=hist->GetZaxis();

  int nbins[3];
  for (int i=0;i<3;i++){
    nbins[i]=ax[i]->GetNbins();//number of bins, not counting under and overflow.
  }
  int nphi=nbins[0]-2;
  float orig[nphi];

  //remember the structure of these histograms, for n 'real' bins:
  //the number of bins reported on the axis will be n+2
  //bin#:         0          1      2     3 ... n   | n+1 |  n+2   |   n+3
  //content: [ underflow | guard |  0  |  1 ... n-2 | n-1 | guard | overflow  ]
  
  for (int r=1;r<=nbins[1];r++){//include guards
    for (int z=1;z<=nbins[2];z++){//include guards
      //get the original values so we can overwrite them:
      for (int p=0;p<nphi;p++){//exclude guards
	orig[p]=hist->GetBinContent(hist->GetBin(2+p,r,z));
      }
      //now fill them with our orig data, rotated:
      for (int p=0;p<nphi;p++){//exclude guards
	hist->SetBinContent(hist->GetBin(2+p,r,z),orig[(p+phiShift)%nphi]);
      }      
      //now match the guard bins up in phi:
      //low guard matches the highest phi bin
      hist->SetBinContent(hist->GetBin(1,r,z),
			  hist->GetBinContent(hist->GetBin(nphi+1,r,z)));
      //and high guard matches the lowest phi bin
      hist->SetBinContent(hist->GetBin(nphi+2,r,z),
			  hist->GetBinContent(hist->GetBin(2,r,z)));
    }
  }

  return;
}

