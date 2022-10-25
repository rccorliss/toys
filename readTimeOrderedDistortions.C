//basic framework to read and compare events in a time ordered distortion file
//these files contain a TTree with xingnum, and three distortion map branches
#include "TTree.h" //this prevents a lazy binding issue and/or is a magic spell.
#include "TCanvas.h" //this prevents a lazy binding issue and/or is a magic spell.
#include "TH3F.h"
#include "TH2F.h"
#include "TH1F.h"

void readTimeOrderedDistortions(char *filename="TimeOrderedDistortions.root"){

  TFile *file=TFile::Open(filename);
  TTree *tree=(TTree*)(file->Get("TimeDists"));

  TH3F *basehist[3];
  TH3F *temphist[3];

  //don't forget you have to make those histogram objects.  It assumes they exist!
  // std::string histname[]={"hIntDistortionZ","hIntDistortionY","hIntDistortionX"};
  // std::string histname[]={"IntDistortionZ_negz","IntDistortionR_negz","IntDistortionP_negz"};
    std::string histname[]={"hIntDistortionZ_negz","hIntDistortionR_negz","hIntDistortionP_negz"};

  for (int i=0;i<3;i++){
    basehist[i]=new TH3F(Form("basehist%d",i),Form("basehist%d",i),10,0,10,20,0,20,30,0,30);
    temphist[i]=new TH3F(Form("temphist%d",i),Form("temphist%d",i),10,0,10,20,0,20,30,0,30);
    tree->SetBranchAddress(histname[i].c_str(),&(temphist[i]));
  }
  printf("histograms built.\n");
// tree->SetBranchAddress("hIntDistortionZ",&(temphist[0]));
//  tree->SetBranchAddress("hIntDistortionY",&(temphist[1]));
//  tree->SetBranchAddress("hIntDistortionX",&(temphist[2]));
  printf("branches set.\n");

  //keep the last histogram of the series as a reference histogram.
  tree->GetEntry(tree->GetEntries()-1);
  float max[3];
  float min[3];
  float edge[3];
  for (int i=0;i<3;i++){
    basehist[i]=(TH3F*)(temphist[i]->Clone());
    min[i]=fabs(basehist[i]->GetBinContent(basehist[i]->GetMinimumBin()));
    max[i]=fabs(basehist[i]->GetBinContent(basehist[i]->GetMaximumBin()));
    edge[i]=2*(min[i]>max[i]?min[i]:max[i]);
  }
  printf("base histograms set.\n");

  
  //now loop over the others and see what the subtraction looks like:
  TH2F *baseVsTempHist[3];
  TH1F *baseRatioHist[3];
  TH1F *histProfile[3];
  
  for (int i=0;i<3;i++){
    baseVsTempHist[i]=new TH2F(Form("baseVsTempHist%d",i),Form("bin(i,xingj) vs bin(i,xing0)(from hist%d)",i),100,-1*edge[i],edge[i],100,-1*edge[i],edge[i]);
    baseRatioHist[i]=new TH1F(Form("baseRatioHist%d",i),Form("Ratio of bin(i,xingj)/bin(i,xing0)(from hist%d)",i),200,-2,3);
    histProfile[i]=new TH1F(Form("histProfile%d",i),Form(" bin(i,xingj)-bin(i,xing0) contents of hist %d across all events",i),500,-2*edge[i],2*edge[i]);
    //histProfile[i]=new TH1F(Form("histProfile%d",i),Form("bin contents of hist %d across all events",i),500,-edge[i],edge[i]);

  }
  printf("comparison histograms built.\n");

  printf("tree has %lld entries.\n", tree->GetEntries());

  for (int i=0;i<tree->GetEntries();i++){
    printf("  working on entry %d\n",i);
    tree->GetEntry(i);
    //loop over all hists and compare:
    for (int j=0;j<3;j++){
      int nbins=temphist[j]->GetNcells();
      printf("  working on entry %d hist %d (%d cells)\n",i,j, nbins);
      for (int k=0;k<nbins;k++){
	float b=basehist[j]->GetBinContent(k);
	float t=temphist[j]->GetBinContent(k);
	baseVsTempHist[j]->Fill(b,t);
	histProfile[j]->Fill(t-b);
	if (b!=0) {
	  baseRatioHist[j]->Fill(t/b);
	} else {
	  baseRatioHist[j]->Fill(-1);
	}
      }
      printf("  done with entry %d hist %d (%d cells)\n",i,j, nbins);

    }
  }
  printf("comparison histograms filled.\n");

  TCanvas* c=new TCanvas();
  c->Divide(3,3);
  for (int i=0;i<3;i++){
    c->cd(1+i);
    baseVsTempHist[i]->Draw("colz");
  }
   for (int i=0;i<3;i++){
    c->cd(4+i);
    baseRatioHist[i]->Draw("colz");
  }
   for (int i=0;i<3;i++){
     c->cd(7+i)->SetLogy();
    histProfile[i]->Draw("colz");
  }

   return;
}

  
  
  
