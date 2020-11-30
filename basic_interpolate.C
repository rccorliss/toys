#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"

float PI=3.1415927;
void basic_interpolate(){
  const float rmin=24.;
  const float activermin=30.;//don't bother plotting averages below this range
  const float rmax=78.;
  float rspan=rmax-rmin;
  int nrsteps=120;
  float rstep=rspan/(1.0*nrsteps);
  const float zmin=0;
  const float zmax=105.5;
  float zspan=zmax-zmin;
  int nzsteps=120;
  float zstep=zspan/(1.0*nzsteps);

  //for normalization studies
  float phimin=0.0;
  float phimax=6.28;
  float phispan=phimax-phimin;
  int nphisteps=60;
  float phistep=phispan/(1.0*nphisteps);


  float IFCrange[]={0.,15.,15.};
  float OFCrange[]={25.,75.,105.};
  float boundM[3];
  float boundB[3];
  for (int i=0;i<3;i++){
    boundM[i]=(OFCrange[i]-IFCrange[i])/(rmax-rmin);
    boundB[i]=IFCrange[i]-rmin*boundM[i];
  }
  
  //load the time-average histogram(s):
  TFile *histfile;
  histfile=TFile::Open("elevatorpitch/fluct_average.rev3.1side.3d.file0.h_negz.real_B1.4_E-400.0.ross_phi1_sphenix_phislice_lookup_r26xp40xz40.distortion_map.hist.root","READ");
  //histfile=TFile::Open("elevatorpitch/average.rev3.1side.3d.file0.h_negz.real_B1.4_E-400.0.ross_phi1_sphenix_phislice_lookup_r26xp40xz40.distortion_map.hist.root","READ");
  TH3F *hData[3];
  hData[0]=(TH3F*)histfile->Get("hIntDistortionR");
  hData[1]=(TH3F*)histfile->Get("hIntDistortionP");
  hData[2]=(TH3F*)histfile->Get("hIntDistortionZ");

  printf("data loaded.\n");

  
  //define the space between coverage regions
  //25cm wide at outer radius of approx 80cm, 12 such tiles
  float deltaPhi=2*PI/12.0;
  float phiplane[3];
  char dirchar[]="rpz";
  phiplane[0]=1;
  phiplane[2]=phiplane[0]+deltaPhi;
  phiplane[1]=0.5*(phiplane[0]+phiplane[2]);
  TH2F *hSlice[3][3];
  for (int i=0;i<3;i++){
    for (int j=0;j<3;j++){
      hSlice[i][j]=new TH2F(Form("hSlice%d%d",i,j),Form("true %c-distortion in r-z plane at phi=%2.2f;z (cm);r (cm)",dirchar[j],phiplane[i]),nzsteps,zmin,zmax,nrsteps,rmin,rmax);
    }
  }
  
  for (int ir=0;ir<nrsteps;ir++){
    float r=rstep*(ir+0.5)+rmin;
    for (int iz=0;iz<nzsteps;iz++){
      float z=zstep*(iz+0.501)+zmin;
       if (iz!=0){
	if ((z<boundM[0]*r+boundB[0])||(z>boundM[1]*r+boundB[1])){
	  continue;
	}
      }
      for (int i=0;i<3;i++){
	for (int j=0;j<3;j++){
	  hSlice[i][j]->Fill(z,r,hData[j]->Interpolate(phiplane[i],r,z));
	}
      }
    }
  }

  //do the cheapest interpolation:
  TH2F *hGuess[3];
  TH2F *hDiff2D[3];
  TH1F *hDiff[3];
  TH2F *hBetterGuess[3];
  TH2F *hBetterDiff2D[3];
  TH1F *hBetterDiff[3];
  for (int j=0;j<3;j++){
    printf("j=%d\n",j);
    hGuess[j]=new TH2F(Form("hGuess%d",j),Form("interpolated %c-distortion in r-z plane at phi=%2.2f;z (cm);r (cm)",dirchar[j],phiplane[1]),nzsteps,zmin,zmax,nrsteps,rmin,rmax);
    hGuess[j]->Add(hSlice[0][j],hSlice[2][j],0.5,0.5);
    hDiff[j]=new TH1F(Form("hDiff%d",j),Form("%c-distortion Difference Between INterpolation and true (r>%1.1f);diff (um)",dirchar[j],activermin),200,-400,400);
    hDiff2D[j]=new TH2F(Form("hDiff2D%d",j),Form("%c-distortion Difference Between INterpolation and true (um);z (cm);r (cm)",dirchar[j]),nzsteps,zmin,zmax,nrsteps,rmin,rmax);
    hDiff2D[j]->Add(hSlice[1][j],hGuess[j],-1.0,1.0);
    printf("j=%d done\n",j);
    
    hBetterGuess[j]=new TH2F(Form("hBetterGuess%d",j),Form("interpolated %c-distortion with CM correction, phi=%2.2f;z (cm);r (cm)",dirchar[j],phiplane[1]),nzsteps,zmin,zmax,nrsteps,rmin,rmax);
    for (int ir=0;ir<nrsteps;ir++){
      float r=rstep*(ir+0.5)+rmin;
      for (int iz=0;iz<nzsteps;iz++){
	float z=zstep*(iz+0.501)+zmin;
      if (iz!=0){
	if ((z<boundM[0]*r+boundB[0])||(z>boundM[1]*r+boundB[1])){
	  continue;
	}
      }
      int cmbin=hGuess[0]->FindBin(zmin+zstep*0.501,r);
	int bin=hGuess[0]->FindBin(z,r);
	float cmguess=0.5*(hSlice[0][j]->GetBinContent(cmbin)+hSlice[2][j]->GetBinContent(cmbin));//average of CM values in two measured points
	float cmtrue=hSlice[1][j]->GetBinContent(cmbin);//true value at CM
	float rescale=cmtrue/cmguess;
	hBetterGuess[j]->Fill(z,r,hGuess[j]->GetBinContent(bin)*rescale);
      }
    }
    hBetterDiff[j]=new TH1F(Form("hBetterDiff%d",j),Form("%c-distortion Difference Between int++ and true (r>%1.1f);diff (um)",dirchar[j],activermin),200,-400,400);


    hBetterDiff2D[j]=new TH2F(Form("hBetterDiff2D%d",j),Form("%c-distortion Difference Between int++ and true (um) ;z (cm);r (cm)",dirchar[j]),nzsteps,zmin,zmax,nrsteps,rmin,rmax);
    hBetterDiff2D[j]->Add(hSlice[1][j],hBetterGuess[j],-1.0,1.0);

    hBetterDiff2D[j]->Scale(1e4);
    //hBetterDiff[j]->Scale(1e4);
    hDiff2D[j]->Scale(1e4);
    //hDiff[j]->Scale(1e4);
    
  }

  for (int ir=0;ir<nrsteps;ir++){
    float r=rstep*(ir+0.5)+rmin;
    for (int iz=0;iz<nzsteps;iz++){
      float z=zstep*(iz+0.501)+zmin;
      if (iz!=0){
	if ((z<boundM[0]*r+boundB[0])||(z>boundM[1]*r+boundB[1])){
	  continue;
	}
      }
      int bin=hGuess[0]->FindBin(z,r);
      if (r>activermin){
	for (int j=0;j<3;j++){
	  hDiff[j]->Fill(1e4*(hGuess[j]->GetBinContent(bin)-hSlice[1][j]->GetBinContent(bin)));
	  hBetterDiff[j]->Fill(1e4*(hBetterGuess[j]->GetBinContent(bin)-hSlice[1][j]->GetBinContent(bin)));
	}
      }
    }
  } //read the true x,y,z distortions at phi[0] and phi[1]
  //create a linear function for all intermediate phi
  //compare the values in these intermediate regions to the true hist values there



  //generate the true slices over the forward z converage for two phi regions covered by OTs.
  //slice 0 represents the slice with full z coverage, slice 1 the one with only partial
   TH2F *hFullSlice[2][3];
  for (int i=0;i<2;i++){
    for (int j=0;j<3;j++){
      hFullSlice[i][j]=new TH2F(Form("hFullSlice%d%d",i,j),Form("true %c-distortion in r-z plane at phi=%2.2f;z (cm);r (cm)",dirchar[j],phiplane[i]),nzsteps,zmin,zmax,nrsteps,rmin,rmax);
    }
  }
  
  for (int ir=0;ir<nrsteps;ir++){
    float r=rstep*(ir+0.5)+rmin;
    for (int iz=0;iz<nzsteps;iz++){
      float z=zstep*(iz+0.501)+zmin;
       if (iz!=0){
	if ((z<boundM[0]*r+boundB[0])||(z>boundM[2]*r+boundB[2])){
	  continue;
	}
      }
      for (int i=0;i<2;i++){
	for (int j=0;j<3;j++){
	  hFullSlice[i][j]->Fill(z,r,hData[j]->Interpolate(phiplane[i*2],r,z));
	}
      }
    }
  }

  //make three guesses:
  //0: Exact copy of phi0 plane
  //1: phi0 plane scaled at each r by the cm ratio
  //2: phi0 plane scaled at each r by the ratio at the coverage edge
  TH2F *hForwardGuess[3][3];
  TH2F *hForwardDiff2D[3][3];
  TH1F *hForwardDiff[3][3];
  TString guessName[]={"copy of phiplane0","phiplane0*CMscale","phiplane0*midZscale"};
  for (int i=0;i<3;i++){//guess type
    for (int j=0;j<3;j++){//direction
      hForwardGuess[i][j]=new TH2F(Form("hForwardGuess%d%d",i,j),Form("%c-distortion in r-z plane at phi=%2.2f using %s;z (cm);r (cm)",dirchar[j],phiplane[2],guessName[i].Data()),nzsteps,zmin,zmax,nrsteps,rmin,rmax);
      hForwardDiff[i][j]=new TH1F(Form("hForwardDiff%d%d",i,j),Form("(Guess%d-True) for %c-hat (r>%1.1f);diff (um)",i,dirchar[j],activermin),200,-150,150);
      hForwardDiff2D[i][j]=new TH2F(Form("hForwardDiff2D%d%d",i,j),Form("(Guess%d-True) for %c-hat (um);z (cm);r (cm)",i,dirchar[j]),nzsteps,zmin,zmax,nrsteps,rmin,rmax);
    }
  }


  float cmrescale[3];
  float midZrescale[3];
  for (int ir=0;ir<nrsteps;ir++){
    float r=rstep*(ir+0.5)+rmin;
    for (int iz=0;iz<nzsteps;iz++){
      float z=zstep*(iz+0.501)+zmin;
      //calculate the CM rescale for this r slice:
      if (iz==0){
	for (int j=0;j<3;j++){
	  int cmbin=hFullSlice[0][j]->FindBin(zmin+zstep*0.501,r);
	  float cmsource=hFullSlice[0][j]->GetBinContent(cmbin);
	  float cmdestination=hFullSlice[1][j]->GetBinContent(cmbin);
	  cmrescale[j]=cmdestination/cmsource;
	  printf("cmrescale%d=%f=(%f/%f)\n",j,cmrescale[j],cmdestination,cmsource);
	}
      }
      //calculate the midZ rescale on the last bin before we're in the forward region
      if (z<(boundM[1]*r+boundB[1]) && (boundM[1]*r+boundB[1]-z)<zstep){
	for (int j=0;j<3;j++){
	  int midzbin=hFullSlice[0][0]->FindBin(z,r);
	  float midzsource=hFullSlice[0][j]->GetBinContent(midzbin);
	  float midzdestination=hFullSlice[1][j]->GetBinContent(midzbin);
	  midZrescale[j]=midzdestination/midzsource;
	}
      }    
      if ((z<boundM[1]*r+boundB[1])||(z>boundM[2]*r+boundB[2])){
	continue;
      }
      //printf("rz=%d,%d\n",ir,iz);

      int bin=hFullSlice[0][0]->FindBin(z,r);
      //compute the guesses:
      for (int j=0;j<3;j++){
	float dist=hFullSlice[0][j]->GetBinContent(bin);
	//printf("%c-distortion=%f\n",dirchar[j],dist);
	//for guess 0, copy all cells from FullSlice that are in the z bounds of the forward segment:
	hForwardGuess[0][j]->Fill(z,r,dist);
	//for guess 1, scale the raw dist by the CM ratio,
	hForwardGuess[1][j]->Fill(z,r,dist*cmrescale[j]);
	//for guess 2, scale the raw dist by the midz ratio, 
	hForwardGuess[2][j]->Fill(z,r,dist*midZrescale[j]);
      }

      //generate the comparisons for the guesses (only covering the regions of guessed in, of course)
      for (int i=0;i<3;i++){//guess
	for (int j=0;j<3;j++){//direction
	  float diff=hForwardGuess[i][j]->GetBinContent(bin)-hFullSlice[1][j]->GetBinContent(bin);
	  if (r>activermin){//only plot the 1D hist for the active tpc region
	    hForwardDiff[i][j]->Fill(diff*1e4);
	  }
	  hForwardDiff2D[i][j]->Fill(z,r,diff);
	}
      }
    }
  }
  //rescale the diffs to be in um instead of cm:
  for (int i=0;i<3;i++){//guess
    for (int j=0;j<3;j++){//direction 
      //hForwardDiff[i][j]->Scale(1e4);
      hForwardDiff2D[i][j]->Scale(1e4);
    }
  }


  
 TCanvas *c=new TCanvas("c","c",1900,800);
 c->Divide(6,3);
 for (int j=0;j<3;j++){//direction
   /*
   c->cd(1+j*6); 
   hSlice[0][j]->Draw("colz");
   c->cd(1+j*6+1); 
   hSlice[2][j]->Draw("colz");
   c->cd(1+j*6+2); 
   hSlice[1][j]->Draw("colz");
   */
   c->cd(1+j*6+0); 
   hGuess[j]->Draw("colz");
   c->cd(1+j*6+1); 
   hDiff2D[j]->Draw("colz");
   c->cd(1+j*6+2); 
   hDiff[j]->Draw("hist");
   c->cd(1+j*6+3); 
   hBetterGuess[j]->Draw("colz");
   c->cd(1+j*6+4); 
   hBetterDiff2D[j]->Draw("colz");
   c->cd(1+j*6+5); 
   hBetterDiff[j]->Draw("hist");
 }

 c=new TCanvas("cforward","cforward",1900,800);
 c->Divide(6,3);
 for (int j=0;j<3;j++){//direction
   for (int i=0;i<3;i++){//guess
     c->cd(1+j*6+i*2);
     //hForwardGuess[i][j]->Draw("colz");
     hForwardDiff2D[i][j]->Draw("colz");
     c->cd(2+j*6+i*2); 
     hForwardDiff[i][j]->Draw("hist");
   }
 }
 return;

 //unroll an r-slice and normalize it by the value at the cm:
 float rplane[]={35,50,75};
  TH2F *hRslice[3][3];
  TH2F *hRsliceNorm[3][3];
  TH2F *hRsliceProfile[3][3];
  TH2F *hRsliceNormProfile[3][3];
  for (int i=0;i<3;i++){
    for (int j=0;j<3;j++){
      hRslice[i][j]=new TH2F(Form("hRslice%d%d",i,j),Form("true %c-distortion in phi-z plane at r=%2.2f;phi (rad);z (cm)",dirchar[j],rplane[i]),nphisteps,phimin,phimax,nzsteps,zmin,zmax);
      hRsliceNorm[i][j]=new TH2F(Form("hRsliceNorm%d%d",i,j),Form("CM-Normalized %c-distortion in phi-z plane at r=%2.2f;phi (rad);z (cm)",dirchar[j],rplane[i]),nphisteps,phimin,phimax,nzsteps,zmin,zmax);
      hRsliceProfile[i][j]=new TH2F(Form("hRsliceProfile%d%d",i,j),Form("true %c-distortion at r=%2.2f (all phi);z (cm);distortion (um)",dirchar[j],rplane[i]),nzsteps,zmin,zmax,200,-200,200);
      hRsliceNormProfile[i][j]=new TH2F(Form("hRsliceNormProfile%d%d",i,j),Form("CM-Normalized %c-distortion at r=%2.2f (all phi);z (cm);distortion (um)",dirchar[j],rplane[i]),nzsteps,zmin,zmax,200,-200,200);
    }
  }

  float distortionatzero[nzsteps];
  float distortionatzeroNorm[nzsteps];
  for (int i=0;i<3;i++){//rplane
    for (int j=0;j<3;j++){//direction
      for (int iphi=0;iphi<nphisteps;iphi++){
	float phi=phistep*(iphi+0.5)+phimin;
	float cmnorm=hData[j]->Interpolate(phi,rplane[i],zmin+zstep*0.5);
	for (int iz=0;iz<nzsteps;iz++){
	  float z=zstep*(iz+0.5)+zmin;
	  float dist=hData[j]->Interpolate(phi,rplane[i],z);
	  if (iphi==0) {
	    distortionatzero[iz]=dist;
	    distortionatzeroNorm[iz]=dist/cmnorm;
	  }
	  hRslice[i][j]->Fill(phi,z,dist);
	  hRsliceNorm[i][j]->Fill(phi,z,dist/cmnorm);
	  hRsliceProfile[i][j]->Fill(z,1e4*(dist-distortionatzero[iz]));
	  hRsliceNormProfile[i][j]->Fill(z,1e4*(dist/cmnorm-distortionatzeroNorm[iz]));
	}
      }
    }
  }

c=new TCanvas("c1","c1",1900,800);
 c->Divide(6,3);
 for (int j=0;j<3;j++){//direction
   c->cd(1+j*6); 
   hRslice[0][j]->Draw("colz");
   c->cd(1+j*6+1); 
   hRsliceNorm[0][j]->Draw("colz");
   c->cd(1+j*6+2); 
   hRslice[1][j]->Draw("colz");
   c->cd(1+j*6+3); 
   hRsliceNorm[1][j]->Draw("colz");
   c->cd(1+j*6+4); 
   hRslice[2][j]->Draw("colz");
   c->cd(1+j*6+5); 
   hRsliceNorm[2][j]->Draw("colz");
 }

c=new TCanvas("c2","c2",1900,800);
 c->Divide(6,3);
 for (int j=0;j<3;j++){//direction
   c->cd(1+j*6); 
   hRsliceProfile[0][j]->Draw("colz");
   c->cd(1+j*6+1); 
   hRsliceNormProfile[0][j]->Draw("colz");
   c->cd(1+j*6+2); 
   hRsliceProfile[1][j]->Draw("colz");
   c->cd(1+j*6+3); 
   hRsliceNormProfile[1][j]->Draw("colz");
   c->cd(1+j*6+4); 
   hRsliceProfile[2][j]->Draw("colz");
   c->cd(1+j*6+5); 
   hRsliceNormProfile[2][j]->Draw("colz");
 }
 
   //do the same
return;
}

/*
void PlotExtrapolationToFullZ(){
  
  //define the space between coverage regions
  //25cm wide at outer radius of approx 80cm, 12 such tiles
  float deltaPhi=2*PI/12.0;
  float phiplane[2];
  char dirchar[]="rpz";
  phiplane[0]=1;
  phiplane[1]=phiplane[0]+deltaPhi;
  int nslices=2;
  TH2F *hSlice[nslices][3];
  TH2F *hKnown[3];
  for (int i=0;i<nslices;i++){
    for (int j=0;j<3;j++){
      hSlice[i][j]=new TH2F(Form("hSlice%d%d",i,j),Form("true %c-distortion in r-z plane at phi=%2.2f;z (cm);r (cm)",dirchar[j],phiplane[i]),nzsteps,zmin,zmax,nrsteps,rmin,rmax);
      if (i==0)
	hKnown[j]=new TH2F(Form("hKnown%d",j),Form("known %c-distortion in r-z plane at phi=%2.2f;z (cm);r (cm)",dirchar[j],phiplane[1]),nzsteps,zmin,zmax,nrsteps,rmin,rmax);
    }
  }
  
  for (int ir=0;ir<nrsteps;ir++){
    float r=rstep*(ir+0.5)+rmin;
    for (int iz=0;iz<nzsteps;iz++){
      float z=zstep*(iz+0.501)+zmin;
       if (iz!=0){
	if ((z<boundM[0]*r+boundB[0])||(z>boundM[2]*r+boundB[2])){
	  continue;
	}
      }
      for (int i=0;i<nslices;i++){
	for (int j=0;j<3;j++){
	  hSlice[i][j]->Fill(z,r,hData[j]->Interpolate(phiplane[i],r,z));
	}
      }
      
      if (z>boundM[2]*r+boundB[2]){//in knowable region for test point
	for (int j=0;j<3;j++){
	  hKnown[j]->Fill(z,r,hData[j]->Interpolate(phiplane[1],r,z));
      }
      }
    }
  }

  //do the cheapest extrapolation:
  TH2F *hGuessFull[3];
  TH2F *hDiff2DFull[3];
  TH1F *hDiffFull[3];
  TH2F *hBetterGuessFull[3];
  TH2F *hBetterDiff2DFull[3];
  TH1F *hBetterDiffFull[3];
  for (int j=0;j<3;j++){
    printf("j=%d\n",j);
    hGuess[j]=new TH2F(Form("hGuess%d",j),Form("extrapolated %c-distortion in r-z plane at phi=%2.2f;z (cm);r (cm)",dirchar[j],phiplane[1]),nzsteps,zmin,zmax,nrsteps,rmin,rmax);
    hGuess[j]->Add(hSlice[0][j],hSlice[2][j],0.5,0.5);
    hDiff[j]=new TH1F(Form("hDiff%d",j),Form("%c-distortion Difference Between INterpolation and true;diff (um)",dirchar[j]),200,-400,400);
    hDiff2D[j]=new TH2F(Form("hDiff2D%d",j),Form("%c-distortion Difference Between INterpolation and true (um);z (cm);r (cm)",dirchar[j]),nzsteps,zmin,zmax,nrsteps,rmin,rmax);
    hDiff2D[j]->Add(hSlice[1][j],hGuess[j],-1.0,1.0);
    printf("j=%d done\n",j);
    
    hBetterGuess[j]=new TH2F(Form("hBetterGuess%d",j),Form("interpolated %c-distortion with CM correction, phi=%2.2f;z (cm);r (cm)",dirchar[j],phiplane[1]),nzsteps,zmin,zmax,nrsteps,rmin,rmax);
    for (int ir=0;ir<nrsteps;ir++){
      float r=rstep*(ir+0.5)+rmin;
      for (int iz=0;iz<nzsteps;iz++){
	float z=zstep*(iz+0.501)+zmin;
      if (iz!=0){
	if ((z<boundM[0]*r+boundB[0])||(z>boundM[1]*r+boundB[1])){
	  continue;
	}
      }
	int bin=hGuess[0]->FindBin(z,r);
	float cmguess=0.5*(hSlice[0][j]->GetBinContent(bin)+hSlice[2][j]->GetBinContent(bin));//average of CM values in two measured points
	float cmtrue=hSlice[2][j]->GetBinContent(bin);//true value at CM
	float rescale=cmtrue/cmguess;
	hBetterGuess[j]->Fill(z,r,hGuess[j]->GetBinContent(bin)*rescale);
      }
    }
    hBetterDiff[j]=new TH1F(Form("hBetterDiff%d",j),Form("%c-distortion Difference Between int++ and true;diff (um)",dirchar[j]),200,-400,400);


    hBetterDiff2D[j]=new TH2F(Form("hBetterDiff2D%d",j),Form("%c-distortion Difference Between int++ and true (um) ;z (cm);r (cm)",dirchar[j]),nzsteps,zmin,zmax,nrsteps,rmin,rmax);
    hBetterDiff2D[j]->Add(hSlice[1][j],hBetterGuess[j],-1.0,1.0);

    hBetterDiff2D[j]->Scale(1e4);
    hBetterDiff[j]->Scale(1e4);
    hDiff2D[j]->Scale(1e4);
    hDiff[j]->Scale(1e4);
    
  }

  for (int ir=0;ir<nrsteps;ir++){
    float r=rstep*(ir+0.5)+rmin;
    for (int iz=0;iz<nzsteps;iz++){
      float z=zstep*(iz+0.501)+zmin;
      if (iz!=0){
	if ((z<boundM[0]*r+boundB[0])||(z>boundM[1]*r+boundB[1])){
	  continue;
	}
      }
      int bin=hGuess[0]->FindBin(z,r);
      for (int j=0;j<3;j++){
	hDiff[j]->Fill(1e4*(hGuess[j]->GetBinContent(bin)-hSlice[1][j]->GetBinContent(bin)));
	hBetterDiff[j]->Fill(1e4*(hBetterGuess[j]->GetBinContent(bin)-hSlice[1][j]->GetBinContent(bin)));
      }
    }
  } //read the true x,y,z distortions at phi[0] and phi[1]
  //create a linear function for all intermediate phi
  //compare the values in these intermediate regions to the true hist values there
 TCanvas *c=new TCanvas("c","c",1900,800);
 c->Divide(6,3);
 for (int j=0;j<3;j++){//direction
   c->cd(1+j*6+0); 
   hGuess[j]->Draw("colz");
   c->cd(1+j*6+1); 
   hDiff2D[j]->Draw("colz");
   c->cd(1+j*6+2); 
   hDiff[j]->Draw("hist");
   c->cd(1+j*6+3); 
   hBetterGuess[j]->Draw("colz");
   c->cd(1+j*6+4); 
   hBetterDiff2D[j]->Draw("colz");
   c->cd(1+j*6+5); 
   hBetterDiff[j]->Draw("hist");
 }
 return;
}
*/
