#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"

float PI=3.1415927;
void basic_interpolate(){
  const float rmin=24;
  const float rmax=78;
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


  float IFCrange[]={0.,15.};
  float OFCrange[]={25.,75.};
  float boundM[2];
  float boundB[2];
  for (int i=0;i<2;i++){
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
      float z=zstep*(iz+0.5)+zmin;
      if ((z<boundM[0]*r+boundB[0])||(z>boundM[1]*r+boundB[1])){
	continue;
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
  for (int j=0;j<3;j++){
    hGuess[j]=new TH2F(Form("hGuess%d",j),Form("interpolated %c-distortion in r-z plane at phi=%2.2f;z (cm);r (cm)",dirchar[j],phiplane[1]),nzsteps,zmin,zmax,nrsteps,rmin,rmax);
    hGuess[j]->Add(hSlice[0][j],hSlice[2][j],0.5,0.5);
    hDiff2D[j]=new TH2F(Form("hDiff2D%d",j),Form("%c-distortion Difference Between INterpolation and true (um) ;z (cm);r (cm)",dirchar[j]),nzsteps,zmin,zmax,nrsteps,rmin,rmax);
    hDiff2D[j]->Add(hSlice[1][j],hGuess[j],-1.0,1.0);
    hDiff[j]=new TH1F(Form("hDiff%d",j),Form("%c-distortion Difference Between INterpolation and true;diff (um)",dirchar[j]),200,-400,400);
  }

  for (int ir=0;ir<nrsteps;ir++){
    float r=rstep*(ir+0.5)+rmin;
    for (int iz=0;iz<nzsteps;iz++){
      float z=zstep*(iz+0.5)+zmin;
      if ((z<boundM[0]*r+boundB[0])||(z>boundM[1]*r+boundB[1])){
	continue;
      }
      int bin=hGuess[0]->FindBin(z,r);
      for (int j=0;j<3;j++){
	hDiff[j]->Fill(1e4*(hGuess[j]->GetBinContent(bin)-hSlice[1][j]->GetBinContent(bin)));
      }
    }
  } //read the true x,y,z distortions at phi[0] and phi[1]
  //create a linear function for all intermediate phi
  //compare the values in these intermediate regions to the true hist values there
 TCanvas *c=new TCanvas("c","c",1900,800);
 c->Divide(6,3);
 for (int j=0;j<3;j++){//direction
   c->cd(1+j*6); 
   hSlice[0][j]->Draw("colz");
   c->cd(1+j*6+1); 
   hSlice[2][j]->Draw("colz");
   c->cd(1+j*6+2); 
   hSlice[1][j]->Draw("colz");
   c->cd(1+j*6+3); 
   hGuess[j]->Draw("colz");
   c->cd(1+j*6+4); 
   hDiff2D[j]->Draw("colz");
   c->cd(1+j*6+5); 
   hDiff[j]->Draw("colz");
 }



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
