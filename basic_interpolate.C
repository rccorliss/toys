#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"

const float PI=3.1415927;

const int ArbitraryCMzBin=1; //which bin we read the 'CM' values from. Should be replaced with a z position eventually.


//z-coverage, in cm, of the outer tracker
//we assume that the full-length OT is centered at phi=0.
const  float IFCrange[]={0.,15.,15.};
const  float OFCrange[]={25.,75.,105.};
const float OFCradius=78.;
const float OTwidth=25.;
const int OTnDetectors=12;
float OTspacing=2*PI/(1.0*OTnDetectors); //phi distance between detector centersoverridden by 2PI/OTnDetectors;
float OTphiWidth=OTwidth/OFCradius;//width in radians
float OTphiHalfWidth=OTphiWidth/2.;//half-width in radians (=distance from center to edge)
float OTseparation=OTspacing-OTphiWidth;//phi distance between end of one and start of the next.

float boundM[]={0.,0.,0.};
float boundB[]={0.,0.,0.};

TH3F* GetDistortionModel(int nr_cm,int np_cm, int nz_cm, int nr_ot,int np_ot, int nz_ot, TH3F *truth);
void PlotClassicInterpolation();
bool CoveredBySingleOT(float r, float phi, float z){
  if (phi<0)phi+=2*PI;
  float distlower=fmod(phi,OTspacing); //distance from first center it is larger than.

  //check to see if we are past the upper edge of the lower tile, and short of the lower edge of the upper:
  if (distlower>OTphiHalfWidth && distlower<OTspacing-OTphiHalfWidth) return false;

  //otherwise check to see if we're in the z coverage range:
  float minz=(boundM[0]*r+boundB[0]);
  float maxz=(boundM[1]*r+boundB[1]);
  //printf("checking for %2.2f< z(%2.2f)< %2.2f\n",minz,z,maxz);
  return ((z>minz) && (z<maxz));
}
bool CoveredByDoubleOT(float r, float phi, float z){
  if (phi<0) phi+=2*PI;
  //static const float tile_phi_margin=(OTwidth*0.5)/OFCradius; //radians from center to edge of tile
  if (phi>OTphiHalfWidth && phi<2*PI-OTphiHalfWidth) return false;
  return true;
}

int GetPhiRegion(float phi){
  if (phi<0) phi+=2*PI;
  return (int)(phi/OTspacing);
}
float GetPhiInterpolationFractionPhi(float phi){
  float distcenter=fmod(phi,OTspacing); //distance from first center it is larger than.
  float distedge=distcenter-OTphiHalfWidth;//distance from upper edge of that tile;
  float fractional=distedge/OTseparation;
  if (fractional<0.) return 0; //not past the upper edge of the lower tile;
  if (fractional>1.) return 1; //past the lower edge of the upper tile;

  return fractional;
}


void basic_interpolate(){

  // PlotClassicInterpolation(); return;

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

  float phimin=0.0;
  float phimax=6.28;
  float phispan=phimax-phimin;
  int nphisteps=60;
  float phistep=phispan/(1.0*nphisteps);


  //float boundM[3];
  //float boundB[3];
  for (int i=0;i<3;i++){
    boundM[i]=(OFCrange[i]-IFCrange[i])/(rmax-rmin);
    boundB[i]=IFCrange[i]-rmin*boundM[i];
    printf("bounds[%d]: m=%2.2f, b=%2.2f\n",i,boundM[i],boundB[i]);
  }
  //return;

  
  
  //load the time-average histogram(s):
  TFile *histfile;
  //histfile=TFile::Open("elevatorpitch/fluct_average.rev3.1side.3d.file0.h_negz.real_B1.4_E-400.0.ross_phi1_sphenix_phislice_lookup_r26xp40xz40.distortion_map.hist.root","READ");
  histfile=TFile::Open("elevatorpitch/average.rev3.1side.3d.file0.h_negz.real_B1.4_E-400.0.ross_phi1_sphenix_phislice_lookup_r26xp40xz40.distortion_map.hist.root","READ");
  TH3F *hData[3];
  hData[0]=(TH3F*)histfile->Get("hIntDistortionR");
  hData[1]=(TH3F*)histfile->Get("hIntDistortionP");
  hData[2]=(TH3F*)histfile->Get("hIntDistortionZ");

  printf("data loaded.\n");


  TH3F* hModel[3];
  TH3F * hComparison[3];
  TH1F * hComparisonDiff[3];
  int i=0;
  hData[i]->Scale(1e4);//shift from cm to um for all that follows.
  hModel[i]=GetDistortionModel(20,12*3, 1, 10,10, 10, hData[i]);
  hComparison[i]=(TH3F*)(hData[i]->Clone(Form("hComparison%d",i)));
  hComparison[i]->Reset();
  return;
  hComparisonDiff[i]=new TH1F(Form("hComparisonDiff%d",i),"Difference between model and truth, r>30;diff (um)",200,-200,200);

  //at every point in the true model plot the difference between truth and the model we've built:

  int nr_tr=hData[0]->GetYaxis()->GetNbins();
  float rmin_tr=hData[0]->GetYaxis()->GetXmin();
  float rmax_tr=hData[0]->GetYaxis()->GetXmax();
  float rstep_tr=(rmax_tr-rmin_tr)/(1.0*nr_tr);

  int np_tr=hData[0]->GetXaxis()->GetNbins();
  float pmin_tr=hData[0]->GetXaxis()->GetXmin();
  float pmax_tr=hData[0]->GetXaxis()->GetXmax();
  float pstep_tr=(pmax_tr-pmin_tr)/(1.0*np_tr);

  int nz_tr=hData[0]->GetZaxis()->GetNbins();
  float zmin_tr=hData[0]->GetZaxis()->GetXmin();
  float zmax_tr=hData[0]->GetZaxis()->GetXmax();
  float zstep_tr=(zmax_tr-zmin_tr)/(1.0*nz_tr);

  for (int ir=0;ir<nr_tr;ir++){
    float r=rmin_tr+(ir+0.5)*rstep_tr;// position of _center_ of the ir-th bin of the source, for histogram filling safety.
    if (r<activermin) continue;
    if (r>rmax) continue;
    for (int ip=0;ip<np_tr;ip++){
      float p=pmin_tr+(ip+0.5)*pstep_tr;// position of _center_ of the ip-th bin of the source, for histogram filling safety.
      for (int iz=0;iz<nz_tr;iz++){
	float z=zmin_tr+(iz+0.5)*zstep_tr;// position of _center_ of the iz-th bin of the source, for histogram filling safety.
	if (z<zmin || z>zmax) continue;
	float valData=hData[i]->GetBinContent(hData[i]->FindBin(p,r,z));
	float valModel=hModel[i]->GetBinContent(hModel[i]->FindBin(p,r,z));
	printf("prz=(%2.2f,%2.2f,%2.2f).  data=%2.2f\tmodel=%2.2f\tdiff=%2.2f\n",p,r,z,valData,valModel,valModel-valData);
	hComparison[i]->Fill(p,r,z,valModel-valData);
	hComparisonDiff[i]->Fill((valModel-valData));
      }
    }
  }
  
  

  
 TCanvas *c=new TCanvas("c","c",1200,400);
 c->Divide(4,1);
 for (int j=0;j<1;j++){//direction
   c->cd(1+j*6+0);
   //hData[i]->GetXaxis()->SetRangeUser(1,1);//slice in phi
   TH2F* hDataProj=(TH2F*)(hData[i]->Project3D("yz"));
   hDataProj->SetTitle("True Distortion in phibin=1 (um);r (cm);z (cm)");
   hDataProj->Draw("colz");//plot in z
   hData[i]->GetXaxis()->SetRange(0,0);//reset the slice
   c->cd(1+j*6+1); 
   //hModel[i]->GetXaxis()->SetRange(1,1);//slice in phi
   TH2F* hModelProj=(TH2F*)(hModel[i]->Project3D("yz"));
   hModelProj->SetTitle("Reco Distortion in phibin=1 (um);r (cm);z (cm)");
   hModelProj->Draw("colz");//plot in z   hModel[i]->Project3D("yz")->Draw("colz");//plot in z
   hModel[i]->GetXaxis()->SetRange(0,0);//reset the slice
    c->cd(1+j*6+2); 
    //hComparison[i]->GetXaxis()->SetRange(1,1);//slice in phi
    hComparison[i]->Project3D("yz")->Draw("colz");//plot in z
   //hComparison[i]->GetXaxis()->SetRange(0,0);//reset the slice
   c->cd(1+j*6+3); 
   hComparisonDiff[i]->Draw();
   c->cd(1+j*6+4); 
   c->cd(1+j*6+5); 
 }

return;
}


void PlotClassicInterpolation(){
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
  float deltaPhi=PI;//2*PI/12.0;
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
    hDiff[j]=new TH1F(Form("hDiff%d",j),Form("%c-distortion Difference Between INterpolation and true (r>%1.1f);diff (um)",dirchar[j],activermin),200,-200,200);
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
    hBetterDiff[j]=new TH1F(Form("hBetterDiff%d",j),Form("%c-distortion Difference Between int++ and true (r>%1.1f);diff (um)",dirchar[j],activermin),200,-200,200);


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
      hForwardGuess[i][j]=new TH2F(Form("hForwardGuess%d%d",i,j),Form("#delta %c-hat in r-z plane at #phi=%2.2f using %s;z (cm);r (cm)",dirchar[j],phiplane[2],guessName[i].Data()),nzsteps,zmin,zmax,nrsteps,rmin,rmax);
      hForwardDiff[i][j]=new TH1F(Form("hForwardDiff%d%d",i,j),Form("(Model%d-True) for %c-hat (r>%1.1f);diff (um)",i,dirchar[j],activermin),200,-150,150);
      hForwardDiff2D[i][j]=new TH2F(Form("hForwardDiff2D%d%d",i,j),Form("(Model%d-True) for %c-hat (um);z (cm);r (cm)",i,dirchar[j]),nzsteps,zmin,zmax,nrsteps,rmin,rmax);
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
	  //printf("cmrescale%d=%f=(%f/%f)\n",j,cmrescale[j],cmdestination,cmsource);
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

 c=new TCanvas("cforwardbest","cforwardbest",950,800);
 c->Divide(3,3);
 for (int j=0;j<3;j++){//direction
   c->cd(1+j*3);
   hForwardGuess[2][j]->Draw("colz");
   c->cd(2+j*3);
   hForwardDiff2D[2][j]->Draw("colz");
   c->cd(3+j*3)->SetLogy(); 
   hForwardDiff[2][j]->Draw("hist");
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


TH3F* GetDistortionModel(int nr_cm,int np_cm, int nz_cm, int nr_ot,int np_ot, int nz_ot, TH3F *truth){
  static int nmod=0;

  const float rmin=24.;
  const float rmax=78.;
  float rspan=rmax-rmin;
  int nrsteps=nr_ot;
  float rstep=rspan/(1.0*nrsteps);
  const float zmin=0;
  const float zmax=105.5;
  float zspan=zmax-zmin;
  int nzsteps=nz_ot;
  float zstep=zspan/(1.0*nzsteps);
  float pmin=0.0;
  float pmax=6.28;
  float pspan=pmax-pmin;
  int npsteps=np_ot;
  float pstep=pspan/(1.0*npsteps);


  
  //define our output histogram, with extra edge steps like for the real thing.
  TH3F *hDistModel=new TH3F(Form("hDistModel%d",nmod),Form("Distortion Model %d",nmod),
		       npsteps+2,pmin-pstep,pmax+pstep,
		       nrsteps+2,rmin-rstep,rmax+rstep,
		       nzsteps+2,zmin-zstep,zmax+zstep);
  TH3F *hDistModelSamples=new TH3F(Form("hDistModelSamples%d",nmod),Form("Distortion Model Samples%d",nmod),
		       npsteps+2,pmin-pstep,pmax+pstep,
		       nrsteps+2,rmin-rstep,rmax+rstep,
		       nzsteps+2,zmin-zstep,zmax+zstep);

  //define our helper slices:
  //the low and high edge of every partial region:
  TH2F *hLowEdgeSlice[OTnDetectors];
  TH1F *hLowEdgeCM[OTnDetectors];
  TH1F *hLowEdgeCMsamples[OTnDetectors];
  TH2F *hLowEdgeSamples[OTnDetectors]; //for normalization
  TH2F *hHighEdgeSlice[OTnDetectors];
  TH2F *hHighEdgeSamples[OTnDetectors]; //for normalization
  TH1F *hHighEdgeCM[OTnDetectors];
  TH1F *hHighEdgeCMsamples[OTnDetectors];

  for (int i=0;i<OTnDetectors;i++){
    hLowEdgeSlice[i]=new TH2F(Form("hLowEdgeSlice%d",i),Form("lowedge slice in r-z plane, sector %d;z (cm);r (cm)",i),nzsteps,zmin,zmax,nrsteps,rmin,rmax);
    hLowEdgeSamples[i]=new TH2F(Form("hLowEdgeSamples%d",i),Form("lowedge samples in r-z plane, sector %d;z (cm);r (cm)",i),nzsteps,zmin,zmax,nrsteps,rmin,rmax);
    hLowEdgeCM[i]=new TH1F(Form("hLowEdgeCM%d",i),Form("CM values at lowedge slice sector %d;r (cm);CMval",i),nrsteps,rmin,rmax);
    hLowEdgeCMsamples[i]=new TH1F(Form("hLowEdgeCMsamples%d",i),Form("CM samples at lowedge slice sector %d;r (cm);CMval",i),nrsteps,rmin,rmax);
    hHighEdgeSlice[i]=new TH2F(Form("hHighEdgeSlice%d",i),Form("highedge slice in r-z plane, sector %d;z (cm);r (cm)",i),nzsteps,zmin,zmax,nrsteps,rmin,rmax);
    hHighEdgeSamples[i]=new TH2F(Form("hHighEdgeSamples%d",i),Form("highedge samples in r-z plane, sector %d;z (cm);r (cm)",i),nzsteps,zmin,zmax,nrsteps,rmin,rmax);
    hHighEdgeCM[i]=new TH1F(Form("hHighEdgeCM%d",i),Form("CM values at highedge slice sector %d;r (cm);CMval",i),nrsteps,rmin,rmax);
    hHighEdgeCMsamples[i]=new TH1F(Form("hHighEdgeCMsamples%d",i),Form("CM samples at highedge slice sector %d;r (cm);CMval",i),nrsteps,rmin,rmax);

    
  }
  
  //and the full region, which we'll average over all such slices:
  TH2F *hFullSlice=new TH2F("hFullSlice","full OT slice in r-z plane, sector 0;z (cm);r (cm)",nzsteps,zmin,zmax,nrsteps,rmin,rmax);
  TH2F *hFullSamples=new TH2F("hFullSamples","full OT samples in r-z plane, sector 0;z (cm);r (cm)",nzsteps,zmin,zmax,nrsteps,rmin,rmax); //for normalization

  //and finally the central membrane:
  TH2F *hCM=new TH2F("hCM","Distortion at CM;phi (rad);r (cm)",np_cm,pmin,pmax,nr_cm,rmin,rmax);
  TH2F *hCMSamples=new TH2F("hCMSamples","Samples at CM;phi (rad);r (cm)",np_cm,pmin,pmax,nr_cm,rmin,rmax);


  //now we loop over all elements in the source histogram and fill the appropriate helper slices:
  
  int nr_tr=truth->GetYaxis()->GetNbins();
  float rmin_tr=truth->GetYaxis()->GetXmin();
  float rmax_tr=truth->GetYaxis()->GetXmax();
  float rstep_tr=(rmax_tr-rmin_tr)/(1.0*nr_tr);

  int np_tr=truth->GetXaxis()->GetNbins();
  float pmin_tr=truth->GetXaxis()->GetXmin();
  float pmax_tr=truth->GetXaxis()->GetXmax();
  float pstep_tr=(pmax_tr-pmin_tr)/(1.0*np_tr);

  int nz_tr=truth->GetZaxis()->GetNbins();
  float zmin_tr=truth->GetZaxis()->GetXmin();
  float zmax_tr=truth->GetZaxis()->GetXmax();
  float zstep_tr=(zmax_tr-zmin_tr)/(1.0*nz_tr);

  for (int ir=0;ir<nr_tr;ir++){
    float r=rmin_tr+(ir+0.5)*rstep_tr;// position of _center_ of the ir-th bin of the source, for histogram filling safety.
    if (r<rmin || r>rmax) continue;
    for (int ip=0;ip<np_tr;ip++){
      float p=pmin_tr+(ip+0.5)*pstep_tr;// position of _center_ of the ip-th bin of the source, for histogram filling safety.
    //int psector=(p+PI/12.)/(PI/6);//position relative to lower edge of zero, in units of center-to-center spacing
    int psector=GetPhiRegion(p+OTphiHalfWidth);//which OT sector we _might_ be in 
    if (psector>=OTnDetectors) psector=0;
      if (p<pmin || p>pmax) continue;
      for (int iz=0;iz<nz_tr;iz++){
	float z=zmin_tr+(iz+0.5)*zstep_tr;// position of _center_ of the iz-th bin of the source, for histogram filling safety.
	if (z<zmin || z>zmax) continue;
	//printf("in detector bulk region\n");

	float val=truth->Interpolate(p,r,z);
	
	if (CoveredByDoubleOT(r,p,z)){ //fill the full-length sample set
	  //printf("in DoubleOT det\n");
	  hFullSlice->Fill(z,r,val);
	  hFullSamples->Fill(z,r,1);
	}
	
	if (iz==ArbitraryCMzBin){ //fill the CM
	  printf("at CM\n");
	  hCM->Fill(p,r,val);
	  hCMSamples->Fill(p,r,1);
	}

	if (CoveredBySingleOT(r,p,z)){ //fill the low or high edges of individual tiles.
	    printf("in singleOT det\n");
	  if (!CoveredBySingleOT(r,p+pstep_tr,z)){//next step is not covered, so we're on a high edge:
	    hHighEdgeSlice[psector]->Fill(z,r,val);
	    hHighEdgeSamples[psector]->Fill(z,r,1.);
	    printf("high edge fill\n");
	    if (iz==ArbitraryCMzBin){
	      hHighEdgeCM[psector]->Fill(r,val);
	      hHighEdgeCMsamples[psector]->Fill(r,1);
	    }
	  }
	  if (!CoveredBySingleOT(r,p-pstep_tr,z)){//previous step is not covered, so we're on a low edge:
	    hLowEdgeSlice[psector]->Fill(z,r,val);
	    hLowEdgeSamples[psector]->Fill(z,r,1.);
	    if (iz==ArbitraryCMzBin){
	      hLowEdgeCM[psector]->Fill(r,val);
	      hLowEdgeCMsamples[psector]->Fill(r,1);
	    }
	  }	    
	}

	if (CoveredBySingleOT(r,p,z) || CoveredByDoubleOT(r,p,z)){//fill directly-measured distortion:
	  hDistModel->Fill(p,r,z,val);
	  hDistModelSamples->Fill(p,r,z,1);
	}
      }
    }
  }

  //normalize our various slices for the number of entries in each:
  //note that Divide() sets bins to zero if the denominator bin is zero.
  hFullSlice->Divide(hFullSamples);
  hCM->Divide(hCMSamples);
  for (int i=0;i<12;i++){
    hLowEdgeSlice[i]->Divide(hLowEdgeSamples[i]);
    hHighEdgeSlice[i]->Divide(hHighEdgeSamples[i]);
    hLowEdgeCM[i]->Divide(hLowEdgeCMsamples[i]);
    hHighEdgeCM[i]->Divide(hHighEdgeCMsamples[i]);
  }
  hDistModel->Divide(hDistModelSamples);


  
  //now that we have defined our partial OT slices, we extend those using the full OT slice
  for (int i=0;i<OTnDetectors;i++){
    float plow=(i*OTphiWidth)-OTphiHalfWidth; //phi position of the low edge slice
    float phigh=(i*OTphiWidth)+OTphiHalfWidth; //phi position of the high edge slice
    if (plow<0) plow+=2*PI;

    float midZrescale[2];
    for (int ir=0;ir<nrsteps;ir++){
      float r=rstep*(ir+0.5)+rmin;
       for (int iz=0;iz<nzsteps;iz++){
	float z=zstep*(iz+0.501)+zmin;
	//calculate the midZ rescale on the last bin before we're in the forward region
	if (CoveredBySingleOT(r,0,z) && !CoveredBySingleOT(r,0,z+zstep)){
	  int midzbin=hFullSlice->FindBin(z,r);
	  float midzsource=hFullSlice->GetBinContent(midzbin);
	  midZrescale[0]=hLowEdgeSlice[i]->GetBinContent(midzbin)/midzsource;
	  midZrescale[1]=hHighEdgeSlice[i]->GetBinContent(midzbin)/midzsource;
	}    
  	if (!CoveredByDoubleOT(r,0,z)){//nothing to extrapolate with)
	  continue;
	}
	//printf("rz=%d,%d\n",ir,iz);

	int bin=hFullSlice->FindBin(z,r);
	float dist=hFullSlice->GetBinContent(bin);
	hLowEdgeSlice[i]->Fill(z,r,dist*midZrescale[0]);
	hHighEdgeSlice[i]->Fill(z,r,dist*midZrescale[1]);
      }
    }
  }

  //extrapolate all cells covered by single OT using the same scheme:
  for (int ip=0;ip<npsteps;ip++){
    float p=pstep*(ip+0.5)+pmin;
    
    if (!CoveredBySingleOT(OFCradius,p,OFCrange[1]-OFCrange[0])) continue;//awkward way to check whether the phi range is okay by picking a known-good r-z combination
    float midZrescale;
    for (int ir=0;ir<nrsteps;ir++){
      float r=rstep*(ir+0.5)+rmin;
      for (int iz=0;iz<nzsteps;iz++){
	float z=zstep*(iz+0.501)+zmin;
	//calculate the midZ rescale on the last bin before we're in the forward region
	if (CoveredBySingleOT(r,p,z) && !CoveredBySingleOT(r,p,z+zstep)){
	  int midzbinSlice=hFullSlice->FindBin(z,r);
	  float midzsource=hFullSlice->GetBinContent(midzbinSlice);
	  int midzbin=hDistModel->FindBin(p,r,z);
	  midZrescale=hDistModel->GetBinContent(midzbin)/midzsource;
	}    
	if (!CoveredByDoubleOT(r,0,z)){//nothing to extrapolate with)
	  continue;
	}
	//printf("rz=%d,%d\n",ir,iz);
	if (CoveredByDoubleOT(r,0,z) && !CoveredBySingleOT(r,p,z)){//fill in the bins for regions covered by the full but not the partial
	  int binSlice=hFullSlice->FindBin(z,r);
	  float dist=hFullSlice->GetBinContent(binSlice);
	  hDistModel->Fill(p,r,z,dist*midZrescale);
	}
      }
    }
  }


  //now all regions covered by a single OT have been extended to the full OT range.
  
  //interpolate all !OT cells we can interpolate from the existing extended slices:
  for (int ip=0;ip<npsteps;ip++){
    float p=pstep*(ip+0.5)+pmin;
    if (p>2*PI) p-=2*PI;

    //skip regions in phi that correspond to detectors
    if (CoveredBySingleOT(OFCradius,p,OFCrange[1]-OFCrange[0])) continue;//awkward way to check whether the phi range is covered by picking a known-good r-z combination
    int psector=GetPhiRegion(p);//gives us which sector we're above, hence which low and high slice to use.
    int phigh=(psector+1)%OTnDetectors;
    float pfraction=GetPhiInterpolationFractionPhi(p);//gives us the fraction of the distance from the high edge we're passed to the low edge we haven't reached yet.
    for (int ir=0;ir<nrsteps;ir++){
      float r=rstep*(ir+0.5)+rmin;
      
      //calculate the CM rescale at this phi-r position:
      int cmBinEdge=hLowEdgeCM[0]->FindBin(r);
      float cmsource=hLowEdgeCM[phigh]->GetBinContent(cmBinEdge)*pfraction+hHighEdgeCM[psector]->GetBinContent(cmBinEdge)*(1-pfraction);
      int cmBin=hCM->FindBin(p,r);
      float cmRescale=hCM->GetBinContent(cmBin)/cmsource;
            
      for (int iz=0;iz<nzsteps;iz++){
	float z=zstep*(iz+0.501)+zmin;
	//printf("rz=%d,%d\n",ir,iz);
	if (CoveredByDoubleOT(r,0,z)){//fill in the bins for regions covered by the full but not the partial
	  int binSlice=hLowEdgeSlice[0]->FindBin(z,r);
	  float dist=hLowEdgeSlice[phigh]->GetBinContent(binSlice)*pfraction+hHighEdgeSlice[psector]->GetBinContent(binSlice)*(1-pfraction);;
	  hDistModel->Fill(p,r,z,dist*cmRescale);
	}
      }
    }
  }




  
  //now all regions have been extended to cover the double-OT range.  The only thing left is to interpolate between those and the known boundaries
  //interpolate all !OT cells we can interpolate from the existing extended slices:
  for (int ir=0;ir<nrsteps;ir++){
    float r=rstep*(ir+0.5)+rmin;

    int zBinLow=0;//z index of the first bin in dOT coverage at this r
    float zLow=0;//z center of that bin
    for (int iz=0;iz<nzsteps;iz++){
      float z=zstep*(iz+0.501)+zmin;
      zBinLow=iz;
      zLow=z;
      if (CoveredByDoubleOT(r,0,z)) {//stop after updating the values, if this is in now in the region.
	break;
      }
    }
    
    int zBinHigh=0;//z index of the last bin in dOT coverage at this r
    float zHigh=0;//z center of that bin
    for (int iz=zBinLow;iz<nzsteps;iz++){
      float z=zstep*(iz+0.501)+zmin;
      if (!CoveredByDoubleOT(r,0,z)) {
	break;
      }//stop before updating the values, if we're no longer covered.
      zBinHigh=iz;
      zHigh=z;
    }
    
    for (int ip=0;ip<npsteps;ip++){
      float p=pstep*(ip+0.5)+pmin;
      if (p>2*PI) p-=2*PI;

      float CMval=hCM->GetBinContent(hCM->FindBin(p,r));
      float lowVal=hDistModel->GetBinContent(hDistModel->FindBin(p,r,zLow));
      float highVal=hDistModel->GetBinContent(hDistModel->FindBin(p,r,zHigh));
      float endVal=0;//no distortion at the endcap itself.

      for (int iz=0;iz<nzsteps;iz++){
	float z=zstep*(iz+0.501)+zmin;
	//printf("rz=%d,%d\n",ir,iz);
	if (iz>=zBinLow || iz<=zBinHigh) continue; //don't need the bins we've already filled.

	float dist=0;//interpolated distortion, no further rescaling since we're only interpolating in z.
	if (iz<zBinLow){
	  float zfraction=(1.0*iz)/(1.0*zBinLow);
	  dist=lowVal*zfraction+CMval*(1-zfraction);
	}
	if (iz>zBinHigh){
	  float zfraction=(1.0*(iz-zBinHigh))/(1.0*(nzsteps-zBinHigh));
	  dist=endVal*zfraction+highVal*(1-zfraction);
	}
	hDistModel->Fill(p,r,z,dist);
      }
    }
  }




  
  //as a last step, we fill in the edge bins:
 
  for (int ir=0;ir<nrsteps;ir++){
    float r=rstep*(ir+0.5)+rmin;
    for (int iz=0;iz<nzsteps;iz++){
      float z=zstep*(iz+0.501)+zmin;
      for (int ip=0;ip<npsteps;ip++){
	float p=pstep*(ip+0.5)+pmin;
	//copy r values at neighbor r into the edge bins:
	if (ir==0) hDistModel->Fill(p,r-rstep,z,hDistModel->GetBinContent(hDistModel->FindBin(p,r,z)));
	if (ir==nrsteps-1) hDistModel->Fill(p,r+rstep,z,hDistModel->GetBinContent(hDistModel->FindBin(p,r,z)));
	//wrap around phi values into the edge bins:
	if (ip==0) hDistModel->Fill(p-pstep,r,z,hDistModel->GetBinContent(hDistModel->FindBin(p-pstep+2.*PI,r,z)));
	if (ip==npsteps-1) hDistModel->Fill(p+pstep,r,z,hDistModel->GetBinContent(hDistModel->FindBin(p+pstep-2.*PI,r,z)));
	//copy z values at neighbor z into the edge bins.
	if (iz==0) hDistModel->Fill(p,r,z-zstep,hDistModel->GetBinContent(hDistModel->FindBin(p,r,z)));
	if (iz==nzsteps-1) hDistModel->Fill(p,r,z+zstep,hDistModel->GetBinContent(hDistModel->FindBin(p,r,z)));
      }
    }
  }

  TCanvas *c=new TCanvas("cModelParts","cModelParts",1000,1000);
  c->Divide(6,4);
  for (int i=0;i<OTnDetectors;i++){
    c->cd(2*i+1);
    hLowEdgeSamples[i]->Draw("colz");
   c->cd(2*i+2);
    hHighEdgeSamples[i]->Draw("colz");
  }
  return hDistModel;
}
