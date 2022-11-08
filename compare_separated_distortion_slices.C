
class Shifter {
public:
  Shifter(const char * filename);
  TVector3 Shift(TVector3 position);
  TVector3 ShiftForward(TVector3 position);
  TVector3 ShiftBack(TVector3 position);
  TVector3 GetDistortionAt(TVector3 position);
  TFile *file;
  TH3F *hX[2], *hY[2], *hZ[2];
};

Shifter::Shifter(const char * filename){
  file=0;
  file=TFile::Open(filename,"READ");
  assert (file && !file->IsZombie());
  hX[0]=(TH3F*)file->Get("hIntDistortionPosX");
  hY[0]=(TH3F*)file->Get("hIntDistortionPosY");
  hZ[0]=(TH3F*)file->Get("hIntDistortionPosZ");
  //need to add a flag to check naming conventions, probably.
  hX[1]=(TH3F*)file->Get("hIntDistortionNegX");
  hY[1]=(TH3F*)file->Get("hIntDistortionNegY");
  hZ[1]=(TH3F*)file->Get("hIntDistortionNegZ");
  return;
}

TVector3 Shifter::GetDistortionAt(TVector3 position){
  TVector3 distortion;
  double z= position.Z();
  double r=position.Perp();
  double phi=position.Phi();
  if(position.Phi() < 0.0){
    phi = position.Phi() + 2.0*TMath::Pi(); 
  }
  int side=(z<0);//sets side to '1'(negative) if z<0.
  //printf("seeking distortion at (rpz)=(%f,%f,%f)\n",r,phi,z);
  
  double xshift=hX[side]->Interpolate(phi,r,z);//coordinate of your stripe
  double yshift=hY[side]->Interpolate(phi,r,z);
  double zshift=hZ[side]->Interpolate(phi,r,z);

  distortion.SetXYZ(xshift,yshift,zshift);

  return distortion;
}

TVector3 Shifter::ShiftForward(TVector3 position){
  TVector3 forwardshift=position+GetDistortionAt(position);
  return forwardshift;
}

TVector3 Shifter::ShiftBack(TVector3 position){
  TVector3 backwardshift=position-GetDistortionAt(position);
  return backwardshift;
}
TVector3 Shifter::Shift(TVector3 position){
  
  return ShiftBack(ShiftForward(position));
}


void compare_separated_distortion_slices(const char* file1, const char* file2, float zin=1e8){
    //reads in a triplet of 3D integral distortion maps, and writes lots of plots, as well as some summary histograms to various output files.

  Shifter *shifter[2];
  shifter[0]=new Shifter(file1);
  shifter[1]=new Shifter(file2);

  //dimensions and bounds of the histogram
  int nph=60;
  int nrh=50;
  int nzh=50;
  float pih=0;
  float pfh=TMath::TwoPi();
  float rih=25;
  float rfh=75;
  float zih=-100;
  float zfh=100;

  float deltar=(rfh-rih)/nrh;
  float deltap=(pfh-pih)/nph;
  float deltaz=(zfh-zih)/(nzh);


    //monitor plots, and the position that plot monitors at:
  float zpos=zin;
  if (zpos>1e6) zpos=zfh/2.0;//(nzh/2+0.5)*deltaz+zih;
  TVector3 pos((nrh/2+0.5)*deltar+rih,0,zpos);
  pos.SetPhi((nph/2+0.5)*deltap+pih);
  float posphi=((nph/2+0.5)*deltap+pih);//make sure we don't have a wrapping issue.

  //and the bin in the output histogram where the monitor position occurs:
  int xi[3];
  //heads up that these cannot be converted to ints in a {} array, so must be done one at a time:
  xi[0]=(pos.Perp()-rih)/(deltar);
  xi[1]=(posphi-pih)/(deltap);
  xi[2]=(pos.Z()-zih)/(deltaz);

  // printf("slicing plots at pos=(%f,%f,%f)(rpz) (side=%d, histbins=(%d,%d,%d)\n",

  
	     //used to be{nrh/2,nph/2,nzh/2};
  const char axname[]="rpzrpz";
  int axn[]={nrh,nph,nzh,nrh,nph,nzh};
  float axval[]={(float)pos.Perp(),(float)pos.Phi(),(float)pos.Z(),(float)pos.Perp(),(float)pos.Phi(),(float)pos.Z()};
  float axbot[]={rih,pih,zih,rih,pih,zih};
  float axtop[]={rfh,pfh,zfh,rfh,pfh,zfh};
  TH1F* hRDist[3];//distortion as a function of R along the midline in z and phi, one for each component
  TH2F* hSummableDist2D[3][3];//distortion as a function of each midline, one for each component.  Stored as a 2D histogram so they can be added together.
  TH2F* hRDist2D[3];//legacy distortion as a function of R along the midline in z and phi, one for each component.  Stored as a 2D histogram so they can be added together.
  TH1F* hRModDist[3];
  float deltaR, deltaP;
  TH2F* hIntDist[3][3];
  TH1F* hRDiffDist[3];
  TH2F* hDiffDist[3][3];
  for (int i=0;i<3;i++){
    //loop over which axis of the distortion to read
    for (int ax=0;ax<3;ax++){
      //loop over which plane to work in
      hDiffDist[ax][i]=new TH2F(TString::Format("hDiffDist%c_%c%c",axname[i],axname[ax+1],axname[ax+2]),
				TString::Format("%c component of diff. distortion in  %c%c plane at %c=%2.3f;%c;%c",
						axname[i],axname[ax+1],axname[ax+2],axname[ax],axval[ax],axname[ax+1],axname[ax+2]),
				axn[ax+1],axbot[ax+1],axtop[ax+1],
				axn[ax+2],axbot[ax+2],axtop[ax+2]);
      hIntDist[ax][i]=new TH2F(TString::Format("hIntDist%c_%c%c",axname[i],axname[ax+1],axname[ax+2]),
			       TString::Format("%c component of int. distortion in  %c%c plane at %c=%2.3f;%c;%c",
					       axname[i],axname[ax+1],axname[ax+2],axname[ax],axval[ax],axname[ax+1],axname[ax+2]),
			       axn[ax+1],axbot[ax+1],axtop[ax+1],
			       axn[ax+2],axbot[ax+2],axtop[ax+2]);

      hSummableDist2D[ax][i]=new TH2F(TString::Format("hSummableDist2D%c_%c",axname[i],axname[ax]),
			       TString::Format("%c component of int. distortion vs %c, for all %c and %c;%c;#delta (um)",
					       axname[i],axname[ax],axname[ax+1],axname[ax+2],axname[ax]),
			       axn[ax],axbot[ax],axtop[ax],
				      100,-500,500);
    }
    hRDist[i]=new TH1F(TString::Format("hRDist%c",axname[i]),
		       TString::Format("%c component of int. distortion vs r with %c=%2.3f and %c=%2.3f;r(cm);#delta (cm)",
				       axname[i],axname[1],axval[1],axname[2],axval[2]),
		       axn[0],axbot[0],axtop[0]);
    hRDist2D[i]=new TH2F(TString::Format("hRDist2D%c",axname[i]),
		       TString::Format("%c component of int. distortion vs r at %c=%2.3f;r(cm);#delta (um)",
				       axname[i],axname[2],axval[2]),
		       axn[0],axbot[0],axtop[0],100,-500,500);
    hRModDist[i]=new TH1F(TString::Format("hRModDist%c",axname[i]),
		       TString::Format("#Delta %c (not %c-hat) from int. distortion vs r with %c=%2.3f and %c=%2.3f;r(cm);#delta (cm)",
				       axname[i],axname[i],axname[1],axval[1],axname[2],axval[2]),
		       axn[0],axbot[0],axtop[0]);
    hRDiffDist[i]=new TH1F(TString::Format("hRDiffDist%c",axname[i]),
		       TString::Format("%c component of diff. distortion vs r with %c=%2.3f and %c=%2.3f;r(cm);#delta (cm)",
				       axname[i],axname[1],axval[1],axname[2],axval[2]),
		       axn[0],axbot[0],axtop[0]);
  }

  

  //note that we apply the adjustment to the particle position (inpart) and not the plotted position (partR etc)
  int ir,ip,iz;
  unsigned long el=0;
  TVector3 distortVec[2],testPosition; //the vectors from each of the shifters.

  float distortR,distortP,distortZ;
  float diffDistortR,diffDistortP,diffDistortZ;
  float partR,partP,partZ;
  float distort[3];//convenient way of holding r,p,z components in loopable form
  float part[3];//convenient way of holding r,p,z coordinates in loopable form
  for (ir=0;ir<nrh;ir++){
    partR=(ir+0.5)*deltar+rih;
    testPosition.SetXYZ(partR,0,0);
    for (ip=0;ip<nph;ip++){
      partP=(ip+0.5)*deltap+pih;
      testPosition.SetPhi(partP);
      for (iz=0;iz<nzh;iz++){
	partZ=(iz)*deltaz+zih; //start us at the EDGE of the z bin, not the center?
	partZ+=0.5*deltaz; //move to center of histogram bin.
	testPosition.SetZ(partZ);
	distortVec[0]=shifter[0]->GetDistortionAt(testPosition);
	distortVec[1]=shifter[1]->GetDistortionAt(testPosition);

	distort[0]=distortR=distortVec[0].Perp()-distortVec[1].Perp();
	distort[1]=distortP=distortVec[0].DeltaPhi(distortVec[1])*testPosition.Perp();
	distort[2]=distortZ=distortVec[0].Z()-distortVec[1].Z();
	//fill the 2D summable plots, which integrate over the two axes that aren't involved:
	part[0]=partR;
	part[1]=partP;
	part[2]=partZ;
	for (int comp=0;comp<3;comp++){//loop over which component:
	  for (int ax=0;ax<3;ax++){//loop over which axis we're filling:
	    hSummableDist2D[ax][comp]->Fill(part[ax],distort[comp]*1e4);
	  }
	}

	
	//now we fill particular slices for integral visualizations:
	if(ir==xi[0]){//r slice
	  //printf("ir=%d, r=%f (pz)=(%d,%d), distortR=%2.2f, distortP=%2.2f\n",ir,partR,ip,iz,distortR,distortP);
	  hIntDist[0][0]->Fill(partP,partZ,distortR);
	  hIntDist[0][1]->Fill(partP,partZ,distortP);
	  hIntDist[0][2]->Fill(partP,partZ,distortZ);
	}
	if(ip==xi[1]){//phi slice
	  //printf("ip=%d, p=%f (rz)=(%d,%d), distortR=%2.2f, distortP=%2.2f\n",ip,partP,ir,iz,distortR,distortP);
	  hIntDist[1][0]->Fill(partZ,partR,distortR);
	  hIntDist[1][1]->Fill(partZ,partR,distortP);
	  hIntDist[1][2]->Fill(partZ,partR,distortZ);
	}
	if(iz==xi[2]){//z slice
	  //printf("iz=%d, z=%f (rp)=(%d,%d), distortR=%2.2f, distortP=%2.2f\n",iz,partZ,ir,ip,distortR,distortP);
		  
	  hIntDist[2][0]->Fill(partR,partP,distortR);
	  hIntDist[2][1]->Fill(partR,partP,distortP);
	  hIntDist[2][2]->Fill(partR,partP,distortZ);
	  hRDist2D[0]->Fill(partR,distortR*1e4,1);//convert cm to um
	  hRDist2D[1]->Fill(partR,distortP*1e4,1);//convert cm to um
	  hRDist2D[2]->Fill(partR,distortZ*1e4,1);//convert cm to um


	  if(ip==xi[1]){//phi slices of z slices = r line at mid phi, mid z:
	    hRDist[0]->Fill(partR,distortR);	    
	    hRDist[1]->Fill(partR,distortP);
	    hRDist[2]->Fill(partR,distortZ);
	  }

	}

	
      }
    }
  }


  TCanvas *canvas=new TCanvas("cdistort","distortion integrals",1200,800);
  //take 10 of the bottom of this for data?
  canvas->cd();
  TPad *c=new TPad("cplots","distortion integral plots",0,0.2,1,1);
  canvas->cd();
  TPad *textpad=new TPad("ctext","distortion integral plots",0,0.0,1,0.2);
  c->Divide(4,3);
  gStyle->SetOptStat();
  for (int i=0;i<3;i++){
    //component
    for (int ax=0;ax<3;ax++){
      //plane
      c->cd(i*4+ax+1);
      gPad->SetRightMargin(0.15);
      hIntDist[ax][i]->SetStats(0);
      hIntDist[ax][i]->Draw("colz");
    }
    c->cd(i*4+4);
    hRDist[i]->SetStats(0);
    hRDist[i]->SetFillColor(kRed);
    hRDist[i]->Draw("hist");
  }
  textpad->cd();
  float texpos=0.9;float texshift=0.12;
  TLatex *tex=new TLatex(0.0,texpos,"Fill Me In");
  tex->SetTextSize(texshift*0.8);
  tex->DrawLatex(0.05,texpos,"slices of DIFFERNCE between two integral distortions");texpos-=texshift;
 //tex->DrawLatex(0.05,texpos,Form("Drift Field = %2.2f V/cm",GetNominalE()));texpos-=texshift;
  tex->DrawLatex(0.05,texpos,Form("Sampling distortion in (rpz)=(%d x %d x %d) grid.", nrh,nph,nzh));texpos-=texshift;
  tex->DrawLatex(0.05,texpos,Form("'A' file: %s",file1));texpos-=texshift;
  tex->DrawLatex(0.05,texpos,Form("'B' file: %s",file2));texpos-=texshift;
  texpos=0.9;
 
  
  
  canvas->cd();
  c->Draw();
  canvas->cd();
  textpad->Draw();
  canvas->SaveAs("LowVsHighRes.distortion_summary.pdf");
  //canvas->SaveAs(gifFilename.Data());



  
  return;
  /*

  outputHistFile->cd();
  for (int i=0;i<3;i++){
    for (int j=0;j<3;j++){
      hIntDist[i][j]->Write();
      hDiffDist[i][j]->Write();
    }
    hRDist[i]->Write();
    hRDist2D[i]->Write();
  }

  for (int comp=0;comp<3;comp++){//loop over which component:
    for (int ax=0;ax<3;ax++){//loop over which axis we're filling:
      hSummableDist2D[ax][comp]->Write();
    }
  }
  
  outputHistFile->Close();
  printf("wrote hists to:%s.\n",histFilename.Data());

  return;
  */
}
