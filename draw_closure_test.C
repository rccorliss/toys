

void draw_closure_test(const char* filebase, float zin=1e8){
  //reads in a triplet of 3D integral distortion maps, and writes lots of plots, as well as some summary histograms to various output files.

  const float MAX_PLOT_DISTORTION=0.02;//2 cm.  huge!  0.0200;//200um.

  //distortion:
  TString distortionFilename;
  distortionFilename=filebase;
  TString basename=filebase;
  basename.ReplaceAll(".distortion_map.hist.root","");

  //outputs:
  TString summaryFilename;
  summaryFilename.Form("%s.z1.0.separate_closure_summary.pdf",basename.Data());
  TString gifFilename=Form("%s.z1.0.separate_closure_summary.gif",basename.Data());
 TString diffFilename;
  diffFilename.Form("%s.z1.0.separate_closure2_summary.pdf",basename.Data());
  TString diffgifFilename=Form("%s.z1.0.separate_closure2_summary.gif",basename.Data());
  TString histFilename=Form("%s.z1.0.separate_closure2_summary.hist.root",basename.Data());
  TFile *outputHistFile=TFile::Open(histFilename.Data(),"RECREATE");

  TFile *distfile=TFile::Open(distortionFilename.Data(),"READ");
  distfile->cd();

  const char fulldistname[]="rpz";
  TString sidename[2];
  sidename[0]="_1"; //neg
  sidename[1]="_2"; //pos
  
  TString histname;
  //load the plots in:
  TH3F *hFullDist[2][3];
  //TH3F *hDiffDist3D[3];
  for (int side=0;side<2;side++){
    for (int i=0;i<3;i++){
      histname=TString::Format("hclosure_%c%s",fulldistname[i],sidename[side].Data());
      printf("seeking %s\n",histname.Data());
      hFullDist[side][i]=(TH3F*)(distfile->Get(histname));
      //hDiffDist3D[side][i]=(TH3F*)(distfile->Get(histname));
    }
  }

  int nph=hFullDist[0][0]->GetXaxis()->GetNbins();
  int nrh=hFullDist[0][0]->GetYaxis()->GetNbins();
  int nzh=hFullDist[0][0]->GetZaxis()->GetNbins()*2;
  float pih=hFullDist[0][0]->GetXaxis()->GetXmin();
  float pfh=hFullDist[0][0]->GetXaxis()->GetXmax();
  float rih=hFullDist[0][0]->GetYaxis()->GetXmin();
  float rfh=hFullDist[0][0]->GetYaxis()->GetXmax();
  float inputzih=hFullDist[0][0]->GetZaxis()->GetXmin();
  float inputzfh=hFullDist[0][0]->GetZaxis()->GetXmax();

  float zih=-fmax(abs(inputzih),abs(inputzfh));
  float zfh=fmax(abs(inputzih),abs(inputzfh));

  float deltar=(rfh-rih)/nrh;
  float deltap=(pfh-pih)/nph;
  float deltaz=(zfh-zih)/(nzh);


    //monitor plots, and the position that plot monitors at:
  float zpos=zin;
  if (zpos>1e6) zpos=zfh/2.0;//(nzh/2+0.5)*deltaz+zih;
  TVector3 pos((nrh/2+0.5)*deltar+rih,0,zpos);
  pos.SetPhi((nph/2+0.5)*deltap+pih);
  float posphi=((nph/2+0.5)*deltap+pih);//make sure we don't have a wrapping issue.

  int selectedzbin=0;
  int selectedzside=0;
 
  int zbin=(pos.Z()-zih)/(deltaz);
  
  //printf("slicing plots at pos=(%f,%f,%f)(rpz) (side=%d, histbins=(%d,%d,%d)\n",

  
  int xi[3]={hFullDist[0][0]->GetYaxis()->FindBin(pos.Perp()),hFullDist[0][0]->GetXaxis()->FindBin(posphi),zbin};
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
  for (int i=0;i<3;i++){
    //loop over which axis of the distortion to read
    for (int ax=0;ax<3;ax++){
      //loop over which plane to work in
      hIntDist[ax][i]=new TH2F(TString::Format("hIntDist%c_%c%c",axname[i],axname[ax+1],axname[ax+2]),
			       TString::Format("%c component of closure residual in  %c%c plane at %c=%2.3f;%c;%c",
					       axname[i],axname[ax+1],axname[ax+2],axname[ax],axval[ax],axname[ax+1],axname[ax+2]),
			       axn[ax+1],axbot[ax+1],axtop[ax+1],
			       axn[ax+2],axbot[ax+2],axtop[ax+2]);

      hSummableDist2D[ax][i]=new TH2F(TString::Format("hSummableDist2D%c_%c",axname[i],axname[ax]),
			       TString::Format("%c component of closure residual vs %c, for all %c and %c;%c;#delta (um)",
					       axname[i],axname[ax],axname[ax+1],axname[ax+2],axname[ax]),
			       axn[ax],axbot[ax],axtop[ax],
				      100,-500,500);
    }
    hRDist[i]=new TH1F(TString::Format("hRDist%c",axname[i]),
		       TString::Format("%c component of closure residual vs r with %c=%2.3f and %c=%2.3f;r(cm);#delta (cm)",
				       axname[i],axname[1],axval[1],axname[2],axval[2]),
		       axn[0],axbot[0],axtop[0]);
    hRDist2D[i]=new TH2F(TString::Format("hRDist2D%c",axname[i]),
		       TString::Format("%c component of closure residual vs r at %c=%2.3f;r(cm);#delta (um)",
				       axname[i],axname[2],axval[2]),
		       axn[0],axbot[0],axtop[0],100,-500,500);
    hRModDist[i]=new TH1F(TString::Format("hRModDist%c",axname[i]),
		       TString::Format("#Delta %c (not %c-hat) from closure residual vs r with %c=%2.3f and %c=%2.3f;r(cm);#delta (cm)",
				       axname[i],axname[i],axname[1],axval[1],axname[2],axval[2]),
		       axn[0],axbot[0],axtop[0]);
  }

  
  //we want to loop over the entire region to be mapped, but we also need to include
  //one additional bin at each edge, to allow the mc drift code to interpolate properly.
  //hence we count from -1 to n+1, and manually adjust the position in those edge cases
  //to avoid sampling nonphysical regions in r and z.  the phi case is free to wrap as
  // normal.

  //note that we apply the adjustment to the particle position (inpart) and not the plotted position (partR etc)
  int ir,ip,iz;
  unsigned long el=0;
  float distortR,distortP,distortZ;
  float diffDistortR,diffDistortP,diffDistortZ;
  float partR,partP,partZ;
  float distort[3];//convenient way of holding r,p,z components in loopable form
  float part[3];//convenient way of holding r,p,z coordinates in loopable form
  for (ir=0;ir<nrh;ir++){
    partR=(ir+0.5)*deltar+rih;
    for (ip=0;ip<nph;ip++){
      partP=(ip+0.5)*deltap+pih;
      for (iz=0;iz<nzh;iz++){
	partZ=(iz)*deltaz+zih; //start us at the EDGE of the z bin, not the center?
	partZ+=0.5*deltaz; //move to center of histogram bin.

	int side=(partZ>0);//if we're below zero, we're on the low side.
	int bin=hFullDist[side][0]->FindBin(partP,partR,partZ);//take advantage of the bin symmetry of the plots.

	//printf("iz=%d, side=%d, zcoord=%2.2f, bin=%d\n",iz,side, partZ,  hFullDist[side][0]->GetZaxis()->FindBin(partZ));

	//printf("part=(rpz)(%f,%f,%f),distortP=%f\n",partP,partR,partZ,distortP);
	distort[0]=distortR=hFullDist[side][0]->GetBinContent(bin);
	distort[1]=distortP=hFullDist[side][1]->GetBinContent(bin);
	distort[2]=distortZ=hFullDist[side][2]->GetBinContent(bin);

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
	    deltaR=sqrt((partR+distortR)*(partR+distortR)+distortP*distortP)-partR;
	    deltaP=atan2(distortP,partR+deltaR)*partR;
	    hRModDist[0]->Fill(partR,deltaR);
	    hRModDist[1]->Fill(partR,deltaP);
	    hRModDist[2]->Fill(partR,distortZ);
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
      if (MAX_PLOT_DISTORTION>0){
	hIntDist[ax][i]->SetMaximum(MAX_PLOT_DISTORTION);
	hIntDist[ax][i]->SetMinimum(-1.*MAX_PLOT_DISTORTION);
      }
      hIntDist[ax][i]->Draw("colz");
    }
    c->cd(i*4+4);
    if (MAX_PLOT_DISTORTION>0){
      hRDist[i]->SetMaximum(MAX_PLOT_DISTORTION);
      hRDist[i]->SetMinimum(-1.*MAX_PLOT_DISTORTION);
    }
    hRDist[i]->SetStats(0);
    hRDist[i]->SetFillColor(kRed);
    hRDist[i]->Draw("hist");
  }
  textpad->cd();
  float texpos=0.9;float texshift=0.12;
  TLatex *tex=new TLatex(0.0,texpos,"Fill Me In");
  tex->SetTextSize(texshift*0.8);
  tex->DrawLatex(0.05,texpos,"Post-hoc slices of integral distortion");texpos-=texshift;
 //tex->DrawLatex(0.05,texpos,Form("Drift Field = %2.2f V/cm",GetNominalE()));texpos-=texshift;
  tex->DrawLatex(0.05,texpos,Form("Drifting grid of (rp)=(%d x %d) electrons with steps per file",nrh,nph));texpos-=texshift;
  tex->DrawLatex(0.05,texpos,Form("Lookup per file: %s",distortionFilename.Data()));texpos-=texshift;
  tex->DrawLatex(0.05,texpos,Form("Gas per file: %s",distortionFilename.Data()));texpos-=texshift;
  texpos=0.9;
 
  
  
  canvas->cd();
  c->Draw();
  canvas->cd();
  textpad->Draw();
  canvas->SaveAs(summaryFilename.Data());
  canvas->SaveAs(gifFilename.Data());



  


  
  
canvas=new TCanvas("cmod","modified differentials",1200,800);
  //take 10 of the bottom of this for data?
  canvas->cd();
  c=new TPad("cplotsmod","distortion modified plots",0,0.2,1,1);
  canvas->cd();
  textpad=new TPad("ctextmod","distortion modified plots",0,0.0,1,0.2);
  c->Divide(3,3);
  gStyle->SetOptStat();
  for (int i=0;i<3;i++){
    //component
    c->cd(i+1);
    hRDist[i]->SetStats(0);
    hRDist[i]->SetFillColor(kRed);
    hRDist[i]->Draw("hist");
    c->cd(i+1+3);
    hRModDist[i]->SetStats(0);
    hRModDist[i]->SetFillColor(kRed);
    hRModDist[i]->Draw("hist");
    c->cd(i+1+6);
    TH1F *ratio=new TH1F(*(hRModDist[i]));
    ratio->Divide(hRDist[i]);
    ratio->SetTitle("Ratio of distortion to #Delta;r;#delta/#Delta");
    ratio->Draw("hist");
  }
  textpad->cd();
  tex->SetTextSize(texshift*0.8);
  tex->DrawLatex(0.05,texpos,"Post-hoc slices of distortion.  Top row: r-hat, #phi-hat, z-hat.  Bottom row #Delta r, #Delta #phi, #Delta z");texpos-=texshift;
  tex->DrawLatex(0.05,texpos,Form("Drifting grid of (rp)=(%d x %d) electrons with steps per file",nrh,nph));texpos-=texshift;
  tex->DrawLatex(0.05,texpos,Form("Settings per file: %s",distortionFilename.Data()));texpos-=texshift;
  tex->DrawLatex(0.05,texpos,Form("#Delta R=sqrt((r+#delta r)^2+#delta #phi^2)"));texpos-=texshift;
  tex->DrawLatex(0.05,texpos,Form("#Delta #phi=r * atan2(#delta #phi,r+#delta r)"));texpos-=texshift;
  texpos=0.9;
 
  
  
  canvas->cd();
  c->Draw();
  canvas->cd();
  textpad->Draw();
  canvas->SaveAs("coordinate.shifts.pdf");

  

  printf("read map from:%s.\n",distortionFilename.Data());
  printf("wrote summary to:%s.\n",summaryFilename.Data());

  outputHistFile->cd();
  for (int i=0;i<3;i++){
    for (int j=0;j<3;j++){
      hIntDist[i][j]->Write();
    }
    hRDist[i]->Write();
    hRDist2D[i]->Write();
    hRModDist[i]->Write();
  }

  for (int comp=0;comp<3;comp++){//loop over which component:
    for (int ax=0;ax<3;ax++){//loop over which axis we're filling:
      hSummableDist2D[ax][comp]->Write();
    }
  }
  
  outputHistFile->Close();
  printf("wrote hists to:%s.\n",histFilename.Data());

  return;
}
