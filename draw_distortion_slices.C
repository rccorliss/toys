

void draw_distortion_slices(const char* filebase, float zpos=1e8){
  TString distortionFilename;
  distortionFilename=filebase;
  TString basename=filebase;
  basename.ReplaceAll(".distortion_map.hist.root","");
  TString summaryFilename;
  summaryFilename.Form("%s.separate_distortion_summary.pdf",basename.Data());
  TString gifFilename=Form("%s.separate_distortion_summary.gif",basename.Data());
 TString diffFilename;
  diffFilename.Form("%s.separate_differential_summary.pdf",basename.Data());
  TString diffgifFilename=Form("%s.separate_differential_summary.gif",basename.Data());

  TFile *distfile=TFile::Open(distortionFilename.Data(),"READ");
  distfile->cd();

  const char fulldistname[]="RPZ";
  TString histname;
  //load the plots in:
  TH3F *hFullDist[3];
  TH3F *hDiffDist3D[3];
  for (int i=0;i<3;i++){
    histname=TString::Format("hIntDistortion%c",fulldistname[i]);
    printf("seeking %s\n",histname.Data());
    hFullDist[i]=(TH3F*)(distfile->Get(histname));
    histname=TString::Format("hDistortion%c",fulldistname[i]);
    printf("seeking %s\n",histname.Data());
    hDiffDist3D[i]=(TH3F*)(distfile->Get(histname));
  }

  int nph=hFullDist[0]->GetXaxis()->GetNbins();
  int nrh=hFullDist[0]->GetYaxis()->GetNbins();
  int nzh=hFullDist[0]->GetZaxis()->GetNbins();
  float pih=hFullDist[0]->GetXaxis()->GetXmin();
  float pfh=hFullDist[0]->GetXaxis()->GetXmax();
  float rih=hFullDist[0]->GetYaxis()->GetXmin();
  float rfh=hFullDist[0]->GetYaxis()->GetXmax();
  float zih=hFullDist[0]->GetZaxis()->GetXmin();
  float zfh=hFullDist[0]->GetZaxis()->GetXmax();

  float deltar=(rfh-rih)/nrh;
  float deltap=(pfh-pih)/nph;
  float deltaz=(zfh-zih)/nzh;


    //monitor plots, and the position that plot monitors at:

  if (zpos>1e6) zpos=(nzh/2+0.5)*deltaz+zih
  TVector3 pos((nrh/2+0.5)*deltar+rih,0,zpos);
  pos.SetPhi((nph/2+0.5)*deltap+pih);
  posphi=((nph/2+0.5)*deltap+pih);//make sure we don't have a wrapping issue.
  int xi[3]={hFullDist[0]->GetYaxis()->FindBin(pos.Perp()),hFullDist[0]->GetXaxis()->FindBin(posphi),hFullDist[0]->GetXaxis()->FindBin(pos.Z())}
	     //used to be{nrh/2,nph/2,nzh/2};
  const char axname[]="rpzrpz";
  int axn[]={nrh,nph,nzh,nrh,nph,nzh};
  float axval[]={(float)pos.Perp(),(float)pos.Phi(),(float)pos.Z(),(float)pos.Perp(),(float)pos.Phi(),(float)pos.Z()};
  float axbot[]={rih,pih,zih,rih,pih,zih};
  float axtop[]={rfh,pfh,zfh,rfh,pfh,zfh};
  TH1F* hRDist[3];
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
    }
    hRDist[i]=new TH1F(TString::Format("hRDist%c",axname[i]),
		       TString::Format("%c component of int. distortion vs r with %c=%2.3f and %c=%2.3f;r(cm);#delta (cm)",
				       axname[i],axname[1],axval[1],axname[2],axval[2]),
		       axn[0],axbot[0],axtop[0]);
    hRDiffDist[i]=new TH1F(TString::Format("hRDist%c",axname[i]),
		       TString::Format("%c component of diff. distortion vs r with %c=%2.3f and %c=%2.3f;r(cm);#delta (cm)",
				       axname[i],axname[1],axval[1],axname[2],axval[2]),
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
  float partR,partP,partZ;
  for (ir=0;ir<nrh;ir++){
    partR=(ir+0.5)*deltar+rih;
    for (ip=0;ip<nph;ip++){
      partP=(ip+0.5)*deltap+pih;
      for (iz=0;iz<nzh;iz++){
	partZ=(iz)*deltaz+zih; //start us at the EDGE of the z bin, not the center?
	partZ+=0.5*deltaz; //move to center of histogram bin.

	//printf("iz=%d, zcoord=%2.2f, bin=%d\n",iz,partZ,  hIntDist[0][0]->GetYaxis()->FindBin(partZ));

	//printf("part=(rpz)(%f,%f,%f),distortP=%f\n",partP,partR,partZ,distortP);
	distortR=hFullDist[0]->GetBinContent(hFullDist[0]->FindBin(partP,partR,partZ));
	distortP=hFullDist[1]->GetBinContent(hFullDist[1]->FindBin(partP,partR,partZ));
	distortZ=hFullDist[2]->GetBinContent(hFullDist[2]->FindBin(partP,partR,partZ));
	diffDistortR=hDiffDist3D[0]->GetBinContent(hDiffDist3D[0]->FindBin(partP,partR,partZ));
	diffDistortP=hDiffDist3D[1]->GetBinContent(hDiffDist3D[1]->FindBin(partP,partR,partZ));
	diffDistortZ=hDiffDist3D[2]->GetBinContent(hDiffDist3D[2]->FindBin(partP,partR,partZ));

	//now we fill particular slices for integral visualizations:
	if(ir==xi[0]){//r slice
	  //printf("ir=%d, r=%f (pz)=(%d,%d), distortR=%2.2f, distortP=%2.2f\n",ir,partR,ip,iz,distortR,distortP);
	  hIntDist[0][0]->Fill(partP,partZ,distortR);
	  hIntDist[0][1]->Fill(partP,partZ,distortP);
	  hIntDist[0][2]->Fill(partP,partZ,distortZ);
	  hDiffDist[0][0]->Fill(partP,partZ,diffDistortR);
	  hDiffDist[0][1]->Fill(partP,partZ,diffDistortP);
	  hDiffDist[0][2]->Fill(partP,partZ,diffDistortZ);
	}
	if(ip==xi[1]){//phi slice
	  //printf("ip=%d, p=%f (rz)=(%d,%d), distortR=%2.2f, distortP=%2.2f\n",ip,partP,ir,iz,distortR,distortP);
	  hIntDist[1][0]->Fill(partZ,partR,distortR);
	  hIntDist[1][1]->Fill(partZ,partR,distortP);
	  hIntDist[1][2]->Fill(partZ,partR,distortZ);
	  hDiffDist[1][0]->Fill(partZ,partR,diffDistortR);
	  hDiffDist[1][1]->Fill(partZ,partR,diffDistortP);
	  hDiffDist[1][2]->Fill(partZ,partR,diffDistortZ);
	}
	if(iz==xi[2]){//z slice
	  //printf("iz=%d, z=%f (rp)=(%d,%d), distortR=%2.2f, distortP=%2.2f\n",iz,partZ,ir,ip,distortR,distortP);
		  
	  hIntDist[2][0]->Fill(partR,partP,distortR);
	  hIntDist[2][1]->Fill(partR,partP,distortP);
	  hIntDist[2][2]->Fill(partR,partP,distortZ);
	  hDiffDist[2][0]->Fill(partR,partP,diffDistortR);
	  hDiffDist[2][1]->Fill(partR,partP,diffDistortP);
	  hDiffDist[2][2]->Fill(partR,partP,diffDistortZ);


	  if(ip==xi[1]){//phi slices of z slices = r line at mid phi, mid z:
	    hRDist[0]->Fill(partR,distortR);	    
	    hRDist[1]->Fill(partR,distortP);
	    hRDist[2]->Fill(partR,distortZ);
	    hRDiffDist[0]->Fill(partR,diffDistortR);	    
	    hRDiffDist[1]->Fill(partR,diffDistortP);
	    hRDiffDist[2]->Fill(partR,diffDistortZ);
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


canvas=new TCanvas("cdifferent","distortion differentials",1200,800);
  //take 10 of the bottom of this for data?
  canvas->cd();
  TPad *c=new TPad("cplots2","distortion differential plots",0,0.2,1,1);
  canvas->cd();
  TPad *textpad=new TPad("ctext2","distortion differential plots",0,0.0,1,0.2);
  c->Divide(4,3);
  gStyle->SetOptStat();
  for (int i=0;i<3;i++){
    //component
    for (int ax=0;ax<3;ax++){
      //plane
      c->cd(i*4+ax+1);
      gPad->SetRightMargin(0.15);
      hDiffDist[ax][i]->SetStats(0);
      hDiffDist[ax][i]->Draw("colz");
    }
    c->cd(i*4+4);
    hRDiffDist[i]->SetStats(0);
    hRDiffDist[i]->SetFillColor(kRed);
    hRDiffDist[i]->Draw("hist");
  }
  textpad->cd();
  float texpos=0.9;float texshift=0.12;
  TLatex *tex=new TLatex(0.0,texpos,"Fill Me In");
  tex->SetTextSize(texshift*0.8);
  tex->DrawLatex(0.05,texpos,"Post-hoc slices of differential distortion");texpos-=texshift;
  tex->DrawLatex(0.05,texpos,Form("Drifting grid of (rp)=(%d x %d) electrons with steps per file",nrh,nph));texpos-=texshift;
  tex->DrawLatex(0.05,texpos,Form("Lookup per file: %s",distortionFilename.Data()));texpos-=texshift;
  tex->DrawLatex(0.05,texpos,Form("Gas per file: %s",distortionFilename.Data()));texpos-=texshift;
  texpos=0.9;
 
  
  
  canvas->cd();
  c->Draw();
  canvas->cd();
  textpad->Draw();
  canvas->SaveAs(diffFilename.Data());
  canvas->SaveAs(diffgifFilename.Data());


  

  printf("read map from:%s.\n",distortionFilename.Data());
  printf("wrote summary to:%s.\n",summaryFilename.Data());
  return;
}
