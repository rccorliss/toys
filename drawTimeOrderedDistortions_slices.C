

void drawTimeOrderedDistortions_slices(const char* filebase, int event=0, float zin=1e8){
  //reads in a triplet of 3D integral distortion maps, and writes lots of plots, as well as some summary histograms to various output files.
  const float MAX_PLOT_DISTORTION=0.05;
  std::vector<TString> histname;
 histname.push_back("hIntDistortionR_negz");
  histname.push_back("hIntDistortionPhi_negz");
  histname.push_back("hIntDistortionZ_negz");
  histname.push_back("hIntDistortionR_posz");
  histname.push_back("hIntDistortionPhi_posz");
  histname.push_back("hIntDistortionZ_posz");
  
  TString distortionFilename;
  distortionFilename=filebase;
  TString basename=filebase;
  basename.ReplaceAll(".root","");
  TString summaryFilename;
  summaryFilename.Form("%s.eve%d.separate_distortion_summary.pdf",basename.Data(),event);
  TString gifFilename=Form("%s.eve%d.separate_distortion_summary.gif",basename.Data(),event);
  TString histFilename=Form("%s.eve%d.separate_distortion_summary.hist.root",basename.Data(),event);
  TFile *outputHistFile=TFile::Open(histFilename.Data(),"RECREATE");

  TFile *distfile=TFile::Open(distortionFilename.Data(),"READ");
  distfile->cd();

  //load the plots in:
  //for the TimeDists ttree:
  //load the ttree:
  TTree * tree=NULL;
  tree=(TTree*)(distfile->Get("TimeDists"));
  if (!tree) {
    printf("could not open TimeDists in file: %s\nExiting.\n",distortionFilename.Data());
    return;
  }
  //create blanks to load the events in:
  TH3F* hFullDist[histname.size()];
  for (int i=0;i<histname.size();i++){
    hFullDist[i]=new TH3F(Form("temphist%d",i),Form("temphist%d",i),10,0,10,20,0,20,30,0,30);
    tree->SetBranchAddress(histname[i].Data(),&(hFullDist[i]));
  }

  //now load the event of interest:
  tree->GetEntry(event);

  //for single event files:
  /*
  TH3F *hFullDist[3];
  for (int i=0;i<3;i++){
    printf("seeking %s\n",histname[i].Data());
    hFullDist[i]=(TH3F*)(distfile->Get(histname[i]));
   }
  */

  //store the number of bins and bin extents for each axis:
  int nph=hFullDist[0]->GetXaxis()->GetNbins();
  int nrh=hFullDist[0]->GetYaxis()->GetNbins();
  int nzh=hFullDist[0]->GetZaxis()->GetNbins()*2;
  float pih=hFullDist[0]->GetXaxis()->GetXmin();
  float pfh=hFullDist[0]->GetXaxis()->GetXmax();
  float rih=hFullDist[0]->GetYaxis()->GetXmin();
  float rfh=hFullDist[0]->GetYaxis()->GetXmax();
  float zih=hFullDist[0]->GetZaxis()->GetXmin();
  float zfh=hFullDist[0]->GetZaxis()->GetXmax();
  //do some quick checks to get a big z axis (really, we know the right answer already, 105.5cm:
  if (fabs(zih)>fabs(zfh)){
    zfh=-1*zih;
  } else {
    zih=-1*zfh;
  }
      

  float deltar=(rfh-rih)/nrh;
  float deltap=(pfh-pih)/nph;
  float deltaz=(zfh-zih)/nzh;


    //monitor plots, and the position that plot monitors at:

  // decide where to draw the slice in z:
  float zpos=zin;
  if (zpos>1e6) zpos=(nzh/2+0.5)*deltaz+zih;

  //store the slice position as a TVector
  TVector3 pos((nrh/2+0.5)*deltar+rih,0,zpos);
  pos.SetPhi((nph/2+0.5)*deltap+pih);
  float posphi=((nph/2+0.5)*deltap+pih);//make sure we don't have a wrapping issue.

  //compute the bin numbers corresponding to the slice position.
  int xi[3]={hFullDist[0]->GetYaxis()->FindBin(pos.Perp()),hFullDist[0]->GetXaxis()->FindBin(posphi),hFullDist[0]->GetZaxis()->FindBin(pos.Z())};
	     //used to be{nrh/2,nph/2,nzh/2};

  //we are about to be permuting over the r, phi, and z axis in various ways.  It's convenient when writing the loops not to have to do modulo arithmetic on every index, so instead of having helper arrays that are 3 elements, we make them with 6 elements, repeating the first three elements:
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
  float distort[3];//convenient way of holding r,p,z components in loopable form
  float part[3];//convenient way of holding r,p,z coordinates in loopable form
  for (ir=0;ir<nrh;ir++){
    partR=(ir+0.5)*deltar+rih;
    for (ip=0;ip<nph;ip++){
      partP=(ip+0.5)*deltap+pih;
      for (iz=0;iz<nzh;iz++){
	partZ=(iz)*deltaz+zih; //start us at the EDGE of the z bin, not the center?
	partZ+=0.5*deltaz; //move to center of histogram bin.

	//printf("iz=%d, zcoord=%2.2f, bin=%d\n",iz,partZ,  hIntDist[0][0]->GetYaxis()->FindBin(partZ));

	//printf("part=(rpz)(%f,%f,%f),distortP=%f\n",partP,partR,partZ,distortP);
	int offset=(partZ>0)?3:0;
	int tempbin=hFullDist[offset+0]->FindBin(partP,partR,partZ);
	distort[0]=distortR=hFullDist[offset+0]->GetBinContent(tempbin);
	distort[1]=distortP=hFullDist[offset+1]->GetBinContent(tempbin);
	distort[2]=distortZ=hFullDist[offset+2]->GetBinContent(tempbin);

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


  return;

  
  
  canvas->cd();
  c->Draw();
  canvas->cd();
  textpad->Draw();
 

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
