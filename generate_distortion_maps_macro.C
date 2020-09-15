/*
digital_current_macro_alice started as code to model the alice-specific TPC properties, to match to work the ALICE group did.  
it turned out to be a convenient point to compare sPHENIX performance, so it became ill-named.

This code loads an AnnularFieldSim model of the TPC, computing spacecharge distortions and using those to propagate particles from arbitrary points in the volume to a specified z-plane.  It assumes propagation or back-propagation based on the sign of (destination-start).  Details of the simulation are set in the first section, and details of the particles to be propagated in the second.

 */



#include "AnnularFieldSim.h"
#include "TTree.h" //this prevents a lazy binding issue and/or is a magic spell.
R__LOAD_LIBRARY(.libs/libfieldsim)

//place test charge  and plot how charges distort.
void TestChargeSign(AnnularFieldSim *t);

void generate_distortion_maps_macro(int reduction=0, bool loadOutputFromFile=false, const char* fname="pre-hybrid_fixed_reduction_0.ttree.root"){

  if (loadOutputFromFile) printf("loading out1 vectors from %s\n",fname);

  TTime now, start;
  start=now=gSystem->Now();
  printf("the time is %lu\n",(unsigned long)now);


  //step 1:  specify the physical parameters, scales, etc of the model, either ALICE or sPHENIXL
  
 
  //load the ALICE TPC space charge model
  const float tpc_rmin=83.5;//cm
  const float tpc_rmax=254.5;//cm
  float tpc_deltar=tpc_rmax-tpc_rmin;
  const float tpc_z=249.7;//cm
  const float tpc_cmVolt=-99930; //V
  const float tpc_driftVel=2.58*1e6;//cm per s
  const float tpc_magField=0.5;//T
  const char detgeoname[]="alice";

  const char scmapfilename[]="HeuristicSc_ALICE.root"; //my ALICE heuristic model matching Carlos
  const char scmaphistname[]="heuristic";
  //const double tpc_chargescale=1e-15;//the heuristic map has charge in fC/cm^3.  Multiply bin*tpc_chargescale*vol to get Q in C.
  
  //const char scmapfilename[]="InputSCDensityHistograms_8000events.root";
  //const char scmaphistname[]="inputSCDensity3D_8000_avg";
  //const double tpc_chargescale=8.85e-14;//their hist. has charge in units of C/cm^3 /eps0.  This is eps0 in (V*cm)/C units so that I can multiple by the volume in cm^3 to get Q in C.  

  /*
  //load the sPHENIX space charge model parameters
  const float tpc_rmin=20.0;
  const float tpc_rmax=78.0;
  float tpc_deltar=tpc_rmax-tpc_rmin;
  const float tpc_z=105.5;
  const float tpc_cmVolt=-400*tpc_z; //V =V_CM-V_RDO -- volts per cm times the length of the drift volume.
  //const float tpc_magField=0.5;//T -- The old value used in carlos's studies.
  //const float tpc_driftVel=4.0*1e6;//cm per s  -- used in carlos's studies
  const float tpc_driftVel=8.0*1e6;//cm per s  -- 2019 nominal value
  //const float tpc_magField=1.4;//T -- 2019 nominal value
  const float tpc_magField=-1.5;//T -- 2019 nominal value in the fieldmap
  const char detgeoname[]="sphenix";

  //heuristic chargescale:
  //const double tpc_chargescale=1e-15;//the heuristic map has charge in fC/cm^3.  Multiply bin*tpc_chargescale*vol to get Q in C.

  //Each unit of '1 ion' is tpc_chargescale coulombs.
  // double tpc_chargescale=1.6e-19/(100*100*100);//evgeny hists have charge in units of ions/bin, so we need to multiply by the electric charge of an electron to get out C/m^3 and then divde by 100^3 to get C/cm^3.  
  //const double tpc_chargescale=1.6e-19;//my older hists have ions/cm^3, so we need to multiply by the electric charge of an electron to get out C/cm^3.  
  */
  
   //step 2: specify the parameters of the field simulation.  Larger numbers of bins will rapidly increase the memory footprint and compute times.
  //there are some ways to mitigate this by setting a small region of interest, or a more parsimonious lookup strategy, specified when AnnularFieldSim() is actually constructed below.
  int nr=20;//26;//10;//24;//159;//159 nominal
  int nr_roi_min=0;
  int nr_roi=nr;//10;
  int nr_roi_max=nr_roi_min+nr_roi;
  int nphi=16;//40;//38;//360;//360 nominal
  int nphi_roi_min=0;
  int nphi_roi=nphi;//38;
  int nphi_roi_max=nphi_roi_min+nphi_roi;
  int nz=20;//40;//62;//62 nominal
  int nz_roi_min=0;
  int nz_roi=nz;
  int nz_roi_max=nz_roi_min+nz_roi;

  

  //step 4:  create the fieldsim object.  different choices of the last few arguments will change how it builds the lookup table spatially, and how it loads the spacecharge.  The various start-up steps are exposed here so they can be timed in the macro.
  //To set rossegger vs trivial greens functions you either call or do not call load_rossegger().
  start=now;
  AnnularFieldSim *tpc=
    new  AnnularFieldSim(tpc_rmin,tpc_rmax,tpc_z,
			 nr, nr_roi_min,nr_roi_max,1,2,
			 nphi,nphi_roi_min, nphi_roi_max,1,2,
			 nz, nz_roi_min, nz_roi_max,1,2,
			 tpc_driftVel, AnnularFieldSim::PhiSlice, AnnularFieldSim::FromFile);
  tpc->UpdateEveryN(10);//show reports every 10%.
  
  now=gSystem->Now();
  printf("created sim obj.  the dtime is %lu\n",(unsigned long)(now-start));
  start=now;






  

    //load the field maps, either flat or actual maps
  //potential is V=(tpc_CM+400V/cm*z) thus grad=+400V/cm and field =-grad=-400V/cm
  tpc->setFlatFields(tpc_magField,tpc_cmVolt/tpc_z);  
  char field_string[200];
  sprintf(field_string,"flat_B%2.1f_E%2.1f",tpc_magField,tpc_cmVolt/tpc_z);

  if (0){
    printf("loading external fieldmaps\n");
    tpc->loadBfield("sPHENIX.2d.root","fieldmap");//soon:, (1.4)/(-1.5));
    tpc->loadEfield("externalEfield.ttree.root","fTree");
    sprintf(field_string,"realE_B%2.1f_E%2.1f",tpc_magField,tpc_cmVolt/tpc_z);
  }
    now=gSystem->Now();
  printf("set fields.  the dtime is %lu\n",(unsigned long)(now-start));
  start=now;




  //load the greens functions:
  char lookup_string[200];
  sprintf(lookup_string,"ross_phi1_%s_phislice_lookup_r%dxp%dxz%d",detgeoname,nr,nphi,nz);
  char lookupFilename[200];
  sprintf(lookupFilename,"%s.root",lookup_string);
  TFile *fileptr=TFile::Open(lookupFilename,"READ");

  if (!fileptr){ //generate the lookuptable
  //to use the full rossegger terms instead of trivial free-space greens functions, uncomment the line below:
    tpc->load_rossegger();
    now=gSystem->Now();
    printf("loaded rossegger greens functions. the dtime is %lu\n",(unsigned long)(now-start));
    start=now;
    tpc->populate_lookup();
    tpc->save_phislice_lookup(lookupFilename);

  } else{ //load it from a file
    fileptr->Close();
    tpc->load_phislice_lookup(lookupFilename);
  }
    now=gSystem->Now();
  printf("populated lookup.  the dtime is %lu\n",(unsigned long)(now-start));
  start=now;



  
 
  
   //set a point we will use to look at the field slices:
  TVector3 pos=0.5*(tpc->GetOuterEdge()+tpc->GetInnerEdge());;
  pos.SetPhi(3.14159);

  
  //filename should be:  spacecharge+fieldtype+greenlookuptype-and-dimensions
  char distortionFilebase[200];
  // Evgeny's stuff:
  /*
 const int nfiles=1;
  const int nhistsper=1;
   const float spacecharge_cm_per_axis_unit=100;
   double tpc_chargescale=1.6e-19;//evgeny hists have charge in units of ions/bin, so we need to multiply by the electric charge of an electron to get out C.
   bool isChargeDensity=false; //and tell the load_ function not to multiply by volume.

   char *scbasename[]={"Cartesian","evgeny/mapFile_bX734587_bias0","Smooth.50kHz","Single.50kHz"};
   char *scfilename[]={"evgeny_sept/Summary_bX1508071_0_10_events.root","outputFile_15kHz_G4Hits_sHijing_0-12fm_002000_008000_bX734587_bias10.root", "evgeny/mapFile_bX734587_bias0.root","Smooth.50kHz.root","BeamXingNBeams.root"};
   char *schistname[]={"h_Charge_evt_0","h_Charge_0",
		      "h_Charge_1",
		      "h_Charge_2",
		      "h_Charge_3",
		      "h_Charge_4",
		      "h_Charge_5",
		      "h_Charge_6",
		      "h_Charge_7",
		      "h_Charge_8",
		      "h_Charge_9"};
  
  
  //char *scbasename[]={"Smooth.50kHz","Single.50kHz"};
  //char *scfilename[]={"Smooth.50kHz.root","BeamXingNBeams.root"};
  //char *schistname[]={"sphenix_minbias_average","sphenix_minbias_charge"};
  //
  */
  //Heuristics:
  
  
   const int nfiles=1;
   const int nhistsper=1;
   const float spacecharge_cm_per_axis_unit=1;
   bool isChargeDensity=true;
   double tpc_chargescale=1e-15;//our plot is in fC, so multiply by 1e-15 to get coulombs

  char *scbasename[]={"HeuristicSc_ALICE"};
  char *scfilename[]={"HeuristicSc_ALICE.root"};
  char *schistname[]={"heuristic"};
  

  //Test charge in flat field, for sanity:
  /*
     const int nfiles=1;
   const int nhistsper=1;
   tpc->add_testcharge(65, -3.102, 50, 4e-12);

  char *scbasename[]={"Testcharge"};
  char *schistname[]={"1e12CoulumbAtMidpoint"};
  */

  double totalQ=0;
  for (int i=0;i<nfiles;i++){
    for (int ihist=0;ihist<nhistsper;ihist++){
      tpc->load_spacecharge(scfilename[i],schistname[ihist],0,tpc_chargescale,spacecharge_cm_per_axis_unit, isChargeDensity);
      printf("Sanity check:  Q has %d elements and dim=%d\n",tpc->q->Length(), tpc->q->dim);
      for (int k=0;k<tpc->q->Length();k++){
	totalQ+=*(tpc->q->GetFlat(k));
      }
      printf("Sanity check:  Total Q in reported region is %E C\n",totalQ);
      tpc->populate_fieldmap();
      //for (int j=0;j<nscales;j++){
	printf("%s file has %s hist.  field=%s, lookup=%s. no scaling.\n",
	       scbasename[i],schistname[i],field_string,lookup_string);
	//tpc->SetDistortionScaleRPZ(scale[j],scale[j],scale[j]);
	//printf("scaled.\n");
	sprintf(distortionFilebase,"%s.hist%d.%s.%s",scbasename[i],ihist,field_string,lookup_string);
	printf("filebase=%s\n",distortionFilebase);
	tpc->GenerateDistortionMaps(distortionFilebase,2,2,2,1,true);
	printf("distortions mapped.\n");
	tpc->PlotFieldSlices(distortionFilebase,pos);
	tpc->PlotFieldSlices(distortionFilebase,pos,'B');
	printf("field mapped.\n");
	//save fieldslices too
	// }
    }
  }

  return;
}

void TestChargeSign(AnnularFieldSim *t){
  printf("Testing Charge Sign\n");
  TVector3 inner=t->GetInnerEdge();
  TVector3 outer=t->GetOuterEdge();
  TVector3 mid=0.5*(inner+outer);
  TVector3 span=(outer-inner);
  float rstart=inner.Perp();
  float rend=outer.Perp();
  float zend=mid.Z()+0.45*span.Z();//travel to higher z
  float zstart=mid.Z()-0.45*span.Z();//start at lower z
  TVector3 cpos(mid.Perp(),0,0);
  float phipos=TMath::Pi();
  cpos.SetPhi(phipos);
  float charge=-1e-12;//16 picocoulomb of charge = 1e8 electrons
  t->add_testcharge(cpos.Perp(),cpos.Phi(),cpos.Z(),charge);
  t->populate_fieldmap();

  //PlotFieldSlices("testcharge",t,cpos);

  int nr=47;
  float deltar=span.Perp()/nr;
  int np=47;
  float deltap=6.28/np;
  int validToStep;
  float z[]={zstart,zend};
  int nDriftSteps=500;
  int nSampleSteps=47;

  TCanvas *c;
  

  TH2F *hDistortR[4];
  TH2F *hDistortP[4];
  TH2F *hDistortC[2];
  TH1F *hProfileSum[2];

  TVector3 inpart(1,0,0);
  TVector3 outpart,distort;
  float partR,partP,distortR,distortP;
  c=new TCanvas("ctestdrift","test drift data",1200,800);
    gStyle->SetOptStat();

  c->Divide(2,1);//raw data on the left, summary on the right
  c->cd(1)->Divide(2,3);//two columns, three rows of raw data
  c->cd(2)->Divide(1,3);//three rows of summary with differing columns
  c->cd(2)->cd(1)->Divide(3,1);//three plots about R
  c->cd(2)->cd(2)->Divide(3,1);//three plots about Phi
  //and cd(3) is text.
  
  for (int iz=0;iz<2;iz++){
    hDistortR[iz]=new TH2F(Form("hDistortR%d",iz),
			   Form("R Distortion Drifting from z=%2.1fcm to z=%2.1fcm;phi;r;#Delta R (cm)",z[iz],z[1-iz]),
			   np,0,6.28,nr,rstart,rend);
    hDistortP[iz]=new TH2F(Form("hDistortP%d",iz),
			   Form("Phi Distortion Drifting from z=%2.1fcm to z=%2.1fcm;phi;r;#Delta#Phi (cm)",z[iz],z[1-iz]),
			   np,0,6.28,nr,rstart,rend);
    hDistortC[iz]=new TH2F(Form("hDistortC%d",iz),
			   Form("%% Successful Steps from z=%2.1fcm to z=%2.1fcm;phi;r;#Delta#Phi (cm)",z[iz],z[1-iz]),
			   np,0,6.28,nr,rstart,rend);
    for (int ir=0;ir<nr;ir++){
      inpart.SetPerp((ir+0.1)*deltar+rstart);
      partR=inpart.Perp();
      for (int ip=0;ip<np;ip++){
	inpart.SetPhi((ip+0.1)*deltap);
	partP=inpart.Phi();
	if (partP<0) partP+=TMath::TwoPi();
      
	inpart.SetZ(z[iz]);
	outpart=t->swimToInSteps(z[1-iz],inpart,nDriftSteps,true, &validToStep);
	distort=outpart-inpart;
	distortR=outpart.Perp()-inpart.Perp();
	distort.RotateZ(-inpart.Phi());//rotate so that that is on the x axis
	distortP=distort.Y();//the phi component is now the y component.
	hDistortR[iz]->Fill(partP,partR,distortR);
	hDistortP[iz]->Fill(partP,partR,distortP);
	hDistortC[iz]->Fill(partP,partR,(1.0*validToStep)/nDriftSteps*100);
      }
    }

    c->cd(1)->cd(iz+1);
    hDistortR[iz]->SetStats(0);
    hDistortR[iz]->Draw("colz");
    c->cd(1)->cd(2+iz+1);
    hDistortP[iz]->SetStats(0);

    hDistortP[iz]->Draw("colz");
    c->cd(1)->cd(4+iz+1);
    hDistortC[iz]->SetStats(0);
    hDistortC[iz]->Draw("colz");
  }


    printf("Finished Test Drifts\n");



  
  c->cd(2)->cd(1)->cd(1)->SetLogz();
  hDistortR[2]=(TH2F*)hDistortR[1]->Clone();
  hDistortR[2]->SetTitle("Ratio of R distortion in Forward and Backward drift;phi,r");
  hDistortR[2]->Divide(hDistortR[0]);
  hDistortR[2]->Scale(-1);
  hDistortR[2]->SetStats(0);
  hDistortR[2]->Draw("colz");
  c->cd(2)->cd(2)->cd(1)->SetLogz();
  hDistortP[2]=(TH2F*)hDistortP[1]->Clone();
  hDistortP[2]->SetTitle("Ratio of Phi distortion in Forward and Backward drift;phi,r");
  hDistortP[2]->Divide(hDistortP[0]);
  hDistortP[2]->Scale(-1);
  hDistortP[2]->SetStats(0);
  hDistortP[2]->Draw("colz");

  c->cd(2)->cd(1)->cd(2);
  hDistortR[3]=(TH2F*)hDistortR[1]->Clone();
  hDistortR[3]->SetTitle("Sum of R distortion in Forward and Backward drift;phi,r");
  hDistortR[3]->Add(hDistortR[0]);
  hDistortR[3]->SetStats(0);
  hDistortR[3]->Draw("colz");
  c->cd(2)->cd(2)->cd(2);
  hDistortP[3]=(TH2F*)hDistortP[1]->Clone();
  hDistortP[3]->SetTitle("Sum of Phi distortion in Forward and Backward drift;phi,r");
  hDistortP[3]->Add(hDistortP[0]);
  hDistortP[3]->SetStats(0);
  hDistortP[3]->Draw("colz");

  float maxdistort[2];
  maxdistort[0]=max(hDistortR[3]->GetMaximum(),-1*hDistortR[3]->GetMinimum());
  maxdistort[1]=max(hDistortP[3]->GetMaximum(),-1*hDistortP[3]->GetMinimum());
  
  for (int i=0;i<2;i++){
    hProfileSum[i]=new TH1F(Form("hProfileSum%d",i),"Histogram of Cell Contents for Drift Steps >99%",100,-1.5*maxdistort[i],1.5*maxdistort[i]);
  }

  for (int i=0;i<hDistortR[3]->GetNbinsX();i++){
    for (int j=0;j<hDistortR[3]->GetNbinsY();j++){
      if (hDistortC[0]->GetBinContent(i+1,j+1)>99 && hDistortC[1]->GetBinContent(i+1,j+1)>99){
      hProfileSum[0]->Fill(hDistortR[3]->GetBinContent(i+1,j+1));
      hProfileSum[1]->Fill(hDistortP[3]->GetBinContent(i+1,j+1));
      }
    }
  }
  
  c->cd(2)->cd(1)->cd(3)->SetLogy();
  hProfileSum[0]->Draw();
  c->cd(2)->cd(2)->cd(3)->SetLogy();
  hProfileSum[1]->Draw();

  
  c->cd(2)->cd(3);
  float texpos=0.9;float texshift=0.08;
  TLatex *tex=new TLatex(0.0,texpos,"Fill Me In");
  tex->DrawLatex(0.0,texpos,Form("Test Charge: %2.2E C @ (rpz)=(%2.1fcm,%2.1frad,%2.1fcm)",charge,cpos.Perp(),cpos.Phi(),cpos.Z()));texpos-=texshift;
  tex->DrawLatex(0.0,texpos,Form("Drift Field = %2.2f V/cm",t->GetNominalE()));texpos-=texshift;
  tex->DrawLatex(0.0,texpos,Form("Drifting grid of (rp)=(%d x %d) electrons in %d steps",nr,np,nDriftSteps));texpos-=texshift;
  tex->DrawLatex(0.0,texpos,t->GetLookupString());texpos-=texshift;
  return;
}


