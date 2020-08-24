//#include "FieldSim.h"

/*
digital_current_macro_alice started as code to model the alice-specific TPC properties, to match to work the ALICE group did.  
it turned out to be a convenient point to compare sPHENIX performance, so it became ill-named.

This code loads an AnnularFieldSim model of the TPC, computing spacecharge distortions and using those to propagate particles from arbitrary points in the volume to a specified z-plane.  It assumes propagation or back-propagation based on the sign of (destination-start).  Details of the simulation are set in the first section, and details of the particles to be propagated in the second.

 */



#include "AnnularFieldSim.h"
R__LOAD_LIBRARY(.libs/libfieldsim)

//place test charge  and plot how charges distort.
void TestChargeSign(AnnularFieldSim *t);

//plot the fields at a certain point.
void PlotFieldSlices(const char* filebase, AnnularFieldSim *t, TVector3 pos);


//place test electrons along a fixed grid and drift them one grid-length in z to calculate local distortions
void GenerateAndSaveDistortionMap(const char* filebase,AnnularFieldSim *t,int nphi,float pi,float pf,int nr,float ri,float rf,int nz,float zi,float zf);
//reads the charge object from the tpc and saves it to file.  Useful to make sure things are being filled correctly.
void SaveChargeAndProjections(const char* filename,AnnularFieldSim *t,int nphi,float pi,float pf,int nr,float ri,float rf,int nz,float zi,float zf);

//reads the field object from the tpc roi and saves it to file.  Useful to make sure the field is sane or to check that propagation makes sense.
void SaveField(const char* filename,AnnularFieldSim *t,int pi,int pf,int ri, int rf, int zi, int zf);



void digital_current_macro_alice(int reduction=0, bool loadOutputFromFile=false, const char* fname="pre-hybrid_fixed_reduction_0.ttree.root"){

  if (loadOutputFromFile) printf("loading out1 vectors from %s\n",fname);

  TTime now, start;
  start=now=gSystem->Now();
  printf("the time is %lu\n",(unsigned long)now);


  //step 1:  specify the physical parameters, scales, etc of the model, either ALICE or sPHENIXL
  
  /*
  //load the ALICE TPC space charge model
  const float tpc_rmin=83.5;
  const float tpc_rmax=254.5;
  float tpc_deltar=tpc_rmax-tpc_rmin;
  const float tpc_z=249.7;
  const float tpc_driftVolt=-99930; //V
  const float tpc_driftVel=2.58*1e6;//cm per s
  const float tpc_magField=0.5;//T
  const double epsilonnaught=8.854e-12;// units of C/(V*m)
  const double eps_in_cm=epsilonnaught/100; //units of C/(V*cm)
  const double tpc_chargescale=8.85e-14;//their hist. has charge in units of C/cm^3 /eps0.  This is eps0 in (V*cm)/C units so that I can multiple by the volume in cm^3 to get Q in C.
  const char scmapfilename[]="InputSCDensityHistograms_8000events.root";
  const char scmaphistname[]="inputSCDensity3D_8000_avg";
  */

  //load the sPHENIX space charge model parameters
  const float tpc_rmin=20.0;
  const float tpc_rmax=78.0;
  float tpc_deltar=tpc_rmax-tpc_rmin;
  const float tpc_z=105.5;
  const float tpc_driftVolt=-400*tpc_z; //V =V_CM-V_RDO -- volts per cm times the length of the drift volume.
  const float tpc_magField=0.5;//T -- The old value used in pedro's studies.
  const float tpc_driftVel=4.0*1e6;//cm per s  -- used in pedro's studies

  

  //  const float tpc_driftVel=8.0*1e6;//cm per s  -- 2019 nominal value
  //const float tpc_magField=1.4;//T -- 2019 nominal value
  const double epsilonnaught=8.854e-12;// units of C/(V*m)
  const double eps_in_cm=epsilonnaught/100; //units of C/(V*cm)
  //old map const double tpc_chargescale=-0.5*1.6e-19*10*1000;//should be -19;//our hist. has charge in units of ions/cm^3, so we need to multiply by the electric charge of an electron to get out C/cm^3.  adding in a factor of 0.5 for the double-counting in the original 10khz sample, a factor of 10 to get to 100khz, and a factor 0f 1k that I don't understand
  //old map const char scmapfilename[]="G4Hits_sHijing_0-12fm_10kHz.rcc_sc.hist.root";
  const double tpc_chargescale=-1.6e-19*1000;//should be -19;//our hist. has charge in units of ions/cm^3, so we need to multiply by the electric charge of an electron to get out C/cm^3.  adding in a factor of 1k that I don't understand

  
  const char scmapfilebase[]="Smooth.50kHz"; //expects suffix '.root' after this string.
  const char scmaphistname[]="sphenix_minbias_average";




  //step 2: specify the parameters of the field simulation.  Larger numbers of bins will rapidly increase the memory footprint and compute times.
  //there are some ways to mitigate this by setting a small region of interest, or a more parsimonious lookup strategy, specified when AnnularFieldSim() is actually constructed below.
  
  //define a region of interest, in units of the intrinsic scale of the tpc histogram:
  //we will reduce these when we call the macro, but keep the full scale here so the calculations for our test grid are not changed.
  int nr=16;//10;//24;//159;//159 nominal
  int nr_roi_min=0;
  int nr_roi=16;//10;
  int nr_roi_max=nr_roi_min+nr_roi;
  int nphi=36;//38;//360;//360 nominal
  int nphi_roi_min=0;
  int nphi_roi=36;//38;
  int nphi_roi_max=nphi_roi_min+nphi_roi;
  int nz=40;//62;//62 nominal
  int nz_roi_min=0;
  int nz_roi=nz;
  int nz_roi_max=nz_roi_min+nz_roi;
  float buffer_fraction=0.5;// to stay clear of the roi bounds, what fraction of the bin size should we stop short by?


  //step 2a:  Calculate some details of the roi dimensions from the variables chosen above.  no user-serviceable parts here:
  float rmin_roi=tpc_rmin+tpc_deltar/(nr*1.0)*nr_roi_min;
  float rmax_roi=rmin_roi+tpc_deltar/nr*nr_roi;
  float phimin_roi=2*TMath::Pi()/(nphi*1.0)*nphi_roi_min;
  float phimax_roi=phimin_roi+2*TMath::Pi()/nphi*nphi_roi;
  float zmin_roi=tpc_z/(nz*1.0)*nz_roi_min;
  float zmax_roi=zmin_roi+tpc_z/nz*nz_roi;
  float rmin_roi_with_buffer=rmin_roi+tpc_deltar/(nr*1.0)*(buffer_fraction);
  float rmax_roi_with_buffer=rmax_roi-tpc_deltar/(nr*1.0)*(buffer_fraction);
  float phimin_roi_with_buffer=phimin_roi+2*TMath::Pi()/(nphi*1.0)*(buffer_fraction);
  float phimax_roi_with_buffer=phimax_roi-2*TMath::Pi()/(nphi*1.0)*(buffer_fraction);
  float zmin_roi_with_buffer=zmin_roi+tpc_z/(nz*1.0)*(buffer_fraction);
  float zmax_roi_with_buffer=zmax_roi-tpc_z/(nz*1.0)*(buffer_fraction);

  printf("r bounds are %f<%f<%f<r<%f<%f<%f\n",tpc_rmin,rmin_roi,rmin_roi_with_buffer,rmax_roi_with_buffer,rmax_roi,tpc_rmax);
  printf("phi bounds are %f<%f<%f<phi<%f<%f<%f\n",0.0,phimin_roi,phimin_roi_with_buffer,phimax_roi_with_buffer,phimax_roi,2*TMath::Pi());


  //step 2b:  The default study I've been doing for a while is to reduce the size of the simulation by a certain factor and see how that impacts the performance.  The code below adjusts the size of the simulation per the specified reduction factor (where the default value of '0' makes this skip):
  if (reduction==0){
    //do nothing
  }else if (reduction>0){
    //reduce all dimensions by the specified amount:
    nr=nr-reduction;
    nphi=nphi-reduction;
    nz=nz-reduction;
  } else if (reduction<0){
    //scale all dimensions by the specified amount in percent:
    nr=ceil(-nr*reduction*0.01);
    nphi=ceil(-nphi*reduction*0.01);
    nz=ceil(-nz*reduction*0.01);
  }
  //make sure we didn't break anything:
  if (nr<1)nr=1;
  if (nphi<1)nphi=1;
  if (nz<1)nz=1;
  if (nr_roi_max>nr)nr_roi_max=nr;
  if (nphi_roi_max>nphi)nphi_roi_max=nphi;
  if (nz_roi_max>nz)nz_roi_max=nz;




   //step 3: load the spacecharge density map from the specified histogram (in the tpc parameters, above).
  TFile *f=TFile::Open(Form("%s.root",scmapfilebase));
  TH3F* tpc_average=(TH3F*)f->Get(scmaphistname);
  now=gSystem->Now();
  printf("loaded hist.  the dtime is %lu\n",(unsigned long)(now-start));




  //step 4:  create the fieldsim object.  different choices of the last few arguments will change how it builds the lookup table spatially, and how it loads the spacecharge.  The various start-up steps are exposed here so they can be timed in the macro.
  //To set rossegger vs trivial greens functions you either call or do not call load_rossegger().
  start=now;
  AnnularFieldSim *tpc=
    new  AnnularFieldSim(tpc_rmin,tpc_rmax,tpc_z,
			 nr, nr_roi_min,nr_roi_max,1,2,
			 nphi,nphi_roi_min, nphi_roi_max,1,2,
			 nz, nz_roi_min, nz_roi_max,1,2,
			 tpc_driftVel, AnnularFieldSim::PhiSlice, AnnularFieldSim::FromFile);
  //  new AnnularFieldSim(tpc_rmin,tpc_rmax,tpc_z,9,120,9,tpc_driftVel);
   
  // dropping half-res for test: new AnnularFieldSim(tpc_rmin,tpc_rmax,tpc_z,53,18,31,tpc_driftVel);
  //full resolution is too big:  new AnnularFieldSim(tpc_rmin,tpc_rmax,tpc_z,159,360,62,tpc_driftVel);
  now=gSystem->Now();
  printf("created sim obj.  the dtime is %lu\n",(unsigned long)(now-start));
  start=now;
  tpc->setFlatFields(tpc_magField,tpc_driftVolt/tpc_z);
  char *field_string=Form("flat_B%2.1f_E%2.1f",tpc_magField,tpc_driftVolt/tpc_z);

  if (1){
    tpc->loadBfield("sPHENIX.2d.root","fieldmap");
    tpc->loadEfield("externalEfield.ttree.root","fTree");
    sprintf(field_string,"real_B%2.1f_E%2.1f",tpc_magField,tpc_driftVolt/tpc_z);
  }
    now=gSystem->Now();
  printf("set fields.  the dtime is %lu\n",(unsigned long)(now-start));
  start=now;


  //load the greens functions:
  char* lookup_string=Form("ross_phislice_lookup_r%dxp%dxz%d",nr,nphi,nz);
  char* lookupFilename=Form("%s.root",lookup_string);
  TFile *fileptr=TFile::Open(lookupFilename,"READ");
  bool lookupFileExists=(fileptr);
  
  if (!fileptr){ //generate the lookuptable
  //to use the full rossegger terms instead of trivial free-space greens functions, uncomment the line below:
    tpc->load_rossegger();
    now=gSystem->Now();
    printf("loaded rossegger greens functions. (phi set to zero) the dtime is %lu\n",(unsigned long)(now-start));
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
 
  
  //load the spacecharge:

  char sc_string[100];
  sprintf(sc_string,"zero_spacecharge");
  
  if (0){
    TestChargeSign(tpc); //adds a test charge, looks at how electrons deflect, and draws some plots.
    return;
  }

  if (0){ //load spacecharge from file
    sprintf(sc_string,scmapfilebase);
    tpc->load_spacecharge(tpc_average,0,tpc_chargescale); //(TH3F charge histogram, float z_shift in cm, float multiplier to local units)
  }
  if (0){ //load spacecharge from analytic formula
    sprintf(sc_string,"analytic_spacecharge");
    //computed the correction to get the same spacecharge as in the tpc histogram:
    //todo: make the analytic scale proportional to the tpc_chargescale.
    double tpc_analytic_scale=1.237320E-06/9.526278E-11;
    tpc->load_analytic_spacecharge(0);//tpc_analytic_scale);
  }
  

  now=gSystem->Now();

  
  printf("loaded spacecharge.  the dtime is %lu\n",(unsigned long)(now-start));
  start=now;
 tpc->populate_fieldmap();
 now=gSystem->Now();
  printf("populated fieldmap.  the dtime is %lu\n",(unsigned long)(now-start));
  start=now;
  //printf("consistency check:  integrate field along IR and OR, confirm V:\n");



  

  //some sanity checks: 
  if (reduction==0){
    //save the charge into a histogram with size and dimensions matching the native tpc object:
    //  SaveChargeAndProjections("last_macro.charge_projections.hist.root",tpc,nphi,0,6.281,nr,tpc_rmin,tpc_rmax,nz,0,tpc_z);
  }

  //SaveField("last_macro.field.hist.root",tpc,nr_roi_min,nr_roi_max,nphi_roi_min,nphi_roi_max,nz_roi_min,nz_roi_max);


  printf("fieldsim ptr=%p\n",(void*)tpc);

  
  //generate the distortions by throwing particles at each z step and letting them propagate to the next:
  //eventually this should be folded into AnnularFieldSim.


  //filename should be:  spacecharge+fieldtype+greenlookuptype-and-dimensions
  char* distortionFilebase=Form("%s.%s.%s",sc_string,field_string,lookup_string);

  
  //char* distortionFilename=Form("zero_sc.realfield.ross_phislice_lookup_r%dxp%dxz%d.distortion_map.hist.root",nr,nphi,nz);


  //save some diagnostic plots about the field:
  TVector3 pos(0.5*(rmax_roi+rmin_roi),0,0.5*(zmax_roi+zmin_roi));
  pos.SetPhi(0.5*(phimax_roi+phimin_roi));
  printf("distortionFilebase before plot=%s\n",distortionFilebase);
  PlotFieldSlices(distortionFilebase,tpc, pos);
  printf("distortionFilebase before generate=%s\n",distortionFilebase);
  //generate distortions based on those fields:
  GenerateAndSaveDistortionMap(distortionFilebase,tpc,
			       nphi_roi,phimin_roi,phimax_roi,
			       nr_roi,rmin_roi,rmax_roi,
			       nz_roi,zmin_roi,zmax_roi);

  printf("all done.\n");
  return;

  
  //define a grid of test points:
  const int divisor=100;
  const int nparticles=divisor*divisor;
  TVector3 testparticle[nparticles];
  TVector3 outparticle[nparticles];
  TVector3 outparticle2[nparticles];
  TVector3 backparticle[nparticles];
  float outx[nparticles], outy[nparticles],outz[nparticles];
  for (int i=0;i<nparticles;i++){
    int rpos=i/divisor;
    int phipos=i%divisor;
    testparticle[i].SetXYZ((rmax_roi_with_buffer-rmin_roi_with_buffer)*(rpos/(divisor*1.0))+rmin_roi_with_buffer,0,zmin_roi);
    testparticle[i].RotateZ((phimax_roi_with_buffer-phimin_roi_with_buffer)*(phipos/(divisor*1.0))+phimin_roi_with_buffer);
  }


  //for this study, we load the 'out' particles from the ttree of the full-scale version, and propagate them backward:
  if (loadOutputFromFile){
    TFile *input=TFile::Open(fname);
    TTree *inTree=(TTree*)input->Get("pTree");
    TVector3 *outFromFile;
    TVector3 *origFromFile;
    inTree->SetBranchAddress("out1",&outFromFile);
    inTree->SetBranchAddress("orig",&origFromFile);
    for (int i=0;i<inTree->GetEntries();i++){
      inTree->GetEntry(i);
      outparticle[i]=(*outFromFile)*(1/(1.0e4)); //convert back to local units.
      testparticle[i]=(*origFromFile)*(1/(1.0e4)); //convert back to local units.
    }
    input->Close();
  }

  TFile *output=TFile::Open("last_macro_output.ttree.root","RECREATE");
  TVector3 orig,outa,out1,out2,back1,back2;
  int goodSteps[2];
  TTree pTree("pTree","Particle Tree");
  pTree.Branch("orig","TVector3",&orig);
  pTree.Branch("outa","TVector3",&outa);
  pTree.Branch("out1","TVector3",&out1);
    pTree.Branch("out1N",&goodSteps[0]);
  pTree.Branch("back1","TVector3",&back1);
  pTree.Branch("back1N",&goodSteps[1]);

  
  int validToStep=-1;
  start=gSystem->Now();
  for (int i=0;i<nparticles;i++){
    if (!(i%((nparticles-1)/20))) {
       now=gSystem->Now();
       printf("(periodic progress...) test[%d]=(%f,%f,%f) to %f\t dtime=%lu\n",
	      i,testparticle[i].X(),testparticle[i].Y(),testparticle[i].Z(),zmax_roi,(unsigned long)(now-start));
      start=now;
    }



      
    orig=testparticle[i];
    if (!loadOutputFromFile)
      {
	outparticle[i]=tpc->swimToInSteps(zmax_roi,testparticle[i],60,true, &validToStep);
	outparticle2[i].SetXYZ(0,0,0);//=tpc->swimToInAnalyticSteps(zmax_roi,testparticle[i],60, &validToStep); //no arg for interpolation in the analytic swim, since we sample the exact position.
      }

    out1=outparticle[i];//not generating from the swim.=tpc->swimToInSteps(zmax_roi,testparticle[i],600,true);
    outa=outparticle2[i];//not generating from the swim.=tpc->swimToInSteps(zmax_roi,testparticle[i],600,true);
    outx[i]=outparticle[i].X();
    outy[i]=outparticle[i].Y();
    outz[i]=outparticle[i].Z();
    goodSteps[0]=validToStep;
    TVector3 delta=out1-testparticle[i];
    //printf("delta1=(%E,%E,%E)\n",delta.X(),delta.Y(),delta.Z());
    //printf("out[%d]=(%f,%f,%f)\n",i,outparticle[i].X(),outparticle[i].Y(),outparticle[i].Z());
    //drift back from the analytic position:
      back1=backparticle[i]=tpc->swimToInSteps(testparticle[i].Z(),out1,600,true,&validToStep);
    goodSteps[1]=validToStep;

    //for convenience of reading, set all of the pTree in microns, not cm:
    orig*=1e4;//10mm/cm*1000um/mm
    out1*=1e4;//10mm/cm*1000um/mm
    outa*=1e4;
    back1*=1e4;//10mm/cm*1000um/mm
    pTree.Fill();
  }
  pTree.Write();



   


  return;
  

  AnnularFieldSim *test;
  AnnularFieldSim *testgen;
  const TVector3 cyl(60,0,100);
  TVector3 ptemp(12.005,45.005,75.99);
  TVector3 ftemp,btemp;
  TTime t[10];
  Int_t isteps=22;
  int imin=8;//can't do interpolation without at least some divisions, otherwise we necessarily go out of bounds.
  int imax=imin+isteps;
  Int_t tCreate[isteps];
  Int_t tSetField[isteps];
  Int_t tLookup[isteps];
  Int_t tGenMap[isteps];
  Int_t tSwim1k[isteps];
  Int_t steps[isteps];
  float rdiff[isteps];
  float rphidiff[isteps];
  
  testgen=new AnnularFieldSim(cyl.Perp(),2*TMath::Pi(),cyl.Z(),isteps+imin-1,isteps+imin-1,isteps+imin-1,100.0/12e-6);
  testgen->setFlatFields(1.4,200);
  testgen->populate_lookup();//2-3
  testgen->q->Set(isteps+imin-2,0,0,1e-12);///that's 10^7 protons in that box.
  testgen->populate_fieldmap();
  ftemp=ptemp;
  for (int j=0;j<isteps+imin;j++){
    ftemp=testgen->swimTo(ptemp.Z()*(isteps+imin-1-j)/(isteps+imin),ftemp,true);
  }  
  for (int i=0;i<isteps;i++){
    printf("create %d\n",i);
    t[0]=gSystem->Now();
    test=new AnnularFieldSim(cyl.Perp(),2*TMath::Pi(),cyl.Z(),i+imin,i+imin,i+imin,100.0/12e-6);
    t[1]=gSystem->Now();
    printf("setFlat %d\n",i);
    test->setFlatFields(1.4,200);
    t[2]=gSystem->Now();
    printf("populate %d\n",i);
    test->populate_lookup();//2-3
    t[3]=gSystem->Now();
    printf("setQ %d\n",i);
    test->q->Set(i+imin-1,0,0,1e-12);///that's 10^7 protons in that box.
    t[4]=gSystem->Now();
    printf("make fieldmap %d\n",i);
   test->populate_fieldmap();
    t[5]=gSystem->Now();
    printf("swim %d\n",i);
    for (int n=0;n<1000;n++){
    btemp=ftemp;
    for (int j=0;j<i+imin;j++){
      btemp=test->swimTo(ptemp.Z()*((j+1)*1.0/(1.0*(i+imin))),btemp,true);
    }
    }
    t[6]=gSystem->Now();
    tCreate[i]=(long)(t[1]-t[0]);
    tSetField[i]=(long)(t[2]-t[1]);
    tLookup[i]=(long)(t[3]-t[2]);
    tGenMap[i]=(long)(t[5]-t[4]);
    tSwim1k[i]=(long)(t[6]-t[5]);
    steps[i]=i+imin;
    rdiff[i]=(btemp.Perp()-ptemp.Perp())*1e4;
    rphidiff[i]=(btemp.Perp()*btemp.Phi()-ptemp.Perp()*ptemp.Phi())*1e4;

  }

  int i=0;

  TGraph *rPhiDiff[3];
  TMultiGraph *rphi=new TMultiGraph();
  rphi->SetTitle("r and r*phi offset for varying reco grid sizes;r(um);rphi(um)");
  rPhiDiff[i]=new TGraph(isteps,rdiff,rphidiff);
  rPhiDiff[i]->SetTitle("rphi diff;size;t(ms)");
  rPhiDiff[i]->SetMarkerColor(kBlack);
  rPhiDiff[i]->SetMarkerStyle(kStar);
  rphi->Add(rPhiDiff[i]);
  i++;
  rPhiDiff[i]=new TGraph(1,rdiff,rphidiff);
  rPhiDiff[i]->SetTitle(Form("rphi diff %d^3 grid;size;t(ms)",steps[0]));
  rPhiDiff[i]->SetMarkerColor(kRed);
  rPhiDiff[i]->SetMarkerStyle(kStar);
  rphi->Add(rPhiDiff[i]);
  i++;
  rPhiDiff[i]=new TGraph(1,rdiff+isteps-1,rphidiff+isteps-1);
  rPhiDiff[i]->SetTitle(Form("rphi diff %d^3 grid;size;t(ms)",steps[isteps-1]));
  rPhiDiff[i]->SetMarkerColor(kGreen);
  rPhiDiff[i]->SetMarkerStyle(kStar);
  rphi->Add(rPhiDiff[i]);
  i++;
  rphi->Draw("AC*");
  
 




 i=0;
  
  TMultiGraph *multi=new TMultiGraph();
  multi->SetTitle("timing for various simulation steps;size;t(ms)");
  TGraph *graph[10];
  graph[i]=new TGraph(isteps,steps,tCreate);
  graph[i]->SetTitle("FieldSim::FieldSim;size;t(ms)");
  graph[i]->SetMarkerStyle(kStar);
  graph[i]->SetMarkerColor(i+1);
  multi->Add(graph[i]);
  i++;
  graph[i]=new TGraph(isteps,steps,tSetField);
  graph[i]->SetTitle("FieldSim::SetField;size;t(ms)");
  graph[i]->SetMarkerStyle(kStar);
  graph[i]->SetMarkerColor(i+1);
  multi->Add(graph[i]);
  i++;
  graph[i]=new TGraph(isteps,steps,tLookup);
  graph[i]->SetTitle("FieldSim::GenLookupTable;size;t(ms)");
  graph[i]->SetMarkerStyle(kStar);
  graph[i]->SetMarkerColor(i+1);
  multi->Add(graph[i]);
  i++;
   graph[i]=new TGraph(isteps,steps,tGenMap);
  graph[i]->SetTitle("FieldSim::GenFieldMap;size;t(ms)");
  graph[i]->SetMarkerStyle(kStar);
  graph[i]->SetMarkerColor(i+1);
  multi->Add(graph[i]);
  i++;
   graph[i]=new TGraph(isteps,steps,tSwim1k);
  graph[i]->SetTitle("FieldSim::Swim*1k;size;t(ms)");
  graph[i]->SetMarkerStyle(kStar);
  graph[i]->SetMarkerColor(i+1);
  multi->Add(graph[i]);
  i++;
  multi->Draw("AC*");  
  return;
}

//rcchere
void GenerateAndSaveDistortionMap(const char* filebase,AnnularFieldSim *t,int np,float pi,float pf,int nr,float ri,float rf,int nz,float zi,float zf){
  //scan over the tpc physical volume in nphi steps from pi to pf, and similar for the other two dimensions.
  //set a particle at those coordinates and drift it to the next coordinate in z, then save the delta in histograms.

  //by doing this per-bin, we can properly apply fractional distortions to particles that are not along a bin boundary, though this may be a minor concern.  If we integrated, we would lose this power.

  
  printf("generating distortion map...\n");
  printf("file=%s\n",filebase);
  printf("Phi:  %d steps from %f to %f\n",np,pi,pf);
  printf("R:  %d steps from %f to %f\n",nr,ri,rf);
  printf("Z:  %d steps from %f to %f\n",nz,zi,zf);
  printf("fieldsim ptr=%p\n",(void*)t);
  TFile *outf=TFile::Open(Form("%s.distortion_map.hist.root",filebase),"RECREATE");
  outf->cd();

  //for interpolation, Henry needs one extra buffer bin on each side.
  TVector3 step((rf-ri)/(1.0*nr),0,(zf-zi)/(1.0*nz));
  step.SetPhi((pf-pi)/(1.0*np));
  //so we define the histogram bounds (the 'h' suffix) to be the full range
  //plus an additional step in each direction so interpolation can work at the edges
  int nph=np+2;
  int nrh=nr+2;
  int nzh=nz+2;
  float rih=ri-step.Perp();
  float rfh=rf+step.Perp();
  float pih=pi-step.Phi();
  float pfh=pf+step.Phi();
  float zih=zi-step.Z();
  float zfh=zf+step.Z();
  

  TH3F* hDistortionR=new TH3F("hDistortionR","Per-z-bin Distortion in the R direction as a function of (r,phi,z) (centered in r,phi, z);phi;r;z",nph,pih,pfh,nrh,rih,rfh,nzh,zih,zfh);
  TH3F* hDistortionP=new TH3F("hDistortionP","Per-z-bin Distortion in the RPhi direction as a function of (r,phi,z)  (centered in r,phi, z);phi;r;z",nph,pih,pfh,nrh,rih,rfh,nzh,zih,zfh);
  TH3F* hDistortionZ=new TH3F("hDistortionZ","Per-z-bin Distortion in the Z direction as a function of (r,phi,z)  (centered in r,phi, z);phi;r;z",nph,pih,pfh,nrh,rih,rfh,nzh,zih,zfh);
  TH3F* hIntDistortionR=new TH3F("hIntDistortionR","Integrated R Distortion from (r,phi,z) to z=0 (centered in r,phi, and z);phi;r;z",nph,pih,pfh,nrh,rih,rfh,nzh,zih,zfh);
  TH3F* hIntDistortionP=new TH3F("hIntDistortionP","Integrated R Distortion from (r,phi,z) to z=0 (centered in r,phi, and z);phi;r;z",nph,pih,pfh,nrh,rih,rfh,nzh,zih,zfh);
  TH3F* hIntDistortionZ=new TH3F("hIntDistortionZ","Integrated R Distortion from (r,phi,z) to z=0  (centered in r,phi, and z);phi;r;z",nph,pih,pfh,nrh,rih,rfh,nzh,zih,zfh);

  TVector3 pos((nrh/2+0.5)*step.Perp()+rih,0,(nzh/2+0.5)*step.Z()+zih);
    pos.SetPhi((nph/2+0.5)*step.Phi()+pih);
  int xi[3]={nrh/2,nph/2,nzh/2};
  const char axname[]="rpzrpz";
  int axn[]={nrh,nph,nzh,nrh,nph,nzh};
  float axval[]={(float)pos.Perp(),(float)pos.Phi(),(float)pos.Z(),(float)pos.Perp(),(float)pos.Phi(),(float)pos.Z()};
  float axbot[]={rih,pih,zih,rih,pih,zih};
  float axtop[]={rfh,pfh,zfh,rfh,pfh,zfh};
  TH2F* hIntDist[3][2];
  for (int ax=0;ax<3;ax++){
    //loop over which plane to work in
    for (int i=0;i<2;i++){
      //loop over which axis of the distortion to read
      hIntDist[ax][i]=new TH2F(Form("hIntDist%c_%c%c",axname[i],axname[ax+1],axname[ax+2]),
			      Form("%c component of distortion in the %c%c plane at %c=%2.3f;%c;%c",
				   axname[i],axname[ax+1],axname[ax+2],axname[ax],axval[ax],axname[ax+1],axname[ax+2]),
			      axn[ax+1],axbot[ax+1],axtop[ax+1],
			      axn[ax+2],axbot[ax+2],axtop[ax+2]);

    }
  }


  
  float deltar=step.Perp();//(rf-ri)/nr;
  float deltap=step.Phi();//(pf-pi)/np;
  float deltaz=step.Z();//(zf-zi)/nz;
  TVector3 inpart,outpart;
  TVector3 distort;
  int validToStep;
  int nSteps=100;

  float partR,partP,partZ;
  int ir,ip,iz;
  float distortR,distortP,distortZ;
  TTree *dTree=new TTree("dTree","Distortion per step z");
  dTree->Branch("r",&partR);
  dTree->Branch("p",&partP);
  dTree->Branch("z",&partZ);
  dTree->Branch("ir",&ir);
  dTree->Branch("ip",&ip);
  dTree->Branch("iz",&iz);
  dTree->Branch("dr",&distortR);
  dTree->Branch("drp",&distortP);
  dTree->Branch("dz",&distortZ);


  printf("generating distortion map with (%dx%dx%d) grid \n",nrh,nph,nzh);
  unsigned long long totalelements=nrh*nph*nzh;
  unsigned long long percent=totalelements/100;
  printf("total elements = %llu\n",totalelements);


  int el=0;

  /*
  for (iz=-1;iz<nz+1;iz++){
    partZ=(iz)*deltaz+zi;
    if (iz<0){
      inpart.SetZ(partZ+deltaz);
    } else if (iz==nz){
      inpart.SetZ(partZ-deltaz);
    } else {
      inpart.SetZ(partZ);
    }
    partZ+=0.5*deltaz;

    printf("iz=%d, zcoord=%2.2f, bin=%d\n",iz,partZ,  hIntDist[0][0]->GetYaxis()->FindBin(partZ));
  }
  assert(1==2);
  return;
  */
  //we want to loop over the entire region to be mapped, but we also need to include
  //one additional bin at each edge, to allow the mc drift code to interpolate properly.
  //hence we count from -1 to n+1, and manually adjust the position in those edge cases
  //to avoid sampling nonphysical regions in r and z.  the phi case is free to wrap as
  // normal.

  //note that we adjust the particle position (inpart) and not the plotted position (partR etc)
  inpart.SetXYZ(1,0,0);
  for (ir=-1;ir<nr+1;ir++){
    partR=(ir+0.5)*deltar+ri;
    if (ir<0){
      inpart.SetPerp(partR+deltar);
    } else if (ir==nr){
      inpart.SetPerp(partR-deltar);
    } else {
      inpart.SetPerp(partR);
    }
    for (ip=-1;ip<np+1;ip++){
      partP=(ip+0.5)*deltap+pi;
      inpart.SetPhi(partP);
      //if (partP<0) partP+=TMath::TwoPi();
      for (iz=-1;iz<nz+1;iz++){
	partZ=(iz)*deltaz+zi;
	if (iz<0){
	  inpart.SetZ(partZ+deltaz);
	} else if (iz==nz){
	  inpart.SetZ(partZ-deltaz);
	} else {
	  inpart.SetZ(partZ);
	}
	partZ+=0.5*deltaz;

	//printf("iz=%d, zcoord=%2.2f, bin=%d\n",iz,partZ,  hIntDist[0][0]->GetYaxis()->FindBin(partZ));

	//differential distortion:
	//be careful with the math of a distortion.  The R distortion is NOT the perp() component of outpart-inpart -- that's the transverse magnitude of the distortion!
	outpart=t->swimToInSteps(inpart.Z()+deltaz,inpart,nSteps,true, &validToStep);
	distort=outpart-inpart;
	distort.RotateZ(-inpart.Phi());//rotate so that that is on the x axis
	distortP=distort.Y();//the phi component is now the y component.
	distortR=distort.X();//and the r component is the z component
	distortZ=0;
	// also works, but not as good:
	//distortR=outpart.Perp()-inpart.Perp();
	//distortP=outpart.Perp()*(outpart.Phi()-inpart.Phi()); //this will have a seam.
	//distortZ=0;
	
	hDistortionR->Fill(partP,partR,partZ,distortR);
	hDistortionP->Fill(partP,partR,partZ,distortP);
	hDistortionZ->Fill(partP,partR,partZ,distortZ);
	dTree->Fill();

	//integral distortion:
	outpart=t->swimToInSteps(zf,inpart,nSteps,true, &validToStep);
	distort=outpart-inpart;
	distort.RotateZ(-inpart.Phi());//rotate so that that is on the x axis
	distortP=distort.Y();//the phi component is now the y component.
	distortR=distort.X();//and the r component is the z component
	distortZ=0;

	//printf("part=(rpz)(%f,%f,%f),distortP=%f\n",partP,partR,partZ,distortP);
	hIntDistortionR->Fill(partP,partR,partZ,distortR);
	hIntDistortionP->Fill(partP,partR,partZ,distortP);
	hIntDistortionZ->Fill(partP,partR,partZ,distortZ);

	//now we fill particular slices for visualizations:
	if(ir==xi[0]){//r slice
	  //printf("ir=%d, r=%f (pz)=(%d,%d), distortR=%2.2f, distortP=%2.2f\n",ir,partR,ip,iz,distortR,distortP);
	    hIntDist[0][0]->Fill(partP,partZ,distortR);
	    hIntDist[0][1]->Fill(partP,partZ,distortP);
	}
	if(ip==xi[1]){//phi slice
	  //printf("ip=%d, p=%f (rz)=(%d,%d), distortR=%2.2f, distortP=%2.2f\n",ip,partP,ir,iz,distortR,distortP);

	    hIntDist[1][0]->Fill(partZ,partR,distortR);
	    hIntDist[1][1]->Fill(partZ,partR,distortP);
	}
	if(iz==xi[2]){//z slice
	  //printf("iz=%d, z=%f (rp)=(%d,%d), distortR=%2.2f, distortP=%2.2f\n",iz,partZ,ir,ip,distortR,distortP);
		  
	    hIntDist[2][0]->Fill(partR,partP,distortR);
	    hIntDist[2][1]->Fill(partR,partP,distortP);
	}

	if(!(el%percent)) {printf("generating distortions %d%%:  ",(int)(el/percent));
	  printf("distortion at (ir=%d,ip=%d,iz=%d) is (%E,%E,%E)\n",
		 ir,ip,iz,distortR,distortP,distortZ);
	}
	el++;
	
      }
    }
  }


  TCanvas *c=new TCanvas("cdistort","distortion integrals",1200,800);
  const char axes[]="XYZXYZ";

  c->Divide(3,2);
  gStyle->SetOptStat();
  // TH2F *hProjection[12];
  TAxis * axis;
  int ih=0;


  for (int ax=0;ax<3;ax++){
    //plane
    for (int i=0;i<2;i++){
      //component
      c->cd(i*3+ax+1);
      gPad->SetRightMargin(0.15);
      hIntDist[ax][i]->SetStats(0);
      hIntDist[ax][i]->Draw("colz");
    }
  }

  /*
  for (int ax=0;ax<3;ax++){
    if (ax==0) axis=hIntDistortionR->GetYaxis();
    if (ax==1) axis=hIntDistortionR->GetXaxis();
    if (ax==2) axis=hIntDistortionR->GetZaxis();
    axis->SetRange(xi[ax],xi[ax]+1);
    c->cd(ax*4+1);
    gPad->SetRightMargin(0.15);
    hProjection[ih]=(TH2F*)hIntDistortionR->Project3D(Form("%c%c%d",axes[ax+1],axes[ax+2],ih));
    hProjection[ih]->SetStats(0);
    hProjection[ih]->Draw("colz");
    ih++;
    gPad->Modified();

      c->cd(ax*4+2);
      gPad->SetRightMargin(0.15);
     hProjection[ih]=(TH2F*)hIntDistortionP->Project3D(Form("%c%c%d",axes[ax+1],axes[ax+2],ih));
      hProjection[ih]->SetStats(0);
      hProjection[ih]->Draw("colz");
      ih++;
      gPad->Modified();
      axis->SetRange(0,0);

  }
  for (int ax=0;ax<3;ax++){
    if (ax==0) axis=hIntDistortionR->GetYaxis();
    if (ax==1) axis=hIntDistortionR->GetXaxis();
    if (ax==2) axis=hIntDistortionR->GetZaxis();

    c->cd(ax*4+3);
    gPad->SetRightMargin(0.15);
    hProjection[ih]=(TH2F*)hDistortionR->Project3D(Form("%c%c%d",axes[ax+1],axes[ax+2],ih));
    hProjection[ih]->SetStats(0);
    hProjection[ih]->Draw("colz");
    ih++; 
    gPad->Modified();
    c->cd(ax*4+4);
    gPad->SetRightMargin(0.15);
    hProjection[ih]=(TH2F*)hDistortionP->Project3D(Form("%c%c%d",axes[ax+1],axes[ax+2],ih));
    hProjection[ih]->SetStats(0);
    hProjection[ih]->Draw("colz");
    ih++; 
    gPad->Modified();
    axis->SetRange(0,0);
	  
  }
  */
  c->SaveAs(TString::Format("%s.distortion_summary.pdf",filebase));
  hDistortionR->Write();
  hDistortionP->Write();
  hDistortionZ->Write();
  hIntDistortionR->Write();
  hIntDistortionP->Write();
  hIntDistortionZ->Write();
  dTree->Write();
  //outf->Close();
  printf("wrote map and summary to %s.\n",filebase);
  return;
}

void SaveChargeAndProjections(const char* filename,AnnularFieldSim *t,int nphi,float pi,float pf,int nr,float ri,float rf,int nz,float zi,float zf){
  //scan over the tpc physical volume in nphi steps from pi to pf, and similar for the other two dimensions,
  //extract the charge from the tpc object at those coordinates, then save the resulting histogram and projections.
  TH3F* hAnCharge=new TH3F("hAnCharge","hAnCharge;phi;r;z",nphi,pi,pf,nr,ri,rf,nz,zi,zf);
  hAnCharge->Fill(0.0,0.0,0.0,0.0); //prime it so it draws correctly
  TH1F *hproj[3];
  int rh,ph,zh;
  for (int i=0;i<hAnCharge->GetNcells();i++){
    hAnCharge->GetBinXYZ(i,ph,rh,zh);
    hAnCharge->SetBinContent(i,t->q->Get(rh%nr,ph%nphi,zh%nz));
  }
  hproj[0]=(TH1F*)(hAnCharge->Project3D("XZ"));
  hproj[1]=(TH1F*)(hAnCharge->Project3D("YZ"));
  hproj[2]=(TH1F*)(hAnCharge->Project3D("YX"));

  TFile *outf=TFile::Open(filename,"RECREATE");
  outf->cd();
  hAnCharge->Write();
  for (int i=0;i<3;i++){
    hproj[i]->Write();
  }
  outf->Close();
  return;
}


void SaveField(const char* filename,AnnularFieldSim *t,int pi,int pf,int ri, int rf, int zi, int zf){
  //scan over the tpc object's field map, from the cell center of (ri,pi,zi) to the cell center of (rf-1,pf-1,zf-1), inclusive,
  //and save the field vectors at each point.
  //note that this may misbehave if you go outside the roi.
  TFile *outf=TFile::Open(filename,"RECREATE");
  outf->cd();
 
  //save data about the Efield:
  TTree fTree("fTree","field Tree");
  TVector3 pos0,pos1,pos2,pos,Efield,Bfield,phihat;
  TVector3 zero(0,0,0);
  TVector3 Eint,EintA;
  bool inroi;
  float charge,eintp, ep;
  fTree.Branch("pos","TVector3",&pos);
  fTree.Branch("phihat","TVector3",&phihat);

  fTree.Branch("E","TVector3",&Efield);
  fTree.Branch("B","TVector3",&Bfield);
  fTree.Branch("Eint","TVector3",&Eint);
  fTree.Branch("EintA","TVector3",&EintA);
  //fTree.Branch("q",&charge);
  fTree.Branch("Eintp",&eintp);
  fTree.Branch("Ep",&ep);
  fTree.Branch("roi",&inroi);
  float delr,delp;//=tpc->GetCellCenter(2,0,0)-tpc->GetCellCenter(1,0,0);
  //TVector3 delp;
  float delz=t->GetCellCenter(0,0,1).Z()-t->GetCellCenter(0,0,0).Z();
  
  bool inr,inp,inz;
  int rl,pl, zl;
  int fieldmap_output_extra_sampling=0;
  if (0)
  for (int ir=ri;ir<rf;ir++){
    for (int ip=pi;ip<pf;ip++){
      for (int iz=zi;iz<zf;iz++){
	pos0=t->GetRoiCellCenter(ir,ip,iz);
	pos1=t->GetRoiCellCenter(ir+1,ip,iz);//shift by one point in r.
	pos2=t->GetRoiCellCenter(ir,ip+1,iz);//shift by one point in phi.
	delr=pos1.Perp()-pos0.Perp();
	delp=pos2.Phi()-pos0.Phi();
	if (pos0.Mag()==0) continue; // we aren't in the roi.
	if (pos1.Mag()==0) continue; // we couldn't get an accurate r spacing
	if (pos2.Mag()==0) continue; // we couldn't get an accurate phi spacing
	Efield=zero;
	Bfield=zero;
	Efield=t->Efield->Get(ir,ip,iz);
	Bfield=t->Bfield->Get(ir,ip,iz);
	int rlocal=0;
	for (int plocal=-fieldmap_output_extra_sampling;plocal<fieldmap_output_extra_sampling+1;plocal++){
	  pos=pos0;
	  pos.SetPerp(pos0.Perp()+delr/(1+2*fieldmap_output_extra_sampling)*rlocal);
	  pos.RotateZ(delp/(1+2*fieldmap_output_extra_sampling)*plocal);//+(rlocal)/20.1*delr+(plocal)/(20.1)*delp;
	  //printf("trying pos=(%f,%f,%f)= rphiz(%f,%f,%f)\n",pos.X(),pos.Y(),pos.Z(),pos.Perp(),pos.Phi(),pos.Z());
	  phihat=pos;// build our phi position by starting with the vector in the rz plane:
	  phihat.SetZ(0);//remove the z component so it points purely in r
	  phihat.SetMag(1.0);//scale it to 1.0;
	  phihat.RotateZ(TMath::Pi()/2);//rotate 90 degrees from the position vector so it now points purely in phi; 
	  Eint=(t->interpolatedFieldIntegral(pos.Z()-delz/(4.0),pos))*(4.0/delz);
	  EintA=(t->analyticFieldIntegral(pos.Z()-delz/(4.0),pos))*(4.0/delz);
	  //charge=tpc->q->Get(ir,ip,iz);
	  eintp=Eint.Dot(phihat);
	  ep=Efield.Dot(phihat);
	  fTree.Fill();
	}
      }
    }
  }
  fTree.Write();
  outf->Close();
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

  PlotFieldSlices("testcharge",t,cpos);

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


void PlotFieldSlices(const char *filebase, AnnularFieldSim *t, TVector3 pos){

  
  printf("plotting field slices...\n");
  printf("file=%s\n",filebase);
  
  TVector3 inner=t->GetInnerEdge();
  TVector3 outer=t->GetOuterEdge();
  TVector3 step=t->GetFieldStep();
  

  int nr=t->GetFieldStepsR();
  int np=t->GetFieldStepsPhi();
  int nz=t->GetFieldStepsZ();

  TCanvas *c;
  TPaletteAxis *pal;

  TH2F *hEfield[3][3];
  TH2F *hCharge[3];
  TH1F *hEfieldComp[3][3];
  char axis[]="rpzrpz";
  float axisval[]={(float)pos.Perp(),(float)pos.Phi(),(float)pos.Z(),(float)pos.Perp(),(float)pos.Phi(),(float)pos.Z()};
  int axn[]={nr,np,nz,nr,np,nz};
  float axtop[]={(float)outer.Perp(),TMath::TwoPi(),(float)outer.Z(),(float)outer.Perp(),TMath::TwoPi(),(float)outer.Z()};
  float axbot[]={(float)inner.Perp(),0,(float)inner.Z(),(float)inner.Perp(),0,(float)inner.Z()};
  float axstep[6];
  for (int i=0;i<6;i++){
    axstep[i]=(axtop[i]-axbot[i])/(1.0*axn[i]);
  }
  TVector3 field;
  TVector3 lpos;
  for (int ax=0;ax<3;ax++){
    //loop over which axis slice to take
    hCharge[ax]=new TH2F(Form("hCharge%c",axis[ax]),
			 Form("Spacecharge Distribution in the %c%c plane at %c=%2.3f;%c;%c",
			      axis[ax+1],axis[ax+2],axis[ax],axisval[ax],axis[ax+1],axis[ax+2]),
			 axn[ax+1],axbot[ax+1],axtop[ax+1],
			 axn[ax+2],axbot[ax+2],axtop[ax+2]);
    for (int i=0;i<3;i++){
      //loop over which axis of the field to read
      /* old version that centers the plotting region by extending half-steps beyond the sampling region:
      hEfield[ax][i]=new TH2F(Form("hEfield%c_%c%c",axis[i],axis[ax+1],axis[ax+2]),
			      Form("%c component of E Field in the %c%c plane at %c=%2.3f;%c;%c",
				   axis[i],axis[ax+1],axis[ax+2],axis[ax],axisval[ax],axis[ax+1],axis[ax+2]),
			      axn[ax+1],axbot[ax+1]-axstep[ax+1]/2,axtop[ax+1]-axstep[ax+1]/2,
			      axn[ax+2],axbot[ax+2]-axstep[ax+2]/2,axtop[ax+2]-axstep[ax+2]/2);
      */
      //new version that centers the sampling region by starting half-steps into the plotting region:
      hEfield[ax][i]=new TH2F(Form("hEfield%c_%c%c",axis[i],axis[ax+1],axis[ax+2]),
			      Form("%c component of E Field in the %c%c plane at %c=%2.3f;%c;%c",
				   axis[i],axis[ax+1],axis[ax+2],axis[ax],axisval[ax],axis[ax+1],axis[ax+2]),
			      axn[ax+1],axbot[ax+1],axtop[ax+1],
			      axn[ax+2],axbot[ax+2],axtop[ax+2]);
      hEfieldComp[ax][i]=new TH1F(Form("hEfieldComp%c_%c%c",axis[i],axis[ax+1],axis[ax+2]),
			      Form("Log Magnitude of %c component of E Field in the %c%c plane at %c=%2.3f;log10(mag)",
				   axis[i],axis[ax+1],axis[ax+2],axis[ax],axisval[ax]),
				  200,-5,5);
    }
  }

  float rpz_coord[3];
  for (int ax=0;ax<3;ax++){
    rpz_coord[ax]=axisval[ax]+axstep[ax]/2;
    for (int i=0;i<axn[ax+1];i++){
      rpz_coord[(ax+1)%3]=axbot[ax+1]+(i+0.5)*axstep[ax+1];
      for (int j=0;j<axn[ax+2];j++){
	rpz_coord[(ax+2)%3]=axbot[ax+2]+(j+0.5)*axstep[ax+2];
	lpos.SetXYZ(rpz_coord[0],0,rpz_coord[2]);
	lpos.SetPhi(rpz_coord[1]);
	if (0 && ax==0){
	  printf("sampling rpz=(%f,%f,%f)=(%f,%f,%f) after conversion to xyz=(%f,%f,%f)\n",
		 rpz_coord[0],rpz_coord[1],rpz_coord[2],
		 lpos.Perp(),lpos.Phi(),lpos.Z(),lpos.X(),lpos.Y(),lpos.Z());
	}
	field=t->GetFieldAt(lpos);
	field.RotateZ(-rpz_coord[1]);//rotate us so we can read the y component as the phi component
	//if (field.Mag()>0) continue; //nothing has mag zero because of the drift field.
	hEfield[ax][0]->Fill(rpz_coord[(ax+1)%3],rpz_coord[(ax+2)%3],field.X());
	hEfield[ax][1]->Fill(rpz_coord[(ax+1)%3],rpz_coord[(ax+2)%3],field.Y());
	hEfield[ax][2]->Fill(rpz_coord[(ax+1)%3],rpz_coord[(ax+2)%3],field.Z());
	hCharge[ax]->Fill(rpz_coord[(ax+1)%3],rpz_coord[(ax+2)%3],t->GetChargeAt(lpos));
	hEfieldComp[ax][0]->Fill((abs(field.X())));
	hEfieldComp[ax][1]->Fill((abs(field.Y())));
	hEfieldComp[ax][2]->Fill((abs(field.Z())));
      }
    }
  }

  c=new TCanvas("cfieldslices","electric field",1200,800);
  c->Divide(4,3);
  gStyle->SetOptStat();
  for (int ax=0;ax<3;ax++){
    for (int i=0;i<3;i++){
      c->cd(ax*4+i+1);
      gPad->SetRightMargin(0.15);
      hEfield[ax][i]->SetStats(0);
      hEfield[ax][i]->Draw("colz");
      hEfield[ax][i]->GetListOfFunctions()->Print();
      	gPad->Modified();

      pal=dynamic_cast<TPaletteAxis*>(hEfield[ax][i]->GetListOfFunctions()->FindObject("palette"));
      if (pal){
	pal->SetX1NDC(0.86);
	pal->SetX2NDC(0.91);
	gPad->Modified();
      }
      //hEfieldComp[ax][i]->Draw();//"colz");
    }
    c->cd(ax*4+4);
    gPad->SetRightMargin(0.15);
    hCharge[ax]->SetStats(0);
    hCharge[ax]->Draw("colz");
    pal=dynamic_cast<TPaletteAxis*>(hCharge[ax]->GetListOfFunctions()->FindObject("palette"));
      if (pal){
	pal->SetX1NDC(0.86);
	pal->SetX2NDC(0.91);
	gPad->Modified();
      }
  }
  c->SaveAs(TString::Format("%s.field_slices.pdf",filebase));
  printf("after plotting field slices...\n");
  printf("file=%s\n",filebase);
  
return;
}
