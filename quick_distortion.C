

#include "AnnularFieldSim.h"
#include "TTree.h" //this prevents a lazy binding issue and/or is a magic spell.
#include "TCanvas.h" //this prevents a lazy binding issue and/or is a magic spell.
R__LOAD_LIBRARY(.libs/libfieldsim)

  char field_string[200];
  char lookup_string[200];

AnnularFieldSim *SetupDefaultSphenixTpc(int nr,int np, int nz);


void quick_distortion(int nr, int np,int nz, const char * chargefile="hist.root", const char *outdir="./quick/"){


  //and some parameters of the files we're loading:
  bool usesChargeDensity=false; //true if source hists contain charge density per bin.  False if hists are charge per bin.
  float tpc_chargescale=1.6e-19;//Coulombs per bin unit.
  float spacecharge_cm_per_axis_unit=100;//cm per histogram axis unit.
  

  TFile *infile;

  TString sourcefilename=chargefile;
  TString basename;
  //get our basename by tokenizing:
  TString tok;
  Ssiz_t from = 0, lastslash=0;
  while (sourcefilename.Tokenize(tok, from, "/")) {
    basename=tok;
    lastslash=from;
  } 
  TString outputfilebase=Form("%s/%s",outdir,basename.Data());
  outputfilebase.ReplaceAll(".hist.root","");
  
  double totalQ=0;

  infile=TFile::Open(sourcefilename.Data(),"READ");
  TList *keys=infile->GetListOfKeys();
  //keys->Print();
  int nKeys=infile->GetNkeys();

  int j;
  bool isHist=false;
  TObject *tobj;
  for (j=0;j<nKeys;j++){
       tobj=infile->Get(keys->At(j)->GetName());
      //if this isn't a 3d histogram, skip it:
      isHist=tobj->InheritsFrom("TH3");
      if (isHist) break;
  }
  if (!isHist) assert(1=2);//no valid th3 in there!

  //don't bother building a tpc until we know we have a chargemap to feed it

  AnnularFieldSim *tpc=SetupDefaultSphenixTpc(nr,np,nz);//loads the lookup, fields, etc.
  //assume this histogram is a charge map.
  tpc->load_spacecharge(sourcefilename.Data(),tobj->GetName(),0,tpc_chargescale,spacecharge_cm_per_axis_unit, usesChargeDensity);
  
  printf("Sanity check:  Q has %d elements and dim=%d\n",tpc->q->Length(), tpc->q->dim);
  totalQ=0;
  for (int k=0;k<tpc->q->Length();k++){
    totalQ+=*(tpc->q->GetFlat(k));
  }
  printf("Sanity check:  Total Q in reported region is %E C\n",totalQ);
  
  tpc->populate_fieldmap();
  TString outputfilename=Form("%s.%s.%s.%s",outputfilebase.Data(),tobj->GetName(),field_string,lookup_string);
  printf("%s file has %s hist.  field=%s, lookup=%s. no scaling.\n",
	 sourcefilename.Data(),tobj->GetName(),field_string,lookup_string);
  tpc->GenerateDistortionMaps(outputfilename.Data(),2,2,2,1,true);
  printf("distortions mapped.\n");


  //add the location to plot the fieldslices about:
 TVector3 pos=0.5*(tpc->GetOuterEdge()+tpc->GetInnerEdge());;
  pos.SetPhi(3.14159);
  
  tpc->PlotFieldSlices(outputfilename.Data(),pos);
  tpc->PlotFieldSlices(outputfilename.Data(),pos,'B');
  printf("fieldslices plotted.\n");     
  printf("obj %d: getname: %s  inherits from TH3D:%d \n",j,tobj->GetName(),tobj->InheritsFrom("TH3"));

  infile->Close();
  return;
  
}



AnnularFieldSim *SetupDefaultSphenixTpc(int nr,int np, int nz){
  //step1:  specify the sPHENIX space charge model parameters
  const float tpc_rmin=20.0;
  const float tpc_rmax=78.0;
  float tpc_deltar=tpc_rmax-tpc_rmin;
  const float tpc_z=105.5;
  const float tpc_cmVolt=-400*tpc_z; //V =V_CM-V_RDO -- volts per cm times the length of the drift volume.
  //const float tpc_magField=0.5;//T -- The old value used in carlos's studies.
  //const float tpc_driftVel=4.0*1e6;//cm per s  -- used in carlos's studies
  const float tpc_driftVel=8.0*1e6;//cm per s  -- 2019 nominal value
  const float tpc_magField=1.5;//T -- 2019 nominal value
  const char detgeoname[]="sphenix";
  
   //step 2: specify the parameters of the field simulation.  Larger numbers of bins will rapidly increase the memory footprint and compute times.
  //there are some ways to mitigate this by setting a small region of interest, or a more parsimonious lookup strategy, specified when AnnularFieldSim() is actually constructed below.
  //nr=nr from the arguments
  int nr_roi_min=0;
  int nr_roi=nr;//10;
  int nr_roi_max=nr_roi_min+nr_roi;
  int nphi=np;
  int nphi_roi_min=0;
  int nphi_roi=nphi;
  int nphi_roi_max=nphi_roi_min+nphi_roi;
  //nz=nz from the arguments.
  int nz_roi_min=0;
  int nz_roi=nz;
  int nz_roi_max=nz_roi_min+nz_roi;

  

  //step 3:  create the fieldsim object.  different choices of the last few arguments will change how it builds the lookup table spatially, and how it loads the spacecharge.  The various start-up steps are exposed here so they can be timed in the macro.
  AnnularFieldSim *tpc=
    new  AnnularFieldSim(tpc_rmin,tpc_rmax,tpc_z,
			 nr, nr_roi_min,nr_roi_max,1,2,
			 nphi,nphi_roi_min, nphi_roi_max,1,2,
			 nz, nz_roi_min, nz_roi_max,1,2,
			 tpc_driftVel, AnnularFieldSim::PhiSlice, AnnularFieldSim::FromFile);
  tpc->UpdateEveryN(25);//show reports every 10%.

    //load the field maps, either flat or actual maps
  tpc->setFlatFields(tpc_magField,-tpc_cmVolt/tpc_z);
  sprintf(field_string,"flat_B%2.1f_E%2.1f",tpc_magField,tpc_cmVolt/tpc_z);

  if (1){
    tpc->loadBfield("sPHENIX.2d.root","fieldmap");
    tpc->loadEfield("externalEfield.ttree.root","fTree");
    sprintf(field_string,"real_B%2.1f_E%2.1f",tpc_magField,tpc_cmVolt/tpc_z);
  }
  printf("set fields.\n");




  //load the greens functions:
  sprintf(lookup_string,"ross_phi1_%s_phislice_lookup_r%dxp%dxz%d",detgeoname,nr,nphi,nz);
  char lookupFilename[200];
  sprintf(lookupFilename,"%s.root",lookup_string);
  TFile *fileptr=TFile::Open(lookupFilename,"READ");

  if (!fileptr){ //generate the lookuptable
  //to use the full rossegger terms instead of trivial free-space greens functions, uncomment the line below:
    tpc->load_rossegger();
    printf("loaded rossegger greens functions.\n");
    tpc->populate_lookup();
    tpc->save_phislice_lookup(lookupFilename);
  } else{ //load it from a file
    fileptr->Close();
    tpc->load_phislice_lookup(lookupFilename);
  }

    printf("populated lookup.\n");

  return tpc;
}
  
