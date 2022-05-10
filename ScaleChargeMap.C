//this is code to scale an input charge map so we can rapidly test different gas scenarios.

void ScaleChargeMap(){
  //for simplicity, I am hardcoding this first pass:

  TFile *infile,*outfile;
  infile=TFile::Open("newAverage.Apr.2022.hist.root","READ");
  outfile=TFile::Open("scaledAverage.Apr.2022.hist.root","RECREATE");
  TH3* original_primaries=(TH3*)infile->Get("_h_SC_prim_0");
  TH3* original_ibf=(TH3*)infile->Get("_h_SC_ibf_0");

  //pertinent gas params:
  int nGas=4;
  float nominal_total_gain=2000;//used for clarity of eqns, but actual value shouldn't matter so long as it's the same across all samples
  float nominal_ibf_fraction=0.0013;//used for clarity, see above.
  TString gasName[]={"Ne:CF4 50:50 (400V/cm)","Ar:CF4 60:40 (350V/cm)","Ar:CF4 60:40 (450V/cm)","Ar:CF4 60:40 (400V/cm)};
  float ion_mobility[]={1.34,1.34,1.34,0};//[cm^2/(V*s)]
  float drift_field[]={400,349.9,451.1,0};//{V/cm]
  float gas_ions_per_cm[]={71.5,96.4,96.4,0};//[#]

  ion_mobility[3]=0.5*(ion_mobility[1]+ion_mobility[2]);
  drift_field[3]=0.5*(drift_field[1]+drift_field[2]);
  gas_ions_per_cm[3]=0.5*(gas_ions_per_cm[1]+gas_ions_per_cm[2]);
  
  //derived quantities:
  float ion_velocity[nGas]; //ion_mobility*drift_field;
  float gem_gain[nGas];//set such that the electrons (ions) per cm and the electron gain yields a fixed value of 2e3 ('nominal_total_gain')
  float ibf_gain[nGas];//number of ions per input electron
  float velocity_ratio[nGas];//ratio of this velocity to the reference velocity [0];
  float ibf_gain_ratio[nGas];//ratio of this gain to the reference gain [0];
  float primary_ratio[nGas];//ratio of the primaries per cm of this setting to the reference primaries per cm, which ought to be the same as the ionizations per keV?
  float ibf_scale[nGas];//factor to multiply the ibf histogram by, including ion_velocity and ibf_gain:  Goes up linearly with the ibf_ratio (more charge if more ibf gain), up linearly with the primary_ratio (more charge if more primaries arrive), and down linearly with the velocity ratio (more charge if lower velocity)
  //NOTE:  since nominal_gain=gem_gain*gas_ions_per_cm is held fixed, ibf_gain_ratio*primary_ratio=(gem_gain[i]/gem_gain[0])*(gas_ions_per_cm[i]/gas_ions_per_cm[0])=nominal_gain[i]/nominal_gain[0] is also fixed.  That means ibf_scale is just the ion_velocity scaling.
  float primary_scale[nGas]; ////factor to multiply the primaries histogram by, including ion_velocity and ions_per_cm:  Goes up linearly with the primary_ratio (more charge if more primaries) and down linearly with the velocity ratio (more charge if lower velocity)

  for (int i=0;i<nGas;i++){
    ion_velocity[i]=ion_mobility[i]*drift_field[i];
    gem_gain[i]=nominal_total_gain/gas_ions_per_cm[i];
    ibf_gain[i]=gem_gain[i]*nominal_ibf_fraction;
    velocity_ratio[i]=ion_velocity[i]/ion_velocity[0];
    ibf_gain_ratio[i]=ibf_gain[i]/ibf_gain[0];
    primary_ratio[i]=gas_ions_per_cm[i]/gas_ions_per_cm[0];
    
    primary_scale[i]=primary_ratio[i]/velocity_ratio[i];
    //most of this cancels out, but it doesn't take long to compute:
    ibf_scale[i]=ibf_gain_ratio[i]*primary_ratio[i]/velocity_ratio[i];
    //the shortcut:
    // ibf_scale[i]=1./velocity_ratio[i];
  }
   
  outfile->cd();
  TH3* chargemap[nGas];
  TH3 *primarycharge[nGas];
  TH3 *ibfcharge[nGas];
  for (int i=0;i<nGas;i++){
    primarycharge[i]=(TH3*)original_primaries->Clone(Form("hprim%d",i));
    primarycharge[i]->SetTitle(Form("Primary Spacecharge with %s",gasName[i].Data()));
    primarycharge[i]->Scale(primary_scale[i]);
    primarycharge[i]->Write();
    
    ibfcharge[i]=(TH3*)original_ibf->Clone(Form("hibf%d",i));
    ibfcharge[i]->SetTitle(Form("IBF Spacecharge with %s",gasName[i].Data()));
    ibfcharge[i]->Scale(ibf_scale[i]);
    ibfcharge[i]->Write();

    chargemap[i]=(TH3*)primarycharge[i]->Clone(Form("htotalcharge%d",i));
    ibfcharge[i]->SetTitle(Form("Total Spacecharge with %s",gasName[i].Data()));
    chargemap[i]->Add(ibfcharge[i]);
    chargemap[i]->Write();

  }
  outfile->Close();
  return;
}
  
 			   
