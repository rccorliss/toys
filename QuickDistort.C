

class DistortionContainer
{
  public:
  // load all the things for dcc
  TH3F* m_hDPint[2] = {nullptr, nullptr};
  TH3F* m_hDRint[2] = {nullptr, nullptr};
  TH3F* m_hDZint[2] = {nullptr, nullptr};
  bool m_phi_hist_in_radians = false;
  int m_dimensions = 3;
  bool m_interpolate_z = false;
  bool m_use_scalefactor = false;
  float m_scalefactor = 1.0;

}
class TpcDistortionCorrectionMini
{
  public:
  DistortionContainer *dcc=nullptr;
  //! constructor
  TpcDistortionCorrectionMini() = default;

  //! destructor
  virtual ~TpcDistortionCorrectionMini() = default;

  //! apply distortion corrections to a cluster and return corrected position
  TVector3 get_corrected_position(const TVector3& source) const;
  void load_distortion_map(const std::string& filename);
}

//generate a vector of N triplets of x,y,z coordinates, corresponding to a straight line that passes through the origin with theta and phi given, and a max radius wrt the z axis of rmax
std::vector<TVector3> GenerateLine(float theta, float phi, float rmax, float rmin, int N)
{
  std::vector<TVector3> result;
  for (int i = 0; i < N; i++) {
    float r = rmin + (rmax - rmin) * (i + 0.5) / N;
    float x = r * sin(theta) * cos(phi);
    float y = r * sin(theta) * sin(phi);
    float z = r * cos(theta);
    result.push_back(TVector3(x, y, z));
  }
  return result;
}

//convert a vector of TVector3 to a triplet std::array<float, 3> 
std::array<float, 3> TVector3ToFloatArray(const TVector3& v)
{
  return {v.x(), v.y(), v.z()};
}

//and go the other way:
TVector3 FloatArrayToTVector3(const std::array<float, 3>& v)
{
  return TVector3(v[0], v[1], v[2]);
}


void QuickDistort(){

  //load the distortion map
  TpcDistortionCorrectionMini *tpcDistort=new TpcDistortionCorrectionMini();
  tpcDistort->load_distortion_map("elevatorpitch/fluct_single.1side.3d.file0.h_Charge_0.real_B1.4_E-400.0.ross_phi1_sphenix_phislice_lookup_r26xp40xz40.distortion_map.hist.root");


//create a set of M lines, each with N points, and distort them, then fit them to a straight line and histogram the distance of closest approach to the origin as a function of the theta parameter used to create that line.
  int M=100;
  int N=10;
  float rmax=75;//in cm
  float rmin=30;//in cm
  float theta=0;
  float dtheta=0.01;
  float phi=0;
  float dphi=0.01;
  TH1D *hDCA=new TH1D("hDCA","Distance of Closest Approach to Origin;DCA (cm);Entries",100,-10,10);
  for (int i=0;i<M;i++){
    theta+=dtheta;
    phi+=dphi;
    std::vector<TVector3> line=GenerateLine(theta,phi,rmax,rmin,N);
    std::vector<TVector3> distortedLine;
    for (int j=0;j<N;j++){
      distortedLine.push_back(tpcDistort->get_corrected_position(line[j]));
    }
    //fit the distorted line to a straight line
    TVector3 mean(0,0,0);
    for (int j=0;j<N;j++){
      mean+=distortedLine[j];
    }
    mean*=1.0/N;
    TMatrixD S(3,3);
    for (int j=0;j<N;j++){
      TVector3 diff=distortedLine[j]-mean;
      for (int k=0;k<3;k++){
          for (int l=0;l<3;l++){
            S[k][l]+=diff[k]*diff[l];
          }
      }
    }
    TVectorD eigenvalues;
    TMatrixD eigenvectors;






  //TFile *dFile=TFile::Open("elevatorpitch/fluct_average.rev3.1side.3d.file0.h_negz.real_B1.4_E-400.0.ross_phi1_sphenix_phislice_lookup_r26xp40xz40.distortion_map.hist.root","READ");
  TFile *dFile=TFile::Open("elevatorpitch/fluct_single.1side.3d.file0.h_Charge_0.real_B1.4_E-400.0.ross_phi1_sphenix_phislice_lookup_r26xp40xz40.distortion_map.hist.root","READ");
  TH3F *hDistortX=(TH3F*)dFile->Get("hIntDistortionX");
  TH3F *hDistortY=(TH3F*)dFile->Get("hIntDistortionY");
  //there's also a z component, but we'll ignore that for now.  It's small.
  
  //load the chargemap for comparison
  //TFile *cFile=TFile::Open("averages/average.rev3.hist.root","READ");
  TFile *cFile=TFile::Open("testfluctcharge/fluct_outputFile_15kHz_G4Hits_sHijing_0-12fm_000000_012000_bX1508071_bias0.root","READ");
  TH3D* hQ3D=(TH3D*)cFile->Get("h_Charge_0");
  TH2D* hQ2D=(TH2D*)hQ3D->Project3D("yx");
  

  //eventually we should probe at every stripe, but for now we'll go once per r-phi bin in the charge map:
  float phibound[2];
  float rbound[2];
  int nphi, nr;
  float phistep,rstep;
  phibound[0]=hQ2D->GetXaxis()->GetXmin();
  phibound[1]=hQ2D->GetXaxis()->GetXmax();
  nphi=hQ2D->GetXaxis()->GetNbins();
  phistep=(phibound[1]-phibound[0])/(nphi*1.0);
  rbound[0]=hQ2D->GetYaxis()->GetXmin();
  rbound[1]=hQ2D->GetYaxis()->GetXmax();
  nr=hQ2D->GetYaxis()->GetNbins();
  rstep=(rbound[1]-rbound[0])/(nr*1.0);
  TH2D* hCmDistortion=new TH2D("hCmDistortion","Radial Shift (r_f-r_i) of CM hits;phi (rad);r (m)",
			       nphi,phibound[0],phibound[1],
			       nr,rbound[0],rbound[1]);
  TH2D* hChargeLow=new TH2D("hChargeLow","Lower-Res Charge Fluctuations;phi (rad);r (m)",
			       nphi/5,phibound[0],phibound[1],
			       nr/5,rbound[0],rbound[1]);
  TH2D* hSanityCheck=new TH2D("hSanityCheck","Sanity Check of how many times each cell is touched;phi (rad);r (m)",
			       nphi,phibound[0],phibound[1],
			       nr,rbound[0],rbound[1]);
  TH2D* hRatio=new TH2D("hRatio","#Delta R / Qtot in the bin;phi (rad);r (m)",
			       nphi,phibound[0],phibound[1],
			       nr,rbound[0],rbound[1]);
  /* for sets that aren't average-subtracted:
  TH2D* hComparison=new TH2D("hComparison","#Delta R vs Qtot in the bin;#Delta R (cm);Qtot (ions/bin)",
			     200,-0.1,1.1,
			     200,0,2e7);
  */
  TH2D* hComparison=new TH2D("hComparison","#Delta R-ave vs Qtot-ave in the bin;#Delta R (cm);Qtot (ions/bin)",
			     200,-0.004,0.006,
			     200,-1e6,3e6);
  TH2D* hDistComp=new TH2D("hDistComp","Fluctuation #Delta R at cm and halfway;#Delta R@0 (um);#Delta R@50 (um)",
			     100,-50,50,
			     100,-50,50);
  vector<float> rdist,qint,rpos,phipos;//radial distortion and charge integral per bin.
  vector<float> rdisthalf;//radial distortion at z=50cm

  TVector3 posBefore, posAfter;
  for (int i=0;i<nr;i++){
    float r=rbound[0]+rstep*(i+0.5);
    posBefore.SetXYZ(r,0,0);//neglect z coordinate for now.
    for (int j=0;j<nphi;j++){
      float phi=phibound[0]+phistep*(j+0.5);
      posBefore.SetPhi(phi);
      //printf("checking phi=%f, r=%f\n",phi,r);
      float deltaX=hDistortX->Interpolate(phi,r*100.0, 1.0);//map is in cm, but charge is in m, so need to convert r here.  set z=1.0cm
      float deltaY=hDistortY->Interpolate(phi,r*100.0,1.0);//map is in cm, but charge is in m, so need to convert r here.
      posAfter.SetX(posBefore.X()+deltaX);
      posAfter.SetY(posBefore.Y()+deltaY);
      float deltaR=posAfter.Perp()-posBefore.Perp();
      float qIntegral=hQ2D->GetBinContent(hQ2D->FindBin(phi,r));

      deltaX=hDistortX->Interpolate(phi,r*100.0, 50.0);//map is in cm, but charge is in m, so need to convert r here.  set z=1.0cm
      deltaY=hDistortY->Interpolate(phi,r*100.0,50.0);//map is in cm, but charge is in m, so need to convert r here.
      posAfter.SetX(posBefore.X()+deltaX);
      posAfter.SetY(posBefore.Y()+deltaY);
      float deltaR2=posAfter.Perp()-posBefore.Perp();
      
      hCmDistortion->Fill(phi,r,deltaR);
      hChargeLow->Fill(phi,r,qIntegral);
      hSanityCheck->Fill(phi,r,1);
      //hRatio->Fill(phi,r,abs(deltaR)/qIntegral);
      hRatio->Fill(phi,r,deltaR/qIntegral);
      hComparison->Fill(deltaR,qIntegral);
      hDistComp->Fill(deltaR*1e4,deltaR2*1e4);
      rdist.push_back(deltaR);
      rdisthalf.push_back(deltaR2);
      qint.push_back(qIntegral);
      rpos.push_back(r);
      phipos.push_back(phi);
    }
  }

  TCanvas *c=new TCanvas("c","distortion comparison to charge",800,800);
  c->Divide(2,2);
  c->cd(1);
  hChargeLow->Draw("colz");
  c->cd(2);
  hCmDistortion->Draw("colz");
  c->cd(3);//->SetLogz();
  TGraph *g;
  g=new TGraph(rdist.size(),&(rdist[0]),&(rdisthalf[0]));
  g->SetTitle("Distortion at z=0 vs at z=50cm;dist@0;dist@50");
  //g->Draw("*A");
  hDistComp->Draw("colz");
  /*
  g=new TGraph(rdist.size(),&(rdist[0]),&(qint[0]));
  g->SetTitle("Distortion vs integrated charge in column;dist;q");
  g->Draw("*A");
  */
  //hRatio->Draw("colz");
  c->cd(4);
  //this looked fine, but if you mess with stuff, look again: hSanityCheck->Draw("colz");
  g=new TGraph(rdist.size(),&(rpos[0]),&(rdist[0]));
  //g->SetTitle("Distortion vs radial position;r (m);#Delta R");
  //g->Draw("*A");
  hComparison->Draw("colz");
  //hqluusinqqtac iie
   

  
  return;
}




namespace
{
  template <class T>
  inline constexpr T square(const T& x)
  {
    return x * x;
  }

  // check boundaries in axis
  /* for the interpolation to work, the value must be within the range of the provided axis, and not into the first and last bin */
  inline bool check_boundaries(const TAxis* axis, double value)
  {
    const auto bin = axis->FindBin(value);
    return (bin >= 2 && bin < axis->GetNbins());
  }

  // check boundaries in histogram, before interpolation
  /* for the interpolation to work, the value must be within the range of the provided axis, and not into the first and last bin */
  inline bool check_boundaries(const TH1* h, double r, double phi, double z)
  {
    return check_boundaries(h->GetXaxis(), r) && check_boundaries(h->GetYaxis(), phi) && check_boundaries(h->GetZaxis(), z);
  }

  // check boundaries in histogram, before interpolation
  /* for the interpolation to work, the value must be within the range of the provided axis, and not into the first and last bin */
  inline bool check_boundaries(const TH1* h, double r, double phi)
  {
    return check_boundaries(h->GetXaxis(), r) && check_boundaries(h->GetYaxis(), phi);
  }

}  // namespace

TpcDistortionCorrectionMini::load_distortion_map(string m_correction_filename)
{
  dcc = new DistortionContainer();
  dcc->m_phi_hist_in_radians = false;
 
    std::cout << "TpcLoadDistortionCorrection::InitRun - reading corrections from " << m_correction_filename << std::endl;
    auto distortion_tfile = TFile::Open(m_correction_filename.c_str());
    if (!distortion_tfile)
    {
      std::cout << "TpcLoadDistortionCorrection::InitRun - cannot open " << m_correction_filename << std::endl;
      exit(1);
    }


    const std::array<const std::string, 2> extension = {{"_negz", "_posz"}};
    for (int j = 0; j < 2; ++j)
    {
      dcc->m_hDPint[j] = dynamic_cast<TH1*>(distortion_tfile->Get((std::string("hIntDistortionP")+extension[j]).c_str()));
      assert(dcc->m_hDPint[j]);
      dcc->m_hDRint[j] = dynamic_cast<TH1*>(distortion_tfile->Get((std::string("hIntDistortionR")+extension[j]).c_str()));
      assert(dcc->m_hDRint[j]);
      dcc->m_hDZint[j] = dynamic_cast<TH1*>(distortion_tfile->Get((std::string("hIntDistortionZ")+extension[j]).c_str()));
      assert(dcc->m_hDZint[j]);
    }
}

//________________________________________________________
TVector3 TpcDistortionCorrectionMini::get_corrected_position(const TVector3& source) const
{
  // get cluster radius, phi and z
  const auto r = std::sqrt(square(source.X()) + square(source.Y()));
  auto phi = std::atan2(source.Y(), source.X());
  if (phi < 0)
  {
    phi += 2 * M_PI;
  }

  const auto z = source.Z();
  const int index = z > 0 ? 1 : 0;

  // apply corrections
  auto phi_new = phi;
  auto r_new = r;
  auto z_new = z;

  // if the phi correction hist units are cm, we must divide by r to get the dPhi in radians
  auto divisor = r;

  if (dcc->m_phi_hist_in_radians)
  {
    // if the phi correction hist units are radians, we must not divide by r.
    divisor = 1.0;
  }

  //set our default corrections to be zero:
  //first inherit the same type
  auto dphi=phi;
  auto dr=r;
  auto dz=z;
  //then set them to zero:
  dphi=0;
  dr=0;
  dz=0;
  
  //get the corrections from the histograms
  if (dcc->m_dimensions == 3)
  {
    if (dcc->m_hDPint[index]  && check_boundaries(dcc->m_hDPint[index], phi, r, z))
    {
      dphi=dcc->m_hDPint[index]->Interpolate(phi, r, z) / divisor;
    }
    if (dcc->m_hDRint[index] && check_boundaries(dcc->m_hDRint[index], phi, r, z))
    {
      dr=dcc->m_hDRint[index]->Interpolate(phi, r, z);
    }
    if (dcc->m_hDZint[index] && check_boundaries(dcc->m_hDZint[index], phi, r, z))
    {
      dz=dcc->m_hDZint[index]->Interpolate(phi, r, z);
    }
  }
  else if (dcc->m_dimensions == 2)
  {
    double zterm = 1.0;

    if (dcc->m_interpolate_z){
      zterm=(1. - std::abs(z) / 105.5);
    }
    if (dcc->m_hDPint[index]  && check_boundaries(dcc->m_hDPint[index], phi, r))
    {
      dphi=dcc->m_hDPint[index]->Interpolate(phi, r) * zterm / divisor;
    }
    if (dcc->m_hDRint[index] && check_boundaries(dcc->m_hDRint[index], phi, r))
    {
      dr=dcc->m_hDRint[index]->Interpolate(phi, r) * zterm;
    }
    if (dcc->m_hDZint[index] && check_boundaries(dcc->m_hDZint[index], phi, r))
    {
      dz=dcc->m_hDZint[index]->Interpolate(phi, r) * zterm;
    }
    
  }

//if we are scaling, apply the scale factor to each correction
if(dcc->m_use_scalefactor)
  {
    dphi *= dcc->m_scalefactor;
    dr *= dcc->m_scalefactor;
    dz *= dcc->m_scalefactor;
  }

  phi_new=phi-dphi;
  r_new=r-dr;
  z_new=z-dz;

  // update cluster
  const auto x_new = r_new * std::cos(phi_new);
  const auto y_new = r_new * std::sin(phi_new);

//convert x_new ,y_new, z_new to TVector3
TVector3 new_pos(x_new, y_new, z_new);

  return new_pos;
}