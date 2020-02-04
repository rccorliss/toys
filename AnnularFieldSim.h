#include "assert.h"
#include "TVector3.h"
#include "AnalyticFieldModel.h"


template <class T> class MultiArray;
class TH3F;

class AnnularFieldSim{
 public:
  enum BoundsCase {InBounds,OnHighEdge, OnLowEdge,OutOfBounds}; //note that 'OnLowEdge' is qualitatively different from 'OnHighEdge'.  Low means there is a non-zero distance between the point and the edge of the bin.  High applies even if that distance is exactly zero.
  enum LookupCase {Full3D,HybridRes, PhiSlice, Analytic};
  //Full3D = uses (nr x nphi x nz)^2 lookup table
  //Hybrid = uses (nr x nphi x nz) x (nr_local x nphi_local x nz_local) + (nr_low x nphi_low x nz_low)^2 set of tables
  //PhiSlice = uses (nr x 1 x nz) x (nr x nphi x nz) lookup table exploiting phi symmetry.
  //Analytic = doesn't use lookup tables -- no memory footprint, uses analytic E field at center of each bin.
  //    Note that this is not the same as analytic propagation, which checks the analytic field integrals in each step.

  //debug items
  //
  int debug_printActionEveryN;
  int debug_printCounter;
  AnalyticFieldModel *aliceModel;
  
  //constants of motion, dimensions, etc:
  //
  //static constexpr float k=8.987e13;//=1/(4*pi*eps0) in N*cm^2/C^2 in a vacuum. N*cm^2/C units, so that we supply space charge in coulomb units.
  static constexpr float k_perm=8.987e11;//=1/(4*pi*eps0) in (V*cm)/C in a vacuum. so that we supply space charge in Coulombs, distance in cm, and fields in V/cm
  double vdrift; //gas drift speed in cm/s
  //double vprime; //first derivative of drift velocity at specific E
  //double vprime2; //second derivative of drift velocity at specific E
  float Bscale;//additional scale factor for debugging B effects. Defaults to 1.0
  float Escale;//additional scale factor for debugging E effects. Defaults to 1.0
  float Enominal;//magnitude of the nominal field on which drift speed is based, in V/cm.
  float phispan;//angular span of the area in the phi direction, since TVector3 is too smart.
  float rmin, rmax;//inner and outer radii of the annulus
  float zmin, zmax;//lower and upper edges of the coordinate system in z (not fully implemented yet)
  //float phimin, phimax;//not implemented at all yet.
  TVector3 dim;//dimensions of simulated region, in cm


  //variables related to the whole-volume tiling:
  //
  int nr,nphi,nz; //number of fundamental bins (f-bins) in each direction = dimensions of 3D array covering entire volume
  TVector3 step; //size of an f-bin in each direction
  LookupCase lookupCase; //which lookup system to instantiate and use.
  
  //variables related to the region of interest:
  //
  int rmin_roi, phimin_roi, zmin_roi; //lower edge of our region of interest, measured in f-bins
  int rmax_roi, phimax_roi, zmax_roi; //excluded upper edge of our region of interest, measured in f-bins
  int nr_roi,nphi_roi, nz_roi; //dimensions of our roi in f-bins

  
  //variables related to the high-res behavior:
  //
  int nr_high, nphi_high, nz_high; //dimensions, in f-bins of neighborhood of a f-bin in which we calculate the field in full resolution

  
  //variables related to the low-res behavior:
  //
  int r_spacing, phi_spacing, z_spacing; //number of f-bins, in each direction, to gang together to make a single low-resolution bin (l-bin)
  int nr_low, nphi_low, nz_low; //dimensions, in l-bins, of the entire volume
  int rmin_roi_low, phimin_roi_low, zmin_roi_low; //lowest l-bin that is at least partly in our region of interest
  int rmax_roi_low, phimax_roi_low, zmax_roi_low; //excluded upper edge l-bin of our region of interest
  int nr_roi_low,nphi_roi_low, nz_roi_low; //dimensions of our roi in l-bins

  
  //3- and 6-dimensional arrays to handle bin and bin-to-bin data
  //
  MultiArray<TVector3> *Efield; //total electric field in each f-bin in the roi for given configuration of charge AND external field.
  MultiArray<TVector3> *Epartial_highres; //electric field in each f-bin in the roi from charge in a given f-bin or summed bin in the high res region.
  MultiArray<TVector3> *Epartial_lowres; //electric field in each l-bin in the roi from charge in a given l-bin anywhere in the volume.
  MultiArray<TVector3> *Epartial; //electric field for the old brute-force model.
  MultiArray<TVector3> *Epartial_phislice; //electric field in a 2D phi-slice from the full 3D region.
  MultiArray<TVector3> *Eexternal; //externally applied electric field in each f-bin in the roi
  MultiArray<TVector3> *Bfield; //magnetic field in each f-bin in the roi
  MultiArray<double> *q; //space charge in each f-bin in the whole volume
  MultiArray<double> *q_local; //temporary holder of space charge in each f-bin and summed bin of the high-res region.
  MultiArray<double> *q_lowres; //space charge in each l-bin. = sums over sets of f-bins.

  
  


 public:
  AnnularFieldSim(float rmin,float rmax, float dz,int r,int phi, int z, float vdr); //abbr. constructor with roi=full region
  AnnularFieldSim(float rin, float rout, float dz,
		  int r, int roi_r0, int roi_r1,
		  int phi, int roi_phi0, int roi_phi1,
		  int z, int roi_z0, int roi_z1,
		  float vdr, LookupCase in_lookupCase=PhiSlice);
  AnnularFieldSim(float in_innerRadius, float in_outerRadius, float in_outerZ,
		  int r, int roi_r0, int roi_r1, int in_rLowSpacing, int in_rHighSize,
		  int phi, int roi_phi0, int roi_phi1, int in_phiLowSpacing, int in_phiHighSize,
		  int z, int roi_z0, int roi_z1,int in_zLowSpacing, int in_zHighSize,
		  float vdr, LookupCase in_lookupCase);

  //debug functions:
  void UpdateEveryN(int n){debug_printActionEveryN=n; return;};
  bool debugFlag(){
    if(debug_printActionEveryN>0 && debug_printCounter++>=debug_printActionEveryN){
      debug_printCounter=0;return true;
    } return false;};

  
  void load_spacecharge(TH3F *hist, float zoffset, float scalefactor);
  void load_analytic_spacecharge(float scalefactor);
  void setScaleFactorB(float x){Bscale=x;return;};
  void setScaleFactorE(float x){Escale=x;return;};
  void setFlatFields(float B, float E);

  TVector3 calc_unit_field(TVector3 at, TVector3 from);
  TVector3 analyticFieldIntegral(float zdest,TVector3 start);
  TVector3 interpolatedFieldIntegral(float zdest,TVector3 start);
  int FilterPhiIndex(int phi,int range); //defaults to using nphi for range.

  TVector3 GetCellCenter(int r, int phi, int z);
  TVector3 GetGroupCellCenter(int r0, int r1, int phi0, int phi1, int z0, int z1);
  TVector3 GetWeightedCellCenter(int r, int phi, int z);
  TVector3 fieldIntegral(float zdest,TVector3 start);
  void populate_fieldmap();
  //now handled by setting 'analytic' lookup:  void populate_analytic_fieldmap();
  void  populate_lookup();
  void  populate_full3d_lookup();
  void  populate_highres_lookup();
  void  populate_lowres_lookup();
  void  populate_phislice_lookup();
  TVector3 sum_field_at(int r,int phi, int z);
  TVector3 sum_full3d_field_at(int r,int phi, int z);
  TVector3 sum_local_field_at(int r,int phi, int z);
  TVector3 sum_nonlocal_field_at(int r,int phi, int z);
  TVector3 sum_phislice_field_at(int r, int phi, int z);
  TVector3 swimToInAnalyticSteps(float zdest,TVector3 start,int steps, int *goodToStep);
  TVector3 swimToInSteps(float zdest,TVector3 start, int steps, bool interpolate, int *goodToStep);
  TVector3 swimTo(float zdest,TVector3 start, bool interpolate=true, bool useAnalytic=false);
 
 private:
  BoundsCase GetRindexAndCheckBounds(float pos, int *r);
  BoundsCase GetPhiIndexAndCheckBounds(float pos, int *phi);
  BoundsCase GetZindexAndCheckBounds(float pos, int *z);

    float RosseggerEterm(int m, int n, TVector3 at, TVector3 from);

};

#ifndef MULTIARRAY
#define MULTIARRAY
template <class T>
class MultiArray : public TObject{
   //class to hold an up-to-six dimensional array of whatever T is.  Any indices not used are flattened.
 public:
   static const int MAX_DIM=6;
   int dim;
   int n[6];
   int length;
   T *field;

   MultiArray(int a=0, int b=0, int c=0, int d=0, int e=0, int f=0){
     int n_[6];
     for (int i=0;i<MAX_DIM;i++)
       n[i]=0;
     n_[0]=a; n_[1]=b; n_[2]=c; n_[3]=d; n_[4]=e; n_[5]=f;
     length=1;
     dim=MAX_DIM;
     for (int i=0;i<dim;i++){
       if (n_[i]<1) {
	 dim=i;
	 break;
       }
       n[i]=n_[i];
       length*=n[i];
     }
     field=static_cast<T*>( malloc(length*sizeof(T) ));
     //field=(T)( malloc(length*sizeof(T) ));
     //for (int i=0;i<length;i++) field[i].SetXYZ(0,0,0);
   }

   void Add(int a, int b, int c, T in){
     Set(a,b,c,0,0,0,in);
     return;
   };
   void Add(int a, int b, int c, int d, int e, int f, T in){
     int n_[6];
     n_[0]=a; n_[1]=b; n_[2]=c; n_[3]=d; n_[4]=e; n_[5]=f;
     int index=n_[0];
     for (int i=1;i<dim;i++){
       index=(index*n[i])+n_[i];
     }
     field[index]=field[index]+in;
     return; 
   }

   
   T Get(int a=0, int b=0, int c=0, int d=0, int e=0, int f=0){
     int n_[6];
     n_[0]=a; n_[1]=b; n_[2]=c; n_[3]=d; n_[4]=e; n_[5]=f;
     int index=0;
     for (int i=0;i<dim;i++){
       if (n[i]<=n_[i] || n_[i]<0){//check bounds
	 printf("asking for el %d %d %d %d %d %d.  %dth element is outside of bounds 0<x<%d\n"
		,n_[0],n_[1],n_[2],n_[3],n_[4],n_[5],n_[i],n[i]);
	 assert(false);
       }
       index=(index*n[i])+n_[i];
     }
     return field[index];
   }
   T* GetPtr(int a=0, int b=0, int c=0, int d=0, int e=0, int f=0){ //faster for repeated access.
     int n_[6];
     n_[0]=a; n_[1]=b; n_[2]=c; n_[3]=d; n_[4]=e; n_[5]=f;
     int index=n_[0];
     for (int i=1;i<dim;i++){
       index=(index*n[i])+n_[i];
     }
     return &(field[index]);
   }

   T* GetFlat(int a=0){
     if (a>=length) assert(false); //check bounds
     return &(field[a]);
   }

   int Length(){
     return length;
   }

   void Set(int a, int b, int c, T in){
     Set(a,b,c,0,0,0,in);
     return;
   };
   void Set(int a, int b, int c, int d, int e, int f, T in){
     int n_[6];
     n_[0]=a; n_[1]=b; n_[2]=c; n_[3]=d; n_[4]=e; n_[5]=f;
     int index=n_[0];
     for (int i=1;i<dim;i++){
       index=(index*n[i])+n_[i];
     }
     field[index]=in;
     return; 
   }
};
#endif //MULTIARRAY
