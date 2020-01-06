#include "assert.h"
#include "TVector3.h"

template <class T> class MultiArray;
class TH3F;
class AnnularFieldSim{
 public:
  enum BoundsCase {InBounds,OnHighEdge, OnLowEdge,OutOfBounds}; //note that 'OnLowEdge' is qualitatively different from 'OnHighEdge'.  Low means there is a non-zero distance between the point and the edge of the bin.  High applies even if that distance is exactly zero.


  //bad form to leave this all exposed, I know.

  int nr,nphi,nz; //dimensions of internal high-res charge grid
  TVector3 step; //step size in each direction. in the high-res grid.
  float phispan;//angular span of the area in the phi direction, since TVector3 is too smart.
  float rmin, rmax;//inner and outer radii of the annulus
  float zmin, zmax;//lower and upper edges of the coordinate system in z (not fully implemented yet)
  //float phimin, phimax;//not implemented at all yet.

  //allow us to track particles only in a subset of the region, to study sections in greater detail
  //particles will only track in indices rmin_roi<=r<rmax_roi etc.
  int rmin_roi, rmax_roi, nr_roi;
  int zmin_roi, zmax_roi, nz_roi;
  int phimin_roi, phimax_roi, nphi_roi;

  //allow us to handle local regions in high-res, and the rest in lower res:
  int r_spacing, phi_spacing, z_spacing;//number of cells in the high-res grid that are ganged together (save for the last bin, which might have fewer)
  int nr_low, nphi_low, nz_low;//dimensions of internal low-res grids.
  int nr_high, nphi_high, nz_high;//dimensions of internal hig-res local grids.
  MultiArray<TVector3> *Epartial_lowres;
  MultiArray<TVector3> *Epartial_highres;
   MultiArray<float> *q_local;

  
  double vdrift; //gas drift speed.
  //double omegatau; //gas propagation constant depends on field and vdrift
  float Bscale;
  float Escale;

  
  static constexpr float k=8.987e13;//gas electric permeability N*cm^2/C^2 in a vacuum.
   TVector3 dim;//dimensions of simulated region, in cm
  
  MultiArray<TVector3> *Efield; //electric field for given configuration of charge AND external field.
  MultiArray<TVector3> *Epartial; //electric field for unit charge in given cell.
  MultiArray<TVector3> *Eexternal; //externally applied electric field
  MultiArray<TVector3> *Bfield; //magnetic field for system.
  MultiArray<float> *q; //space charge

  MultiArray<int> *rUpperBounds; //the (excluded) upper bounds of the effective rbins for each r point in the high-res grid
  MultiArray<int> *phiUpperBounds; //the (excluded) upper bounds of the effective phibins for each phi point in the high-res grid
  MultiArray<int> *zUpperBounds; //the (excluded) upper bounds of the effective zbins for each z point in the high-res grid
  
 public:
  AnnularFieldSim(float rmin,float rmax, float dz,int r,int phi, int z, float vdr); //abbr. constructor with roi=full region
  AnnularFieldSim(float rin, float rout, float dz,
		  int r, int roi_r0, int roi_r1,
		  int phi, int roi_phi0, int roi_phi1,
		  int z, int roi_z0, int roi_z1,
		  float vdr);
  void load_spacecharge(TH3F *hist, float zoffset, float scalefactor);
  void setScaleFactorB(float x){Bscale=x;return;};
  void setScaleFactorE(float x){Escale=x;return;};
  void setFlatFields(float B, float E);

  TVector3 calc_unit_field(TVector3 at, TVector3 from);
  TVector3 interpolatedFieldIntegral(float zdest,TVector3 start);
  int FilterPhiIndex(int phi,int range); //defaults to using nphi for range.

  TVector3 GetCellCenter(int r, int phi, int z);
  TVector3 GetGroupCellCenter(int r0, int r1, int phi0, int phi1, int z0, int z1);
  TVector3 GetWeightedCellCenter(int r, int phi, int z);
  TVector3 fieldIntegral(float zdest,TVector3 start);
  void populate_fieldmap();
  void  populate_lookup();
  TVector3 sum_field_at(int r,int phi, int z);
  TVector3 swimToInSteps(float zdest,TVector3 start, int steps, bool interpolate);
  TVector3 swimTo(float zdest,TVector3 start, bool interpolate);
 
 private:
  BoundsCase GetRindexAndCheckBounds(float pos, int *r);
  BoundsCase GetPhiIndexAndCheckBounds(float pos, int *phi);
  BoundsCase GetZindexAndCheckBounds(float pos, int *z);
  int FilterPhiIndex(int phi);

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
