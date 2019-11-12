#include "assert.h"
#include "TVector3.h"

template <class T> class MultiArray;
class FieldSim{
 public: //bad form to leave this all exposed, I know.
  int nx,ny,nz; //dimensions of internal grid
  TVector3 step; //step size in each direction.
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


 public:
  FieldSim(float dx,float dy, float dz,int x,int y, int z, float vdr);
  void setScaleFactorB(float x){Bscale=x;return;};
  void setScaleFactorE(float x){Escale=x;return;};
  void setFlatFields(float B, float E);

  TVector3 calc_unit_field(TVector3 at, TVector3 from);
  TVector3 fieldIntegral(float zdest,TVector3 start);
  void populate_fieldmap();
  void  populate_lookup();
  TVector3 sum_field_at(int x,int y, int z);
  TVector3 swimTo(float zdest,TVector3 start);
 
 private:
  
};

template <class T>
class MultiArray{
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
