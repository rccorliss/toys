#include "FieldSim.h"

FieldSim::FieldSim(int nx, int ny, int nz, float omtau){
  //create a grid with the specified dimensions
  Efield=new MultiArray<TVector3>(nx,ny,nz);
  for (int i=0;i<Efield->GetLength();i++)
    Efield->GetFlat(i)->SetXYZ(0,0,0);

  //and a grid to compute the impact form each cell on each cell.
  Epartial=new MultiArray<TVector3>(nx,ny,nz,nx,ny,nz);
  for (int i=0;i<Epartial->GetLength();i++)
    Epartial->GetFlat(i)->SetXYZ(0,0,0);

  //and to hold the external fieldmap
  Eexternal=new MultiArray<TVector3>(nx,ny,nz);
  for (int i=0;i<Eexternal->GetLength();i++)
    Eexternal->GetFlat(i)->SetXYZ(0,0,0);

  //ditto the external magnetic fieldmap
  Bfield=new MultiArray<TVector3>(nx,ny,nz);
  for (int i=0;i<Bfield->GetLength();i++)
    Bfield->GetFlat(i)->SetXYZ(0,0,0);

  //and a grid of charges in those volumes.
  q=new MultiArray<float>(nx,ny,nz);
  for (int i=0;i<q->GetLength();i++)
    *(q->GetFlat(i))=0;
}


