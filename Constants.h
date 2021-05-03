#ifndef __CONSTANTS_H__
#define __CONSTANTS_H__
#include <string>

namespace Constants {
  //TPC dimensions
  //const float kInnerRadius = 83.5; //cm
  //const float kOutterRadius = 254.51; //cm
  //const float kHalfLength = 250.0; //cm
  const float kInnerRadius = 20.0; //cm
  const float kOutterRadius = 78.0; //cm
  const float kHalfLength = 105.5; //cm

  //number of steps for volume discretization
  const int kNRadialSteps = 159;
  const int kNAzimuthalSteps = 360;
  const int kNLongitudinalSteps = 62;

  //Running conditions
  const float kE0 = 400; //V/cm
  const float kB0 = 14; //kGauss

  //Langevin gas parameters
  // P9 gas
  //const float kT1 = 1.34; //NIM Phys Res A235,296 (1985)
  //const float kT2 = 1.11; //NIM Phys Res A235,296 (1985)
  // P10 gas
  //const float kT1 = 1.36; //measured by STAR
  //const float kT2 = 1.11; //measured by STAR
  // 85.7%Ne / 9.5%CO2 / 4.5%N2
  const float kT1 = 1.0; // measured by ALICE
  const float kT2 = 1.0; // measured by ALICE
  const float kV0 = 4.0; // [cm/us] @ E=400V/cm (need fine tunning)

  //info exchanger root file
  const std::string kFileNameRoot = "RccDistortions";
}

#endif /* __CONSTANTS_H__ */
