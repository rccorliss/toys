#!/usr/bin/env wolframscript
(* ::Package:: *)

(*"mirror 1 position and home orientation"*)
homeN1={1,0,0};
(*true value:  M1={2.53,1.03,4.19};*)
(*for testing: *)M1={5.03,0,3.97};
(*"mirror 2 home position and orientation"*)
homeN2={0,1,0};
homeM2={5.03,7.97,3.97};
(*"light pipe position"*)
(*true: P={4.87,8.13,4.19};*)
(*for testing:*) P=homeM2+{0,0,2.5};
(*"incoming laser direction:"*)
L0={-1,0,0};


mirrorNorm[in_,out_]:=Normalize[-in+out]
reflectBeam[in_,norm_]:=in-2(in.norm)norm
rotatedNorm[\[Alpha]_,\[Beta]_,\[Gamma]_,norm0_]:=RotationMatrix[\[Beta],{0,0,1}].RotationTransform[\[Alpha],{1,0,0}].RotationTransform[\[Gamma]{0,1,0}].norm0
positionAtZ[vector_,start_,z_]:=start+(z-start[[3]])/vector[[3]]*vector
\[Theta]ofBeam[L_]:=ArcCos[L[[3]]]
\[Phi]ofBeam[L_]:=ArcTan[L[[1]],L[[2]]](*note that the args are denom, numerator, not the other way around!*)
safe\[Phi]ofBeam[L_]:=If[L[[1]]==0,If[L[[2]]==0,0,\[Pi]/2],\[Phi]ofBeam[L]]
pitchOfMirror[N_]:=ArcSin[N[[3]]]
\[Beta]ofMirror[N_,Nhome_]:=ArcCos[Normalize[N-{0,0,N[[3]]}].Nhome]


(*"forward derivation:"*)
forwardN1[\[Gamma]_,\[Beta]_]:=rotatedNorm[0,\[Gamma],\[Beta],homeN1]
forwardN2[\[Alpha]_,\[Beta]_]:=rotatedNorm[\[Alpha],\[Gamma],\[Beta],homeN2]
forwardL1[N1_]:=reflectBeam[L0,N1]
forwardL1c[\[Gamma]_,\[Beta]_]:=forwardL1[forwardN1[\[Gamma],\[Beta]]]
forwardM2[L1_]:=positionAtZ[L1,M1,homeM2[[3]]]
forwardM2c[\[Gamma]_,\[Beta]_]:=forwardM2[forwardL1c[\[Gamma],\[Beta]]]
forwardL2[L1_,N2_]:=reflectBeam[L1,N2]
forwardL2c[\[Gamma]1_,\[Beta]1_,\[Alpha]2_,\[Beta]2_]:=forwardL2[forwardL1c[\[Gamma]1,\[Beta]1],forwardN2[\[Alpha]2,\[Beta]2]]

forward\[Theta][\[Gamma]1_,\[Beta]1_,\[Alpha]2_,\[Beta]2_]:=\[Theta]ofBeam[forwardL2c[\[Gamma]1,\[Beta]1,\[Alpha]2,\[Beta]2]]
forward\[Phi][\[Gamma]1_,\[Beta]1_,\[Alpha]2_,\[Beta]2_]:=\[Phi]ofBeam[forwardL2c[\[Gamma]1,\[Beta]1,\[Alpha]2,\[Beta]2]]



(*"backward derivation:"*)
backwardL2[\[Theta]_,\[Phi]_]:={Sin[\[Theta]]Cos[\[Phi]],Sin[\[Theta]]Sin[\[Phi]],Cos[\[Theta]]}
backwardM2[L2_]:=positionAtZ[L2,P,homeM2[[3]]]
backwardM2c[\[Theta]_,\[Phi]_]:=backwardM2[backwardL2[\[Theta],\[Phi]]]
backwardL1[M2_]:=Normalize[M2-M1]
backwardL1c[\[Theta]_,\[Phi]_]:=backwardL1[backwardM2c[\[Theta],\[Phi]]]
backwardN2[L1in_,L2out_]:=mirrorNorm[L1in,L2out]
backwardN2c[\[Theta]_,\[Phi]_]:=backwardN2[backwardL1c[\[Theta],\[Phi]],backwardL2[\[Theta],\[Phi]]]
backwardN1[L0in_,L1out_]:=mirrorNorm[L0in,L1out]
backwardN1c[\[Theta]_,\[Phi]_]:=backwardN1[L0,backwardL1c[\[Theta],\[Phi]]]

backward\[Beta]1[\[Theta]_,\[Phi]_]:=\[Beta]ofMirror[backwardN1c[\[Theta],\[Phi]],homeN1]
backward\[Gamma]1[\[Theta]_,\[Phi]_]:=pitchOfMirror[backwardN1c[\[Theta],\[Phi]]]
backward\[Beta]2[\[Theta]_,\[Phi]_]:=\[Beta]ofMirror[backwardN1c[\[Theta],\[Phi]],homeN2]
backward\[Alpha]2[\[Theta]_,\[Phi]_]:=pitchOfMirror[backwardN2c[\[Theta],\[Phi]]]


"calculations for home:"
L2ymax=Normalize[P-(homeM2)];
positionAtZ[L2ymax,(homeM2),P[[3]]]-P;
positionAtZ[L2ymax,P,homeM2[[3]]]-homeM2;
\[Theta]ymax=\[Theta]ofBeam[L2ymax];
\[Theta]ymax*180/\[Pi]"\[Degree]=\[Theta]"
\[Phi]ymax=safe\[Phi]ofBeam[L2ymax];
\[Phi]ymax*180/\[Pi]"\[Degree]=\[Phi]"
"M2-homeM2="
backwardM2c[\[Theta]ymax,\[Phi]ymax]-homeM2
"\[Degree]=\[Beta]1"backward\[Beta]1[\[Theta]ymax,\[Phi]ymax]*180/\[Pi]
"\[Degree]=\[Gamma]1"backward\[Gamma]1[\[Theta]ymax,\[Phi]ymax]*180/\[Pi]
"\[Degree]=\[Beta]2"backward\[Beta]2[\[Theta]ymax,\[Phi]ymax]*180/\[Pi]
"\[Degree]=\[Alpha]2"backward\[Alpha]2[\[Theta]ymax,\[Phi]ymax]*180/\[Pi]


"calculations for YMax:"
L2ymax=Normalize[P-(homeM2+{0,+2.50,0})];
positionAtZ[L2ymax,(homeM2+{0,+2.52,0}),P[[3]]]-P;
positionAtZ[L2ymax,P,homeM2[[3]]]-homeM2;
\[Theta]ymax=\[Theta]ofBeam[L2ymax];
\[Theta]ymax*180/\[Pi]"\[Degree]=\[Theta]"
\[Phi]ymax=\[Phi]ofBeam[L2ymax];
\[Phi]ymax*180/\[Pi]"\[Degree]=\[Phi]"
"M2-homeM2="
backwardM2c[\[Theta]ymax,\[Phi]ymax]-homeM2
"\[Degree]=\[Beta]1"backward\[Beta]1[\[Theta]ymax,\[Phi]ymax]*180/\[Pi]
"\[Degree]=\[Gamma]1"backward\[Gamma]1[\[Theta]ymax,\[Phi]ymax]*180/\[Pi]
"\[Degree]=\[Beta]2"backward\[Beta]2[\[Theta]ymax,\[Phi]ymax]*180/\[Pi]
"\[Degree]=\[Alpha]2"backward\[Alpha]2[\[Theta]ymax,\[Phi]ymax]*180/\[Pi]


"calculations for XYMax:"
L2xymax=Normalize[P-(homeM2+{-2.50,+2.50,0})];
\[Theta]xymax=\[Theta]ofBeam[L2xymax];
\[Theta]xymax*180/\[Pi]"\[Degree]=\[Theta]"
\[Phi]xymax=\[Phi]ofBeam[L2xymax];
\[Phi]xymax*180/\[Pi]"\[Degree]=\[Phi]"
"M2-homeM2="
backwardM2c[\[Theta]xymax,\[Phi]xymax]-homeM2
"\[Degree]=\[Beta]1"backward\[Beta]1[\[Theta]xymax,\[Phi]xymax]*180/\[Pi]
"\[Degree]=\[Gamma]1"backward\[Gamma]1[\[Theta]xymax,\[Phi]xymax]*180/\[Pi]
"\[Degree]=\[Beta]2"backward\[Beta]2[\[Theta]xymax,\[Phi]xymax]*180/\[Pi]
"\[Degree]=\[Alpha]2"backward\[Alpha]2[\[Theta]xymax,\[Phi]xymax]*180/\[Pi]


"calculations for XMax:"
L2xmax=Normalize[P-(homeM2+{-2.50,0,0})];
\[Theta]xmax=\[Theta]ofBeam[L2xmax];
\[Theta]xmax*180/\[Pi]"\[Degree]=\[Theta]"
\[Phi]xmax=\[Phi]ofBeam[L2xmax];
\[Phi]xmax*180/\[Pi]"\[Degree]=\[Phi]"
"M2-homeM2="
backwardM2c[\[Theta]ymax,\[Phi]xmax]-homeM2
"\[Degree]=\[Beta]1"backward\[Beta]1[\[Theta]xmax,\[Phi]xmax]*180/\[Pi]
"\[Degree]=\[Gamma]1"backward\[Gamma]1[\[Theta]xmax,\[Phi]xmax]*180/\[Pi]
"\[Degree]=\[Beta]2"backward\[Beta]2[\[Theta]xmax,\[Phi]xmax]*180/\[Pi]
"\[Degree]=\[Alpha]2"backward\[Alpha]2[\[Theta]xmax,\[Phi]xmax]*180/\[Pi]



