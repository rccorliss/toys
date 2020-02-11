#!/usr/bin/env wolframscript
(* ::Package:: *)

Limu[mu_,x_]:=1/2*(BesselI[-I*mu,x]+BesselI[I*mu,x]) (**eqn 5.45**)
(*Rnk[munk_,n_,a_,b_]:=(Limu[munk,n*\[Pi]/L*a]*BesselK[I*munk,n*\[Pi]/L*b]
                      -BesselK[I*munk,n*\[Pi]/L*a]*Limu[munk,n*\[Pi]/L*b])
Rnk[munk,n,rmin,rmax]\[Equal]0 *)(**eqn 5.46**)
(* if we define aL=a/L and bL=b/L, this becomes independent of scale: *)
Rnk[munk_,n_,aL_,bL_]:=(Limu[munk,n*\[Pi]*aL]*BesselK[I*munk,n*\[Pi]*bL]
                      -BesselK[I*munk,n*\[Pi]*aL]*Limu[munk,n*\[Pi]*bL])
(*Rnk[munk,n,aL,bL]\[Equal]0*)
(*So we need the 'x' argument to go from 0 to, NumberOfOrders (15) times the aspect ratio, in steps of epsilon.
  and we need the mu argument to be...?  *)
  (*higher mu orders fluctuate more rapidly, so I will get zeros faster.  all of these have zeroes within 0<x<1.5*)


(* run the previous block, then this block to export limu and kimu tables with the specified range and number of entries *)
With[{mumin=0.0001,mumax=80,muN=400,
      xmin=0.2,xmax=80.0,xN=100},
     {Export[StringJoin[NotebookDirectory[],"limu_table.csv"],
             StringJoin[StringJoin["mu,",ToString[mumin],",",ToString[mumax],",",ToString[muN],"\n"],
                       StringJoin["x,",ToString[xmin],",",ToString[xmax],",",ToString[xN],"\n"],
                       ExportString[Flatten[Table[{mu,x,Re[Limu[mu,x]]},
                                                   {x,xmin,xmax,(xmax-xmin)/(xN-1)},
                                                   {mu,mumin,mumax,(mumax-mumin)/(muN-1)}],
                                            1](*close Flatten arguments*),
                                    "CSV"](*close ExportString arguments*)
                      ](*close StringJoin arguments*)
             ,"CSV","TextDelimiters"->""],(*close Export arguments -- limu stable now written*)
      Export[StringJoin[NotebookDirectory[],"kimu_table.csv"],
            StringJoin[StringJoin["mu,",ToString[mumin],",",ToString[mumax],",",ToString[muN],"\n"],
                       StringJoin["x,",ToString[xmin],",",ToString[xmax],",",ToString[xN],"\n"],
                       ExportString[Flatten[Table[{mu,x,Re[BesselK[I*mu,x]]},
                                                   {x,xmin,xmax,(xmax-xmin)/(xN-1)},
                                                   {mu,mumin,mumax,(mumax-mumin)/(muN-1)}],
                                            1](*close Flatten arguments*),
                                    "CSV"](*close ExportString arguments*)
                      ](*close StringJoin arguments*)
             ,"CSV","TextDelimiters"->""](*close Export arguments -- kimu stable now written*)
      }
    ]


With[{mumin=0.0001,mumax=80,muN=400,
      xmin=0.2,xmax=80.0,xN=100},
     {limutable=OpenWrite[StringJoin[NotebookDirectory[],"limu_table.csv"]];
      Write[limutable,StringJoin["mu,",ToString[mumin],",",ToString[mumax],",",ToString[muN]]];
      Write[limutable,StringJoin["x,",ToString[xmin],",",ToString[xmax],",",ToString[xN]];
                       ExportString[Flatten[Table[{mu,x,Re[Limu[mu,x]]},
                                                   {x,xmin,xmax,(xmax-xmin)/(xN-1)},
                                                   {mu,mumin,mumax,(mumax-mumin)/(muN-1)}],
                                            1](*close Flatten arguments*),
                                    "CSV"](*close ExportString arguments*)
                      ](*close StringJoin arguments*)
             ,"CSV","TextDelimiters"->""],(*close Export arguments -- limu stable now written*)
      Export[StringJoin[NotebookDirectory[],"kimu_table.csv"],
            StringJoin[StringJoin["mu,",ToString[mumin],",",ToString[mumax],",",ToString[muN],"\n"],
                       StringJoin["x,",ToString[xmin],",",ToString[xmax],",",ToString[xN],"\n"],
                       ExportString[Flatten[Table[{mu,x,Re[BesselK[I*mu,x]]},
                                                   {x,xmin,xmax,(xmax-xmin)/(xN-1)},
                                                   {mu,mumin,mumax,(mumax-mumin)/(muN-1)}],
                                            1](*close Flatten arguments*),
                                    "CSV"](*close ExportString arguments*)
                      ](*close StringJoin arguments*)
             ,"CSV","TextDelimiters"->""](*close Export arguments -- kimu stable now written*)
      }
    ]


(*remaining items are for debugging.*)



Plot[Evaluate[
Table[Rnk[mu,n,1.2,2.1],{n,5,10,1}]],{mu,30,100}]


Export["limu_table.csv",ExportString[Table[{mu,x,Re[Limu[mu,x]]},{x,0.1,1.5,0.1},{mu,1,10,1}],"CSV"]]


mumin=1
mumax=2
muN=3
muStep=(mumax-mumin)/muN
StringJoin["mu,",ToString[mumin],",",ToString[mumax],",",ToString[muN]]
ExportString[Flatten[Table[{"mu",mumin,mumax,muN},{muN,1,1,1}],1],"CSV"]
StringJoin[StringJoin["mu,",ToString[mumin],",",ToString[mumax],",",ToString[muN],"\n"],ExportString[Flatten[Table[{mu,x,Re[Limu[mu,x]]},{x,0.5,1.5,0.5},{mu,2,10,4}],1],"CSV"]]


Print[NotebookDirectory[],"hello"]


Export[StringJoin[NotebookDirectory[],"limu_table.csv"],
       ExportString[Flatten[Table[{mu,x,Re[Limu[mu,x]]},{mu,0.5,20,0.5},{x,0.05,3.0,0.05}],1],"CSV"],
       "CSV"]
Export[StringJoin[NotebookDirectory[],"kimu_table.csv"],
       ExportString[Flatten[Table[{mu,x,Re[BesselK[I*mu,x]]},{mu,0.5,20,0.5},{x,0.05,3.0,0.05}],1],"CSV"],
       "CSV"]


With[{mumin=0.5,mumax=4,muN=5,
      xmin=0.5,xmax=3.0,xN=6},
      StringJoin[StringJoin["mu,",ToString[mumin],",",ToString[mumax],",",ToString[muN],"\n"],
                 StringJoin["x,",ToString[xmin],",",ToString[xmax],",",ToString[xN],"\n"],
           ExportString[Flatten[Table[{mu,x,Re[Limu[mu,x]]},{x,xmin,xmax,(xmax-xmin)/(xN-1)},{mu,mumin,mumax,(mumax-mumin)/(muN-1)}],1],"CSV"]]]
