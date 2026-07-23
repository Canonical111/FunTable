(* ::Package:: *)

SetDirectory[DirectoryName[$InputFileName]];

Get[FileNameJoin[{"src", "OneDfunctional.wl"}]];

Get[FileNameJoin[{"src", "OneDAsymp.wl"}]];

ClearAll@generate1ptmx;
generate1ptmx[filename_String,dphi1_,orderMax1_,\[CapitalDelta]max_:250,normedFlag_:False]:=Block[{
dphi=dphi1,
orderMax=orderMax1,
table=Range[0,\[CapitalDelta]max,1/202],
\[CapitalDelta]table=Range[0,\[CapitalDelta]max,1/202],
tableLength,
lowseuil=2*dphi1+2*orderMax1+25+EulerGamma,
seuilind,
tablelow,
tablehigh,
prec=50,
precMax=0,
minushigh,
plushigh,
plus,
minus,
minuslow,
pluslow,
normed=normedFlag
},
Print["1d \[CapitalDelta]\[Phi]="<>ToString[N@dphi]];
seuilind=First@FirstPosition[table,x_/;x>lowseuil];
tablelow=table[[;;seuilind]];
tablehigh=table[[seuilind+1;;]];
tableLength=table//Length;
While[precMax<27,
prec=prec+10;
precMax=Precision@functionalpair[dphi,orderMax,tablelow[[-1]],PrecisionGoal->prec];
];
Print["Start"];
functionalbos[\[CapitalDelta]\[Phi]_,order_,\[CapitalDelta]_]:=(If[IntegerQ[\[CapitalDelta]/10],Print[\[CapitalDelta]]];bosonicminusfuncList[\[CapitalDelta]\[Phi],order,\[CapitalDelta],PrecisionGoal->prec]);     
minuslow=ParallelMap[functionalbos[dphi,orderMax,#]&,tablelow,Method->"FinestGrained"];
Print["Start Fer"];
functionalfer[\[CapitalDelta]\[Phi]_,order_,\[CapitalDelta]_]:=(If[IntegerQ[\[CapitalDelta]/10],Print[\[CapitalDelta]]];TimeConstrained[fermionicplusfuncList[\[CapitalDelta]\[Phi],order,\[CapitalDelta],PrecisionGoal->prec],2,Print[\[CapitalDelta]]]);pluslow=ParallelMap[functionalfer[dphi,orderMax,#]&,tablelow,Method->"FinestGrained"];
Print["Asymp"];
{minushigh,plushigh}=asymptablelist[dphi,orderMax,tablehigh];
plus=Join[pluslow,plushigh];
minus=Join[minuslow,minushigh];
dphi=N@dphi1;
If[normed,
(* exp-only rescaling in ARBITRARY PRECISION before converting to machine reals.
   This is essential: minus/plus at large \[CapitalDelta] grow as 4^\[CapitalDelta], which overflows doubles (~10^308)
   for \[CapitalDelta]>256. Rescaling first (while still in ~50-digit arb. prec.) gives
   2^(-2\[CapitalDelta]) \[Times] 4^\[CapitalDelta] \[Times] poly = poly, which is then safely converted by N[]. *)
minus=N[2^(-2 \[CapitalDelta]table)*minus];
plus=N[2^(-2 \[CapitalDelta]table)*plus];
,
minus=N@minus;
plus=N@plus;
];
DumpSave[filename,{dphi,precMax,orderMax,\[CapitalDelta]table,minus,plus,normed}];
];

Get[FileNameJoin[{"src", "Table.wl"}]];
