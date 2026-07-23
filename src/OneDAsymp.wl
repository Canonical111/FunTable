(* ::Package:: *)

If[ValueQ[$DistributedContexts],AppendTo[$DistributedContexts,"OneDAsymp`"],$DistributedContexts={"Global`","OneDAsymp`"}];


BeginPackage["OneDAsymp`"]

bosonicminusAsymp::usage = "bosonicminusAsymp[\!\(\*SubscriptBox[\(\[CapitalDelta]\), \(\[Phi]\)]\), m, \[CapitalDelta]] returns {betafunctionallist[\!\(\*SubscriptBox[\(\[CapitalDelta]\), \(\[Phi]\)]\), m, \[CapitalDelta]],  alphafunctionallist[\!\(\*SubscriptBox[\(\[CapitalDelta]\), \(\[Phi]\)]\), m, \[CapitalDelta]]}. for m from 0 to order";

fermionicplusAsymp::usage = "fermionicplusAsymp[\!\(\*SubscriptBox[\(\[CapitalDelta]\), \(\[Phi]\)]\), m, \[CapitalDelta]] returns {betafunctionallist[\!\(\*SubscriptBox[\(\[CapitalDelta]\), \(\[Phi]\)]\), m, \[CapitalDelta]],  alphafunctionallist[\!\(\*SubscriptBox[\(\[CapitalDelta]\), \(\[Phi]\)]\), m, \[CapitalDelta]]}. for m from 0 to order";

functionalpairAsymp::usage = "functionalpairAsymp[\[CapitalDelta]\[Phi],order,\[CapitalDelta]] returns the pair \[VeryThinSpace]{bosonicminusAsymp, fermionicplusAsymp} which is a 2*2*(order+1) dimensional list";
asymp1dfast;
asymp2dfast;
bosonicminussereva;
fermionicplussereva;
bosonicminusserevalist;
fermionicplusserevalist;
asymptablelist;


Begin["`Private`"]


\[Mu][n_]:=-n^2;
\[Nu][\[CapitalDelta]\[Phi]_,\[CapitalDelta]_][n_]:=(-1+\[CapitalDelta]) \[CapitalDelta]+\[CapitalDelta]\[Phi]+1/2 n (-1+n+4 \[CapitalDelta]\[Phi]);
\[Rho][\[CapitalDelta]\[Phi]_][n_]:=-(((n+2 \[CapitalDelta]\[Phi])^2 (-1+n+4 \[CapitalDelta]\[Phi])^2)/(4 (-1+2 n+4 \[CapitalDelta]\[Phi]) (1+2 n+4 \[CapitalDelta]\[Phi])));
\[Mu]p[n_]:=-2 n;
\[Nu]p[\[CapitalDelta]\[Phi]_][n_]:=1/2 (-1+n+2 \[CapitalDelta]\[Phi])+1/2 (n+2 \[CapitalDelta]\[Phi]);
\[Rho]p[\[CapitalDelta]\[Phi]_][n_]:=-(((n+2 \[CapitalDelta]\[Phi]) (-1+n+4 \[CapitalDelta]\[Phi]) (1+4 n^3-6 \[CapitalDelta]\[Phi]+24 n^2 \[CapitalDelta]\[Phi]+32 \[CapitalDelta]\[Phi]^3+n (-2+48 \[CapitalDelta]\[Phi]^2)))/(2 (-1+2 n+4 \[CapitalDelta]\[Phi])^2 (1+2 n+4 \[CapitalDelta]\[Phi])^2));
(* Same fix as in OneDfunctional.wl: the n=-1 coefficients always multiply b[-1]=0 / a[-1]=0,
   but their denominators contain (-3+4\[CapitalDelta]\[Phi]), which vanishes at \[CapitalDelta]\[Phi]=3/4 and turns
   0*ComplexInfinity into Indeterminate. The poisoned series then fall through the
   series2value[x_,...]:=0 fallback as SILENT ZEROS in the asymp region. *)
\[Rho][\[CapitalDelta]\[Phi]_][-1]=0;
\[Rho]p[\[CapitalDelta]\[Phi]_][-1]=0;

bs1[\[CapitalDelta]\[Phi]_,\[CapitalDelta]_][n_]=(2^(-2 (n+\[CapitalDelta]\[Phi])) Sqrt[\[Pi]] (n+\[CapitalDelta]\[Phi]) Gamma[n+\[CapitalDelta]\[Phi]]^2 Gamma[1+n+\[CapitalDelta]\[Phi]] Gamma[1/2+n+2 \[CapitalDelta]\[Phi]]^2)/((1+2 n-\[CapitalDelta]+2 \[CapitalDelta]\[Phi]) (2 n+\[CapitalDelta]+2 \[CapitalDelta]\[Phi]) Gamma[1+n]^2 Gamma[3/2+n+\[CapitalDelta]\[Phi]] Gamma[1/2+2 n+2 \[CapitalDelta]\[Phi]]);
as1[\[CapitalDelta]\[Phi]_,\[CapitalDelta]_][n_]=(16^(n+\[CapitalDelta]\[Phi]) Gamma[n+\[CapitalDelta]\[Phi]]^2 Gamma[1+n+\[CapitalDelta]\[Phi]]^2 Gamma[1/2+n+2 \[CapitalDelta]\[Phi]]^2 (-((-1+\[CapitalDelta]) \[CapitalDelta])-8 (n+\[CapitalDelta]\[Phi])^2 (1+2 n+2 \[CapitalDelta]\[Phi])+2 (-1-2 n+\[CapitalDelta]-2 \[CapitalDelta]\[Phi]) (n+\[CapitalDelta]\[Phi]) (1+2 n+2 \[CapitalDelta]\[Phi]) (2 n+\[CapitalDelta]+2 \[CapitalDelta]\[Phi]) (-(1/(n+\[CapitalDelta]\[Phi]))+HarmonicNumber[n]-2 HarmonicNumber[-1+n+\[CapitalDelta]\[Phi]]+2 HarmonicNumber[4 (n+\[CapitalDelta]\[Phi])]-HarmonicNumber[-(1/2)+n+2 \[CapitalDelta]\[Phi]]-Log[4])))/(Sqrt[\[Pi]] (1+2 n+2 \[CapitalDelta]\[Phi])^2 (1+2 n-\[CapitalDelta]+2 \[CapitalDelta]\[Phi])^2 (2 n+\[CapitalDelta]+2 \[CapitalDelta]\[Phi])^2 Gamma[1+n]^2 Gamma[1+4 n+4 \[CapitalDelta]\[Phi]]);
bs0[\[CapitalDelta]\[Phi]_,\[CapitalDelta]_][n_]=((-1)^(2 n) Sqrt[\[Pi]] Gamma[n+\[CapitalDelta]\[Phi]]^4 Gamma[-(1/2)+n+2 \[CapitalDelta]\[Phi]]^2)/((2 n-\[CapitalDelta]+2 \[CapitalDelta]\[Phi]) (-1+2 n+\[CapitalDelta]+2 \[CapitalDelta]\[Phi]) (n!)^2 Gamma[\[CapitalDelta]\[Phi]]^4 Gamma[2 (n+\[CapitalDelta]\[Phi])] Gamma[-(1/2)+2 n+2 \[CapitalDelta]\[Phi]]);
as0[\[CapitalDelta]\[Phi]_,\[CapitalDelta]_][n_]=((-1)^(1+2 n) Sqrt[\[Pi]] Gamma[n+\[CapitalDelta]\[Phi]]^4 Gamma[-(1/2)+n+2 \[CapitalDelta]\[Phi]]^2 (-1+4 n+4 \[CapitalDelta]\[Phi]+(2 n-\[CapitalDelta]+2 \[CapitalDelta]\[Phi]) (-1+2 n+\[CapitalDelta]+2 \[CapitalDelta]\[Phi]) (HarmonicNumber[n]-2 HarmonicNumber[-1+n+\[CapitalDelta]\[Phi]]-HarmonicNumber[-(3/2)+n+2 \[CapitalDelta]\[Phi]]+2 HarmonicNumber[-2+4 n+4 \[CapitalDelta]\[Phi]]-Log[4])))/((-2 n+\[CapitalDelta]-2 \[CapitalDelta]\[Phi])^2 (-1+2 n+\[CapitalDelta]+2 \[CapitalDelta]\[Phi])^2 (n!)^2 Gamma[\[CapitalDelta]\[Phi]]^4 Gamma[2 (n+\[CapitalDelta]\[Phi])] Gamma[-(1/2)+2 n+2 \[CapitalDelta]\[Phi]]);


R[\[CapitalDelta]\[Phi]_][n_]:=If[Mod[2 n,2]==0,(\[Pi]^(1/2) (-1)^(1-2n) Gamma[n+\[CapitalDelta]\[Phi]]^4 Gamma[-1/2+n+2\[CapitalDelta]\[Phi]]^2)/((n!)^2 Gamma[\[CapitalDelta]\[Phi]]^4 Gamma[2(n+\[CapitalDelta]\[Phi])]Gamma[2(n+\[CapitalDelta]\[Phi])-1/2]),0]
S[\[CapitalDelta]\[Phi]_][n_]:=If[Mod[2 n,2]==0,(\[Pi]^(1/2) (-1)^(-2n) Gamma[n+\[CapitalDelta]\[Phi]]^4 Gamma[-1/2+n+2\[CapitalDelta]\[Phi]]^2)/((n!)^2 Gamma[\[CapitalDelta]\[Phi]]^4 Gamma[2(n+\[CapitalDelta]\[Phi])]Gamma[2(n+\[CapitalDelta]\[Phi])-1/2]) (-2HarmonicNumber[n+\[CapitalDelta]\[Phi]-1]-HarmonicNumber[n+2\[CapitalDelta]\[Phi]-3/2]+2HarmonicNumber[4n+4\[CapitalDelta]\[Phi]-2]+HarmonicNumber[n]-Log[4]),0]
Rev[\[CapitalDelta]\[Phi]_][n_]:=-(((-1)^(2 n) 2^(-1-2 n-2 \[CapitalDelta]\[Phi]) Sqrt[\[Pi]] (6 n^2+3 n (-1+4 \[CapitalDelta]\[Phi])+\[CapitalDelta]\[Phi] (-1+4 \[CapitalDelta]\[Phi])) Gamma[n+\[CapitalDelta]\[Phi]]^3 Gamma[-(1/2)+n+2 \[CapitalDelta]\[Phi]]^2)/(Gamma[1+n]^2 Gamma[1/2+n+\[CapitalDelta]\[Phi]] Gamma[-(1/2)+2 n+2 \[CapitalDelta]\[Phi]]))
Rod[\[CapitalDelta]\[Phi]_][n_]:=-(((-1)^(2 n) 2^(-2+4 n+4 \[CapitalDelta]\[Phi]) Gamma[n+\[CapitalDelta]\[Phi]]^2 Gamma[1+n+\[CapitalDelta]\[Phi]]^2 Gamma[1/2+n+2 \[CapitalDelta]\[Phi]]^2)/(Sqrt[\[Pi]] (1+2 n+2 \[CapitalDelta]\[Phi]) Gamma[1+n]^2 Gamma[4 (n+\[CapitalDelta]\[Phi])]))
Sev[\[CapitalDelta]\[Phi]_][n_]:=If[n==0,-(1/2) \[CapitalDelta]\[Phi] Gamma[\[CapitalDelta]\[Phi]]^4 ((2 \[CapitalDelta]\[Phi] Gamma[\[CapitalDelta]\[Phi]]^2 Gamma[1/2+2 \[CapitalDelta]\[Phi]])/(Gamma[2 \[CapitalDelta]\[Phi]] Gamma[1+\[CapitalDelta]\[Phi]]^2)+((-1+4 \[CapitalDelta]\[Phi]) Gamma[-(1/2)+2 \[CapitalDelta]\[Phi]] (EulerGamma+2 PolyGamma[0,\[CapitalDelta]\[Phi]]-PolyGamma[0,2 \[CapitalDelta]\[Phi]]))/Gamma[2 \[CapitalDelta]\[Phi]]-(2 \[CapitalDelta]\[Phi] Gamma[1/2+2 \[CapitalDelta]\[Phi]] (EulerGamma+PolyGamma[0,\[CapitalDelta]\[Phi]]+PolyGamma[0,1+\[CapitalDelta]\[Phi]]-PolyGamma[0,1+2 \[CapitalDelta]\[Phi]]))/Gamma[1+2 \[CapitalDelta]\[Phi]]),-(1/2) (-1)^(2 n) \[CapitalDelta]\[Phi] Gamma[n+\[CapitalDelta]\[Phi]]^2 (-((2 Gamma[\[CapitalDelta]\[Phi]]^2 Gamma[n+\[CapitalDelta]\[Phi]]^2 Gamma[-(1/2)+n+2 \[CapitalDelta]\[Phi]] Gamma[1/2+n+2 \[CapitalDelta]\[Phi]] (2+\[CapitalDelta]\[Phi] PolyGamma[0,n]-2 \[CapitalDelta]\[Phi] PolyGamma[0,1+\[CapitalDelta]\[Phi]]+\[CapitalDelta]\[Phi] PolyGamma[0,1/2+n+2 \[CapitalDelta]\[Phi]]))/((-1+n)! n! Gamma[1+\[CapitalDelta]\[Phi]]^2 Gamma[2 (n+\[CapitalDelta]\[Phi])] Gamma[-(1/2)+2 (n+\[CapitalDelta]\[Phi])]))-(\[CapitalDelta]\[Phi] (-1+n+\[CapitalDelta]\[Phi]) Gamma[\[CapitalDelta]\[Phi]]^2 Gamma[-1+n+\[CapitalDelta]\[Phi]]^2 Gamma[-(1/2)+n+2 \[CapitalDelta]\[Phi]]^2 (1/(1-6 (n+\[CapitalDelta]\[Phi])+8 (n+\[CapitalDelta]\[Phi])^2) 2 (-(1/2)+n+\[CapitalDelta]\[Phi]) (-1+2 n+2 \[CapitalDelta]\[Phi]) (-3+4 n+4 \[CapitalDelta]\[Phi]) (3-4 n-4 \[CapitalDelta]\[Phi]-2 (-1+n+\[CapitalDelta]\[Phi]) (-(1/2)+n+\[CapitalDelta]\[Phi]) PolyGamma[0,n]+4 (-1+n+\[CapitalDelta]\[Phi]) (-(1/2)+n+\[CapitalDelta]\[Phi]) PolyGamma[0,\[CapitalDelta]\[Phi]]-2 (-1+n+\[CapitalDelta]\[Phi]) (-(1/2)+n+\[CapitalDelta]\[Phi]) PolyGamma[0,-1+n+\[CapitalDelta]\[Phi]]+2 (-1+n+\[CapitalDelta]\[Phi]) (-(1/2)+n+\[CapitalDelta]\[Phi]) PolyGamma[0,n+\[CapitalDelta]\[Phi]]-2 (-1+n+\[CapitalDelta]\[Phi]) (-(1/2)+n+\[CapitalDelta]\[Phi]) PolyGamma[0,-(1/2)+n+2 \[CapitalDelta]\[Phi]])-1/(1-4 (n+\[CapitalDelta]\[Phi]))^2 (4 (-1+n+\[CapitalDelta]\[Phi]) (-1+2 n+2 \[CapitalDelta]\[Phi])^2 (-3+4 n+4 \[CapitalDelta]\[Phi])-(-1+4 n+4 \[CapitalDelta]\[Phi]) (-11+48 n^3+54 \[CapitalDelta]\[Phi]-88 \[CapitalDelta]\[Phi]^2+48 \[CapitalDelta]\[Phi]^3+8 n^2 (-11+18 \[CapitalDelta]\[Phi])+2 n (27-88 \[CapitalDelta]\[Phi]+72 \[CapitalDelta]\[Phi]^2))+(-1+n+\[CapitalDelta]\[Phi]) (-1+2 n+2 \[CapitalDelta]\[Phi])^2 (-3+4 n+4 \[CapitalDelta]\[Phi]) (-1+4 n+4 \[CapitalDelta]\[Phi]) PolyGamma[0,n]+2 (-1+n+\[CapitalDelta]\[Phi]) (-1+2 n+2 \[CapitalDelta]\[Phi])^2 (-3+4 n+4 \[CapitalDelta]\[Phi]) (-1+4 n+4 \[CapitalDelta]\[Phi]) PolyGamma[0,\[CapitalDelta]\[Phi]]-3 (-1+n+\[CapitalDelta]\[Phi]) (-1+2 n+2 \[CapitalDelta]\[Phi])^2 (-3+4 n+4 \[CapitalDelta]\[Phi]) (-1+4 n+4 \[CapitalDelta]\[Phi]) PolyGamma[0,-1+n+\[CapitalDelta]\[Phi]]-(-1+n+\[CapitalDelta]\[Phi]) (-1+2 n+2 \[CapitalDelta]\[Phi])^2 (-3+4 n+4 \[CapitalDelta]\[Phi]) (-1+4 n+4 \[CapitalDelta]\[Phi]) PolyGamma[0,n+\[CapitalDelta]\[Phi]]-3 (-1+n+\[CapitalDelta]\[Phi]) (-1+2 n+2 \[CapitalDelta]\[Phi])^2 (-3+4 n+4 \[CapitalDelta]\[Phi]) (-1+4 n+4 \[CapitalDelta]\[Phi]) PolyGamma[0,-(1/2)+n+2 \[CapitalDelta]\[Phi]]+2 (-1+n+\[CapitalDelta]\[Phi]) (-1+2 n+2 \[CapitalDelta]\[Phi])^2 (-3+4 n+4 \[CapitalDelta]\[Phi]) (-1+4 n+4 \[CapitalDelta]\[Phi]) PolyGamma[0,-(3/2)+2 n+2 \[CapitalDelta]\[Phi]]+2 (-1+n+\[CapitalDelta]\[Phi]) (-1+2 n+2 \[CapitalDelta]\[Phi])^2 (-3+4 n+4 \[CapitalDelta]\[Phi]) (-1+4 n+4 \[CapitalDelta]\[Phi]) PolyGamma[0,-1+2 n+2 \[CapitalDelta]\[Phi]])))/((-1+2 n+2 \[CapitalDelta]\[Phi])^2 (-3+4 n+4 \[CapitalDelta]\[Phi])^2 ((-1+n)!)^2 Gamma[1+\[CapitalDelta]\[Phi]]^2 Gamma[-(3/2)+2 n+2 \[CapitalDelta]\[Phi]] Gamma[-1+2 n+2 \[CapitalDelta]\[Phi]])-(2 Gamma[\[CapitalDelta]\[Phi]]^2 Gamma[n+\[CapitalDelta]\[Phi]]^2 Gamma[-(1/2)+n+2 \[CapitalDelta]\[Phi]] Gamma[1/2+n+2 \[CapitalDelta]\[Phi]] (-2+\[CapitalDelta]\[Phi] PolyGamma[0,1+n]+2 \[CapitalDelta]\[Phi] PolyGamma[0,1+\[CapitalDelta]\[Phi]]-4 \[CapitalDelta]\[Phi] PolyGamma[0,n+\[CapitalDelta]\[Phi]]+2 \[CapitalDelta]\[Phi] PolyGamma[0,2 (n+\[CapitalDelta]\[Phi])]-\[CapitalDelta]\[Phi] PolyGamma[0,-(1/2)+n+2 \[CapitalDelta]\[Phi]]-2 \[CapitalDelta]\[Phi] PolyGamma[0,1/2+n+2 \[CapitalDelta]\[Phi]]+2 \[CapitalDelta]\[Phi] PolyGamma[0,-(1/2)+2 (n+\[CapitalDelta]\[Phi])]))/((-1+n)! n! Gamma[1+\[CapitalDelta]\[Phi]]^2 Gamma[2 (n+\[CapitalDelta]\[Phi])] Gamma[-(1/2)+2 (n+\[CapitalDelta]\[Phi])])-1/(n!)^2 ((2 \[CapitalDelta]\[Phi] Gamma[\[CapitalDelta]\[Phi]]^2 Gamma[1+n+\[CapitalDelta]\[Phi]]^2 Gamma[1/2+n+2 \[CapitalDelta]\[Phi]]^2 (-PolyGamma[0,1+n]+PolyGamma[0,n+\[CapitalDelta]\[Phi]]+PolyGamma[0,1+n+\[CapitalDelta]\[Phi]]+PolyGamma[0,1/2+n+2 \[CapitalDelta]\[Phi]]-PolyGamma[0,1/2+2 n+2 \[CapitalDelta]\[Phi]]-PolyGamma[0,1+2 n+2 \[CapitalDelta]\[Phi]]))/(Gamma[1+\[CapitalDelta]\[Phi]]^2 Gamma[1/2+2 n+2 \[CapitalDelta]\[Phi]] Gamma[1+2 n+2 \[CapitalDelta]\[Phi]])+((-1+4 \[CapitalDelta]\[Phi]) Gamma[n+\[CapitalDelta]\[Phi]]^2 Gamma[-(1/2)+n+2 \[CapitalDelta]\[Phi]]^2 (PolyGamma[0,1+n]-2 PolyGamma[0,n+\[CapitalDelta]\[Phi]]+PolyGamma[0,2 (n+\[CapitalDelta]\[Phi])]-PolyGamma[0,-(1/2)+n+2 \[CapitalDelta]\[Phi]]+PolyGamma[0,-(1/2)+2 (n+\[CapitalDelta]\[Phi])]))/(Gamma[2 (n+\[CapitalDelta]\[Phi])] Gamma[-(1/2)+2 (n+\[CapitalDelta]\[Phi])])))]

Sod[\[CapitalDelta]\[Phi]_][n_]:=((-1)^(2 n) 2^(-1+4 n+4 \[CapitalDelta]\[Phi]) Gamma[n+\[CapitalDelta]\[Phi]]^2 Gamma[1+n+\[CapitalDelta]\[Phi]]^2 Gamma[1/2+n+2 \[CapitalDelta]\[Phi]]^2 (-1+2 (n+2 n^2+\[CapitalDelta]\[Phi]+4 n \[CapitalDelta]\[Phi]+2 \[CapitalDelta]\[Phi]^2) PolyGamma[0,1+n]-2 (n+2 n^2+\[CapitalDelta]\[Phi]+4 n \[CapitalDelta]\[Phi]+2 \[CapitalDelta]\[Phi]^2) PolyGamma[0,n+\[CapitalDelta]\[Phi]]-2 n PolyGamma[0,1+n+\[CapitalDelta]\[Phi]]-4 n^2 PolyGamma[0,1+n+\[CapitalDelta]\[Phi]]-2 \[CapitalDelta]\[Phi] PolyGamma[0,1+n+\[CapitalDelta]\[Phi]]-8 n \[CapitalDelta]\[Phi] PolyGamma[0,1+n+\[CapitalDelta]\[Phi]]-4 \[CapitalDelta]\[Phi]^2 PolyGamma[0,1+n+\[CapitalDelta]\[Phi]]-2 n PolyGamma[0,1/2+n+2 \[CapitalDelta]\[Phi]]-4 n^2 PolyGamma[0,1/2+n+2 \[CapitalDelta]\[Phi]]-2 \[CapitalDelta]\[Phi] PolyGamma[0,1/2+n+2 \[CapitalDelta]\[Phi]]-8 n \[CapitalDelta]\[Phi] PolyGamma[0,1/2+n+2 \[CapitalDelta]\[Phi]]-4 \[CapitalDelta]\[Phi]^2 PolyGamma[0,1/2+n+2 \[CapitalDelta]\[Phi]]+2 n PolyGamma[0,1/2+2 n+2 \[CapitalDelta]\[Phi]]+4 n^2 PolyGamma[0,1/2+2 n+2 \[CapitalDelta]\[Phi]]+2 \[CapitalDelta]\[Phi] PolyGamma[0,1/2+2 n+2 \[CapitalDelta]\[Phi]]+8 n \[CapitalDelta]\[Phi] PolyGamma[0,1/2+2 n+2 \[CapitalDelta]\[Phi]]+4 \[CapitalDelta]\[Phi]^2 PolyGamma[0,1/2+2 n+2 \[CapitalDelta]\[Phi]]+2 n PolyGamma[0,1+2 n+2 \[CapitalDelta]\[Phi]]+4 n^2 PolyGamma[0,1+2 n+2 \[CapitalDelta]\[Phi]]+2 \[CapitalDelta]\[Phi] PolyGamma[0,1+2 n+2 \[CapitalDelta]\[Phi]]+8 n \[CapitalDelta]\[Phi] PolyGamma[0,1+2 n+2 \[CapitalDelta]\[Phi]]+4 \[CapitalDelta]\[Phi]^2 PolyGamma[0,1+2 n+2 \[CapitalDelta]\[Phi]]))/(Sqrt[\[Pi]] (1+2 n+2 \[CapitalDelta]\[Phi])^2 Gamma[1+n]^2 Gamma[1+4 n+4 \[CapitalDelta]\[Phi]]);





f1spin1sqexp[s_,\[CapitalDelta]_,\[CapitalDelta]\[Phi]_,qq_]:=(2^(1-2 s) Sqrt[\[Pi]] (-1+\[CapitalDelta]) (-1+4 s+\[CapitalDelta]-2 \[CapitalDelta]\[Phi]) Gamma[s]^3 Gamma[1/2+2 \[CapitalDelta]\[Phi]] Pochhammer[1,qq] Pochhammer[s,qq]^2 Pochhammer[-(1/2)+2 s-2 \[CapitalDelta]\[Phi],qq] Pochhammer[1/2+\[CapitalDelta]/2-\[CapitalDelta]\[Phi],qq] Pochhammer[3/4+s+\[CapitalDelta]/4-\[CapitalDelta]\[Phi]/2,qq])/((-1+2 s+\[CapitalDelta]-2 \[CapitalDelta]\[Phi])^2 (\[CapitalDelta]+2 \[CapitalDelta]\[Phi]) qq! Gamma[1/2+s] Pochhammer[2 s,qq] Pochhammer[1/2+s+\[CapitalDelta]/2-\[CapitalDelta]\[Phi],qq]^2 Pochhammer[-(1/4)+s+\[CapitalDelta]/4-\[CapitalDelta]\[Phi]/2,qq] Pochhammer[1+\[CapitalDelta]/2+\[CapitalDelta]\[Phi],qq])+(2^(3-2 s) Sqrt[\[Pi]] (-1+4 s+\[CapitalDelta]-2 \[CapitalDelta]\[Phi]) (s-\[CapitalDelta]\[Phi]) Gamma[s]^3 Gamma[1/2+2 \[CapitalDelta]\[Phi]] Pochhammer[1,qq] Pochhammer[s,qq]^2 Pochhammer[-(1/2)+2 s-2 \[CapitalDelta]\[Phi],qq] Pochhammer[1/2+\[CapitalDelta]/2-\[CapitalDelta]\[Phi],qq] Pochhammer[3/4+s+\[CapitalDelta]/4-\[CapitalDelta]\[Phi]/2,qq])/((-1+2 s+\[CapitalDelta]-2 \[CapitalDelta]\[Phi])^2 (\[CapitalDelta]+2 \[CapitalDelta]\[Phi]) qq! Gamma[1/2+s] Pochhammer[2 s,qq] Pochhammer[1/2+s+\[CapitalDelta]/2-\[CapitalDelta]\[Phi],qq]^2 Pochhammer[-(1/4)+s+\[CapitalDelta]/4-\[CapitalDelta]\[Phi]/2,qq] Pochhammer[1+\[CapitalDelta]/2+\[CapitalDelta]\[Phi],qq])+(2^(-2 s) Sqrt[\[Pi]] (1+\[CapitalDelta]-2 \[CapitalDelta]\[Phi])^2 (1+4 s+\[CapitalDelta]-2 \[CapitalDelta]\[Phi]) Gamma[s]^3 Gamma[-(1/2)+2 \[CapitalDelta]\[Phi]] Pochhammer[1,qq] Pochhammer[s,qq]^2 Pochhammer[1/2+2 s-2 \[CapitalDelta]\[Phi],qq] Pochhammer[3/2+\[CapitalDelta]/2-\[CapitalDelta]\[Phi],qq] Pochhammer[5/4+s+\[CapitalDelta]/4-\[CapitalDelta]\[Phi]/2,qq])/((1+2 s+\[CapitalDelta]-2 \[CapitalDelta]\[Phi])^2 (\[CapitalDelta]+2 \[CapitalDelta]\[Phi]) qq! Gamma[1/2+s] Pochhammer[2 s,qq] Pochhammer[3/2+s+\[CapitalDelta]/2-\[CapitalDelta]\[Phi],qq]^2 Pochhammer[1/4+s+\[CapitalDelta]/4-\[CapitalDelta]\[Phi]/2,qq] Pochhammer[1+\[CapitalDelta]/2+\[CapitalDelta]\[Phi],qq]);

f1spin1sqexpD[s_,\[CapitalDelta]_,\[CapitalDelta]\[Phi]_,qq_]=D[f1spin1sqexp[s,\[CapitalDelta],\[CapitalDelta]\[Phi],qq],s];

f1sqexp[s_,\[CapitalDelta]_,\[CapitalDelta]\[Phi]_,qq_]:=((2 Gamma[s]^4  Gamma[-(1/2)+2 \[CapitalDelta]\[Phi]] Pochhammer[1,qq] Pochhammer[s,qq]^2 Pochhammer[1/2+2 s-2 \[CapitalDelta]\[Phi],qq] Pochhammer[1+\[CapitalDelta]/2-\[CapitalDelta]\[Phi],qq] Pochhammer[1+s+\[CapitalDelta]/4-\[CapitalDelta]\[Phi]/2,qq])/(qq! Gamma[2 s]  Pochhammer[2 s,qq] Pochhammer[1+s+\[CapitalDelta]/2-\[CapitalDelta]\[Phi],qq]^2 Pochhammer[s+\[CapitalDelta]/4-\[CapitalDelta]\[Phi]/2,qq] Pochhammer[1/2+\[CapitalDelta]/2+\[CapitalDelta]\[Phi],qq])) (4 (4 s+\[CapitalDelta]-2 \[CapitalDelta]\[Phi]))/((2 s+\[CapitalDelta]-2 \[CapitalDelta]\[Phi])^2 (-1+\[CapitalDelta]+2 \[CapitalDelta]\[Phi]));

f1sqexpD[s_,\[CapitalDelta]_,\[CapitalDelta]\[Phi]_,qq_]=D[f1sqexp[s,\[CapitalDelta],\[CapitalDelta]\[Phi],qq],s];

\[Beta]0s1qexp[\[CapitalDelta]_,\[CapitalDelta]\[Phi]_,qqmax_]:=(2^(-2 \[CapitalDelta]\[Phi]-1) Sqrt[\[Pi]] Gamma[\[CapitalDelta]\[Phi]]^3 Gamma[-(1/2)+2 \[CapitalDelta]\[Phi]])/Gamma[1/2+\[CapitalDelta]\[Phi]]-1/2Sum[f1spin1sqexp[\[CapitalDelta]\[Phi],\[CapitalDelta],\[CapitalDelta]\[Phi],qq],{qq,0,qqmax}];
\[Alpha]0s1qexp[\[CapitalDelta]_,\[CapitalDelta]\[Phi]_,qqmax_]:=-((4^(-\[CapitalDelta]\[Phi]-1) Sqrt[\[Pi]] Gamma[\[CapitalDelta]\[Phi]]^3 Gamma[-(1/2)+2 \[CapitalDelta]\[Phi]] (Log[4]-3 PolyGamma[0,\[CapitalDelta]\[Phi]]+PolyGamma[0,1/2+\[CapitalDelta]\[Phi]]))/Gamma[1/2+\[CapitalDelta]\[Phi]])-1/4*Sum[f1spin1sqexpD[\[CapitalDelta]\[Phi],\[CapitalDelta],\[CapitalDelta]\[Phi],qq],{qq,0,qqmax}]+EulerGamma \[Beta]0s1qexp[\[CapitalDelta],\[CapitalDelta]\[Phi],qqmax];

\[Beta]0s0qexp[\[CapitalDelta]_,\[CapitalDelta]\[Phi]_,qqmax_]:=-(Sqrt[\[Pi]]/(8 Gamma[\[CapitalDelta]\[Phi]]^4))Sum[f1sqexp[\[CapitalDelta]\[Phi],\[CapitalDelta],\[CapitalDelta]\[Phi],qq],{qq,0,qqmax}]
\[Alpha]0s0qexp[\[CapitalDelta]_,\[CapitalDelta]\[Phi]_,qqmax_]:=-((1/2 Sqrt[\[Pi]]/(8 Gamma[\[CapitalDelta]\[Phi]]^4) Sum[f1sqexpD[\[CapitalDelta]\[Phi],\[CapitalDelta],\[CapitalDelta]\[Phi],qq],{qq,0,qqmax}])+EulerGamma Sqrt[\[Pi]]/(8 Gamma[\[CapitalDelta]\[Phi]]^4) Sum[f1sqexp[\[CapitalDelta]\[Phi],\[CapitalDelta],\[CapitalDelta]\[Phi],qq],{qq,0,qqmax}]);




\[Beta]mB[\[CapitalDelta]\[Phi]_,\[CapitalDelta]_][n_]:=bs0[\[CapitalDelta]\[Phi],\[CapitalDelta]][n]+2 b0[\[CapitalDelta]\[Phi],\[CapitalDelta]][2n]-(2^(1-2 n) Gamma[n+\[CapitalDelta]\[Phi]]^3 ((4^-\[CapitalDelta]\[Phi] \[Pi])/(\[CapitalDelta]-\[CapitalDelta]^2-2 \[CapitalDelta]\[Phi]+4 \[CapitalDelta]\[Phi]^2)+(b0[\[CapitalDelta]\[Phi],\[CapitalDelta]][0] Gamma[\[CapitalDelta]\[Phi]] Gamma[1/2+\[CapitalDelta]\[Phi]])/Gamma[-(1/2)+2 \[CapitalDelta]\[Phi]]) Gamma[-(1/2)+n+2 \[CapitalDelta]\[Phi]]^2)/(Gamma[1+n]^2 Gamma[\[CapitalDelta]\[Phi]]^4 Gamma[1/2+n+\[CapitalDelta]\[Phi]] Gamma[-(1/2)+2 n+2 \[CapitalDelta]\[Phi]]);
\[Alpha]mB[\[CapitalDelta]\[Phi]_,\[CapitalDelta]_][n_]:=(as0[\[CapitalDelta]\[Phi],\[CapitalDelta]][n]+2 a0[\[CapitalDelta]\[Phi],\[CapitalDelta]][2n]+(2^(-3-2 n+6 \[CapitalDelta]\[Phi]) Gamma[1/2+\[CapitalDelta]\[Phi]]^2 Gamma[n+\[CapitalDelta]\[Phi]]^3 Gamma[-(1/2)+n+2 \[CapitalDelta]\[Phi]]^2 (2 HarmonicNumber[n]-3 HarmonicNumber[-1+n+\[CapitalDelta]\[Phi]]+HarmonicNumber[-(1/2)+n+\[CapitalDelta]\[Phi]]-2 HarmonicNumber[-(3/2)+n+2 \[CapitalDelta]\[Phi]]+2 HarmonicNumber[-(3/2)+2 n+2 \[CapitalDelta]\[Phi]]+Log[4]))/(\[Pi] Gamma[1+n]^2 Gamma[\[CapitalDelta]\[Phi]]^2 Gamma[1/2+n+\[CapitalDelta]\[Phi]] Gamma[-(1/2)+2 n+2 \[CapitalDelta]\[Phi]] Gamma[-1+4 \[CapitalDelta]\[Phi]])b0[\[CapitalDelta]\[Phi],\[CapitalDelta]][0]+(2^(-2 (n+\[CapitalDelta]\[Phi])) \[Pi] Gamma[n+\[CapitalDelta]\[Phi]]^3 Gamma[-(1/2)+n+2 \[CapitalDelta]\[Phi]]^2 (2 HarmonicNumber[n]-3 HarmonicNumber[-1+n+\[CapitalDelta]\[Phi]]+HarmonicNumber[-(1/2)+n+\[CapitalDelta]\[Phi]]-2 HarmonicNumber[-(3/2)+n+2 \[CapitalDelta]\[Phi]]+2 HarmonicNumber[-(3/2)+2 n+2 \[CapitalDelta]\[Phi]]+Log[4]))/((\[CapitalDelta]-\[CapitalDelta]^2-2 \[CapitalDelta]\[Phi]+4 \[CapitalDelta]\[Phi]^2) Gamma[1+n]^2 Gamma[\[CapitalDelta]\[Phi]]^4 Gamma[1/2+n+\[CapitalDelta]\[Phi]] Gamma[-(1/2)+2 n+2 \[CapitalDelta]\[Phi]]));
\[Beta]pF[\[CapitalDelta]\[Phi]_,\[CapitalDelta]_][n_]:=bs1[\[CapitalDelta]\[Phi],\[CapitalDelta]][n]-2 b1[\[CapitalDelta]\[Phi],\[CapitalDelta]][2n+1];
\[Alpha]pF[\[CapitalDelta]\[Phi]_,\[CapitalDelta]_][n_]:=as1[\[CapitalDelta]\[Phi],\[CapitalDelta]][n]-2 a1[\[CapitalDelta]\[Phi],\[CapitalDelta]][2n+1];

defaultPrec=150;
Options[bosonicminusAsymp]=Options[fermionicplusAsymp]={WorkingPrecision->defaultPrec};

bosonicminusAsymp[ext_, order_,qmax_:46,  OptionsPattern[]]:=Block[{
\[CapitalDelta]\[Phi] = SetPrecision[ext,OptionValue[WorkingPrecision]], 
\[Beta]0s0Asymp,
\[Alpha]0s0Asymp,
b0,
a0
}, 
ClearAll[\[CapitalDelta]];
\[Beta]0s0Asymp=Series[\[Beta]0s0qexp[\[CapitalDelta],\[CapitalDelta]\[Phi],qmax],{\[CapitalDelta], \[Infinity], 2qmax+2}];\[Alpha]0s0Asymp=Series[\[Alpha]0s0qexp[\[CapitalDelta],\[CapitalDelta]\[Phi],qmax],{\[CapitalDelta], \[Infinity], 2qmax+2}];
a0[\[CapitalDelta]\[Phi]_,\[CapitalDelta]1_][-1]=0;
a0[\[CapitalDelta]\[Phi]_,\[CapitalDelta]1_][0]=\[Alpha]0s0Asymp;a0[\[CapitalDelta]\[Phi]_,\[CapitalDelta]1_][n_] := a0[\[CapitalDelta]\[Phi],\[CapitalDelta]1][n] = (S[\[CapitalDelta]\[Phi]][(n - 1)/2] - (\[Rho][\[CapitalDelta]\[Phi]][n-2] a0[\[CapitalDelta]\[Phi],\[CapitalDelta]1][n-2] +\[Nu][\[CapitalDelta]\[Phi],\[CapitalDelta]1][n-1] a0[\[CapitalDelta]\[Phi],\[CapitalDelta]1][n-1]+\[Rho]p[\[CapitalDelta]\[Phi]][n-2] b0[\[CapitalDelta]\[Phi],\[CapitalDelta]1][n-2]+\[Nu]p[\[CapitalDelta]\[Phi]][n-1] b0[\[CapitalDelta]\[Phi],\[CapitalDelta]1][n-1]+\[Mu]p[n] b0[\[CapitalDelta]\[Phi],\[CapitalDelta]1][n]))/\[Mu][n];b0[\[CapitalDelta]\[Phi]_,\[CapitalDelta]1_][-1]=0;
b0[\[CapitalDelta]\[Phi]_,\[CapitalDelta]1_][0]=\[Beta]0s0Asymp;b0[\[CapitalDelta]\[Phi]_,\[CapitalDelta]1_][n_]:=b0[\[CapitalDelta]\[Phi],\[CapitalDelta]1][n]=(R[\[CapitalDelta]\[Phi]][(n - 1)/2] -(\[Rho][\[CapitalDelta]\[Phi]][n - 2]*b0[\[CapitalDelta]\[Phi],\[CapitalDelta]1][n - 2] + \[Nu][\[CapitalDelta]\[Phi],\[CapitalDelta]1][n - 1]*b0[\[CapitalDelta]\[Phi],\[CapitalDelta]1][n - 1]))/\[Mu][n];{{0}~Join~( \[Beta]mB[\[CapitalDelta]\[Phi], \[CapitalDelta]]/@Range[order]),\[Alpha]mB[\[CapitalDelta]\[Phi], \[CapitalDelta]]/@Range[0,order]}
];

fermionicplusAsymp[ext_, order_,qmax_:46,  OptionsPattern[]]:=Block[{
\[CapitalDelta]\[Phi] = SetPrecision[ext,OptionValue[WorkingPrecision]], 
\[Beta]01Asymp,
\[Alpha]01Asymp,
b1,
a1
}, 
ClearAll[\[CapitalDelta]];
\[Beta]01Asymp=Series[\[Beta]0s1qexp[\[CapitalDelta],\[CapitalDelta]\[Phi],qmax],{\[CapitalDelta], \[Infinity], 2qmax}];\[Alpha]01Asymp=Series[\[Alpha]0s1qexp[\[CapitalDelta],\[CapitalDelta]\[Phi],qmax],{\[CapitalDelta], \[Infinity], 2qmax}];

a1[\[CapitalDelta]\[Phi]_,\[CapitalDelta]2_][-1] = 0; 
a1[\[CapitalDelta]\[Phi]_,\[CapitalDelta]2_][0]=\[Alpha]01Asymp; 
a1[\[CapitalDelta]\[Phi]_,\[CapitalDelta]2_][n_?OddQ]:= a1[\[CapitalDelta]\[Phi],\[CapitalDelta]2][n] = (Sev[\[CapitalDelta]\[Phi]][(n - 1)/2] -(\[Rho][\[CapitalDelta]\[Phi]][n - 2]*a1[\[CapitalDelta]\[Phi],\[CapitalDelta]2][n - 2] + \[Nu][\[CapitalDelta]\[Phi],\[CapitalDelta]2][n - 1]*a1[\[CapitalDelta]\[Phi],\[CapitalDelta]2][n - 1]+\[Rho]p[\[CapitalDelta]\[Phi]][n-2] b1[\[CapitalDelta]\[Phi],\[CapitalDelta]2][n-2]+\[Nu]p[\[CapitalDelta]\[Phi]][n-1] b1[\[CapitalDelta]\[Phi],\[CapitalDelta]2][n-1]+\[Mu]p[n] b1[\[CapitalDelta]\[Phi],\[CapitalDelta]2][n]))/\[Mu][n]; a1[\[CapitalDelta]\[Phi]_,\[CapitalDelta]2_][n_?EvenQ]:= a1[\[CapitalDelta]\[Phi],\[CapitalDelta]2][n] = (Sod[\[CapitalDelta]\[Phi]][(n - 2)/2] -(\[Rho][\[CapitalDelta]\[Phi]][n - 2]*a1[\[CapitalDelta]\[Phi],\[CapitalDelta]2][n - 2] + \[Nu][\[CapitalDelta]\[Phi],\[CapitalDelta]2][n - 1]*a1[\[CapitalDelta]\[Phi],\[CapitalDelta]2][n - 1]+\[Rho]p[\[CapitalDelta]\[Phi]][n-2] b1[\[CapitalDelta]\[Phi],\[CapitalDelta]2][n-2]+\[Nu]p[\[CapitalDelta]\[Phi]][n-1] b1[\[CapitalDelta]\[Phi],\[CapitalDelta]2][n-1]+\[Mu]p[n] b1[\[CapitalDelta]\[Phi],\[CapitalDelta]2][n]))/\[Mu][n]; 
b1[\[CapitalDelta]\[Phi]_,\[CapitalDelta]2_][-1]=0;
b1[\[CapitalDelta]\[Phi]_,\[CapitalDelta]2_][0]:=\[Beta]01Asymp;b1[\[CapitalDelta]\[Phi]_,\[CapitalDelta]2_][n_?OddQ]:= b1[\[CapitalDelta]\[Phi],\[CapitalDelta]2][n] = (Rev[\[CapitalDelta]\[Phi]][(n - 1)/2] -(\[Rho][\[CapitalDelta]\[Phi]][n - 2]*b1[\[CapitalDelta]\[Phi],\[CapitalDelta]2][n - 2] + \[Nu][\[CapitalDelta]\[Phi],\[CapitalDelta]2][n - 1]*b1[\[CapitalDelta]\[Phi],\[CapitalDelta]2][n - 1]))/\[Mu][n]; b1[\[CapitalDelta]\[Phi]_,\[CapitalDelta]2_][n_?EvenQ]:= b1[\[CapitalDelta]\[Phi],\[CapitalDelta]2][n] = (Rod[\[CapitalDelta]\[Phi]][(n - 2)/2] -(\[Rho][\[CapitalDelta]\[Phi]][n - 2]*b1[\[CapitalDelta]\[Phi],\[CapitalDelta]2][n - 2] + \[Nu][\[CapitalDelta]\[Phi],\[CapitalDelta]2][n - 1]*b1[\[CapitalDelta]\[Phi],\[CapitalDelta]2][n - 1]))/\[Mu][n];

{( \[Beta]pF[\[CapitalDelta]\[Phi], \[CapitalDelta]]/@Range[0,order]),\[Alpha]pF[\[CapitalDelta]\[Phi], \[CapitalDelta]]/@Range[0,order]}
];

bosonicminusAsymp[ext_, order_]:=Block[{
qmax=27+2*order,
prec
}, 
prec=2*qmax-10;
ClearAll[\[CapitalDelta]];
bosonicminusAsymp[ext,order,qmax,WorkingPrecision->prec]
];

fermionicplusAsymp[ext_, order_]:=Block[{
qmax=27+2*order,
prec
}, 
prec=2*qmax-10;
ClearAll[\[CapitalDelta]];
fermionicplusAsymp[ext,order,qmax,WorkingPrecision->prec]
];

\[Beta]\[Alpha]list={{1,1},{2,1},{1,2},{2,2}};
functionalindCut[\[CapitalLambda]_]:=Table[{nm,\[CapitalLambda]+2-nm},{nm,1,\[CapitalLambda]+1}];
functionalindCutOff[\[CapitalLambda]_]:=Flatten[#,1]&@Table[functionalindCut[i],{i,0,\[CapitalLambda]}];
functionallistCutOff[\[CapitalLambda]_]:=functionallistCutOff[\[CapitalLambda]]=Developer`ToPackedArray@DeleteCases[Flatten[Table[Flatten[Outer[{Flatten@{#1[[1]],#2[[1]]},Flatten@{#1[[2]],#2[[2]]}}&,\[Beta]\[Alpha]list,functionalindCut[i],1],1],{i,0,\[CapitalLambda]}],1],{{1,1},{___}},{1}];

asymp1dfast[\[CapitalDelta]\[Phi]N_,\[CapitalLambda]_,\[CapitalDelta]sermax_:28]:=Block[{bm=bosonicminusAsymp[\[CapitalDelta]\[Phi]N,\[CapitalLambda]],fp=fermionicplusAsymp[\[CapitalDelta]\[Phi]N,\[CapitalLambda]]},
Developer`ToPackedArray@N@{{PadLeft[bm[[1,2;;,3,;;\[CapitalDelta]sermax-5]],{\[CapitalLambda]+1,\[CapitalDelta]sermax-5}],bm[[2,All,3,;;\[CapitalDelta]sermax-5]]},fp[[All,All,3,;;\[CapitalDelta]sermax-5]]}
];
asymp2dfast[\[CapitalDelta]\[Phi]_,\[CapitalLambda]_,\[CapitalDelta]sermax_:28]:=Module[{tmp=asymp1dfast[\[CapitalDelta]\[Phi],\[CapitalLambda],\[CapitalDelta]sermax],minusasymp,plusasymp},
minusasymp=tmp[[1]];
plusasymp=tmp[[2]];
{Extract[minusasymp,functionallistCutOff[\[CapitalLambda]][[All,1]]],Extract[plusasymp,functionallistCutOff[\[CapitalLambda]][[All,2]]]}//Transpose];

bosonicminusfactor2prod[\[CapitalDelta]\[Phi]_,\[CapitalDelta]_]:=-((4^(2-\[CapitalDelta]) Gamma[2 \[CapitalDelta]] Gamma[\[CapitalDelta]\[Phi]]^4)/(Gamma[\[CapitalDelta]/2]^4 Gamma[1/2 (-1+\[CapitalDelta])+\[CapitalDelta]\[Phi]]^2 Gamma[-(\[CapitalDelta]/2)+\[CapitalDelta]\[Phi]]^2));
fermionicplusfactor2prod[\[CapitalDelta]\[Phi]_,\[CapitalDelta]_]:=-((2 \[CapitalDelta] Gamma[1/2+\[CapitalDelta]] Gamma[\[CapitalDelta]](-1+\[CapitalDelta]))/(Gamma[1/2 (1+\[CapitalDelta])]^4  Gamma[1/2 (1-\[CapitalDelta]+2 \[CapitalDelta]\[Phi])]^2 Gamma[1/2 (\[CapitalDelta]+2 \[CapitalDelta]\[Phi])]^2));
series2value[x_,\[CapitalDelta]_]:=0;
series2value[ser_SeriesData,\[CapitalDelta]_]:=# . Table[\[CapitalDelta]^(-5-i),{i,1,Length[#]}]&[ser[[3]]];
bosonicminussereva[\[CapitalDelta]\[Phi]_,order_,\[CapitalDelta]_]:=Map[series2value[#,\[CapitalDelta]]&,bosonicminusfactor2prod[\[CapitalDelta]\[Phi],\[CapitalDelta]]*bosonicminusAsymp[\[CapitalDelta]\[Phi],order],{2}];
fermionicplussereva[\[CapitalDelta]\[Phi]_,order_,\[CapitalDelta]_]:=Map[series2value[#,\[CapitalDelta]]&,fermionicplusfactor2prod[\[CapitalDelta]\[Phi],\[CapitalDelta]]*fermionicplusAsymp[\[CapitalDelta]\[Phi],order],{2}];
(*bosonicminusserevalist[\[CapitalDelta]\[Phi]_,order_,\[CapitalDelta]list_]:=Block[{bmlist=Map[series2value[#,\[CapitalDelta]]&,bosonicminusAsymp[\[CapitalDelta]\[Phi],order],{2}]},
ParallelTable[bmlist,{\[CapitalDelta],\[CapitalDelta]list}]Table[N[bosonicminusfactor2prod[\[CapitalDelta]\[Phi],\[CapitalDelta]],150],{\[CapitalDelta],\[CapitalDelta]list}]];
fermionicplusserevalist[\[CapitalDelta]\[Phi]_,order_,\[CapitalDelta]list_]:=Block[{fplist=Map[series2value[#,\[CapitalDelta]]&,fermionicplusAsymp[\[CapitalDelta]\[Phi],order],{2}]},
Table[fplist,{\[CapitalDelta],N@\[CapitalDelta]list}]Table[N[fermionicplusfactor2prod[\[CapitalDelta]\[Phi],\[CapitalDelta]],150],{\[CapitalDelta],\[CapitalDelta]list}]];*)

series2valuelist[x_,\[CapitalDelta]mat_]:=Table[0.,Last@Dimensions[\[CapitalDelta]mat]];
series2valuelist[ser_SeriesData,\[CapitalDelta]mat_]:=# . \[CapitalDelta]mat[[;;Length[#]]]&[ser[[3]]];

asymptablelist[\[CapitalDelta]\[Phi]_,order_,\[CapitalDelta]list_]:=Block[{
bm=bosonicminusAsymp[\[CapitalDelta]\[Phi],order],
fp=fermionicplusAsymp[\[CapitalDelta]\[Phi],order],
serlength,
dmatrix,
bmmat,
fpmat,
factormh=ParallelMap[bosonicminusfactor2prod[\[CapitalDelta]\[Phi],#]&,N[\[CapitalDelta]list,32]],
factorph=ParallelMap[fermionicplusfactor2prod[\[CapitalDelta]\[Phi],#]&,N[\[CapitalDelta]list,32]]
},
serlength=Length[bm[[2,1,3]]];
dmatrix=Table[Table[\[CapitalDelta]^(-5-i),{\[CapitalDelta],\[CapitalDelta]list}],{i,1,serlength}]; (* exact exponent keeps factormh in arb. prec., avoiding machine-precision overflow at large \[CapitalDelta] *)
bmmat=Transpose[Map[series2valuelist[#,dmatrix]&,bm,{2}],{2,3,1}]*factormh;
fpmat=Transpose[Map[series2valuelist[#,dmatrix]&,fp,{2}],{2,3,1}]*factorph;
{bmmat,fpmat}
];


End[]
EndPackage[]
