(* ::Package:: *)

If[ValueQ[$DistributedContexts],AppendTo[$DistributedContexts,"OneDfunctional`"],$DistributedContexts={"Global`","OneDfunctional`"}];

BeginPackage["OneDfunctional`"]

betaminus::usage = "betaminus[\!\(\*SubscriptBox[\(\[CapitalDelta]\), \(\[Phi]\)]\), m, \[CapitalDelta]] returns \!\(\*SubscriptBox[\(\[Beta]\), \(m\)]\)(\[CapitalDelta]) at external dimension \!\(\*SubscriptBox[\(\[CapitalDelta]\), \(\[Phi]\)]\).";
alphaminus::usage = "alphaminus[\!\(\*SubscriptBox[\(\[CapitalDelta]\), \(\[Phi]\)]\), m, \[CapitalDelta]] returns \!\(\*SubscriptBox[\(\[Alpha]\), \(m\)]\)(\[CapitalDelta]) at external dimension \!\(\*SubscriptBox[\(\[CapitalDelta]\), \(\[Phi]\)]\).";

betaplus::usage = "betaplus[\!\(\*SubscriptBox[\(\[CapitalDelta]\), \(\[Phi]\)]\), m, \[CapitalDelta]] returns \!\(\*SubscriptBox[\(\[Beta]\), \(m\)]\)(\[CapitalDelta]) at external dimension \!\(\*SubscriptBox[\(\[CapitalDelta]\), \(\[Phi]\)]\).";
alphaplus::usage = "alphaplus[\!\(\*SubscriptBox[\(\[CapitalDelta]\), \(\[Phi]\)]\), m, \[CapitalDelta]] returns \!\(\*SubscriptBox[\(\[Alpha]\), \(m\)]\)(\[CapitalDelta]) at external dimension \!\(\*SubscriptBox[\(\[CapitalDelta]\), \(\[Phi]\)]\).";

bosonicminusfuncList::usage = "bosonicminusfuncList[\!\(\*SubscriptBox[\(\[CapitalDelta]\), \(\[Phi]\)]\), m, \[CapitalDelta]] returns {betafunctionallist[\!\(\*SubscriptBox[\(\[CapitalDelta]\), \(\[Phi]\)]\), m, \[CapitalDelta]],  alphafunctionallist[\!\(\*SubscriptBox[\(\[CapitalDelta]\), \(\[Phi]\)]\), m, \[CapitalDelta]]}. for m from 0 to order";

fermionicplusfuncList::usage = "bosonicminusfuncList[\!\(\*SubscriptBox[\(\[CapitalDelta]\), \(\[Phi]\)]\), m, \[CapitalDelta]] returns {betafunctionallist[\!\(\*SubscriptBox[\(\[CapitalDelta]\), \(\[Phi]\)]\), m, \[CapitalDelta]],  alphafunctionallist[\!\(\*SubscriptBox[\(\[CapitalDelta]\), \(\[Phi]\)]\), m, \[CapitalDelta]]}. for m from 0 to order";

functionalpair::usage = "functionalpair[\[CapitalDelta]\[Phi],order,\[CapitalDelta]] returns the pair \[VeryThinSpace]{bosonicminusfuncList, fermionicplusfuncList} which is a 2*2*(order+1) dimensional list";

betatildezero::usage = "betatildezero[\!\(\*SubscriptBox[\(\[CapitalDelta]\), \(\[Phi]\)]\), \[CapitalDelta]] returns \!\(\*SubscriptBox[OverscriptBox[\(\[Beta]\), \(~\)], \(0\)]\)[\!\(\*SubscriptBox[\(\[CapitalDelta]\), \(\[Phi]\)]\), \[CapitalDelta]].";

agff::usage = "agff[\[CapitalDelta]\[Phi]_,\[CapitalDelta]_] gives the OPE^2 of \[CapitalDelta]\[Phi] at \[CapitalDelta]";

betaminusSym::usage="";
(*tilde0D2N::usage = "tilde0D2N[a,n] returns \!\(\*SuperscriptBox[SubscriptBox[OverscriptBox[\(\[Beta]\), \(~\)], \(0\)], \(\*\"\\\"\<\>\"\)]\) at extdim a and dimension 2a+2n for n strict positive integer";
betafunctionalderi2::usage="betafunctionalderi2[a,m,n] returns \!\(\*SuperscriptBox[SubscriptBox[\(\[Beta]\), \(m\)], \(\*\"\\\"\<\>\"\)]\) at extdim a and dimension 2a+2n for n strict positive integer";
alphatildefunctional::usage="alphatildefunctional[a,m,\[CapitalDelta]] returns \!\(\*SubscriptBox[OverscriptBox[\(\[Alpha]\), \(~\)], \(m\)]\)(\[CapitalDelta]).";
betatildefunctional::usage="betatildefunctional[a,m,\[CapitalDelta]] returns \!\(\*SubscriptBox[OverscriptBox[\(\[Beta]\), \(~\)], \(m\)]\)(\[CapitalDelta]).";
cn::usage="cn[\!\(\*SubscriptBox[\(\[CapitalDelta]\), \(\[Phi]\)]\)][n] return \!\(\*SubscriptBox[\(c\), \(n\)]\).";
dn::usage="dn[\!\(\*SubscriptBox[\(\[CapitalDelta]\), \(\[Phi]\)]\)][n] return \!\(\*SubscriptBox[\(d\), \(n\)]\).";
afermion::usage="afermion[\!\(\*SubscriptBox[\(\[CapitalDelta]\), \(\[Phi]\)]\)][m] return generalized free fermion OPE^2 \!\(\*SubscriptBox[\(a\), \(m\)]\) at 2\!\(\*SubscriptBox[\(\[CapitalDelta]\), \(\[Phi]\)]\)+2m+1.";
am::usage="am[\!\(\*SubscriptBox[\(\[CapitalDelta]\), \(\[Phi]\)]\)][m] return generalized free boson OPE^2 \!\(\*SubscriptBox[\(a\), \(m\)]\) at 2\!\(\*SubscriptBox[\(\[CapitalDelta]\), \(\[Phi]\)]\)+2m.";
cb1::usage="cb1[\!\(\*SubscriptBox[\(\[CapitalDelta]\), \(\[Phi]\)]\)][\[CapitalDelta],z] return 1d Conformal block.";
functionallistPara::usage="";
betaminusSym::usage="";*)
cn::usage="cn[\!\(\*SubscriptBox[\(\[CapitalDelta]\), \(\[Phi]\)]\)][n] return \!\(\*SubscriptBox[\(c\), \(n\)]\).";
dn::usage="dn[\!\(\*SubscriptBox[\(\[CapitalDelta]\), \(\[Phi]\)]\)][n] return \!\(\*SubscriptBox[\(d\), \(n\)]\).";


Begin["`Private`"]


\[Beta]0s0func[\[CapitalDelta]\[Phi]_,\[CapitalDelta]_]:=-(1/2) Gamma[\[CapitalDelta]/2]^2 Gamma[-(1/2)+\[CapitalDelta]/2+\[CapitalDelta]\[Phi]] Gamma[-(1/2)+2 \[CapitalDelta]\[Phi]] Gamma[-1+\[CapitalDelta]+2 \[CapitalDelta]\[Phi]] HypergeometricPFQRegularized[{-(1/2)+\[CapitalDelta]/2+\[CapitalDelta]\[Phi],1/2+\[CapitalDelta]/2,1/2+\[CapitalDelta]/2,-(1/2)+2 \[CapitalDelta]\[Phi],-1+\[CapitalDelta]+2 \[CapitalDelta]\[Phi]},{1/2+\[CapitalDelta],1/2+\[CapitalDelta]/2+\[CapitalDelta]\[Phi],-(1/2)+\[CapitalDelta]/2+2 \[CapitalDelta]\[Phi],-(1/2)+\[CapitalDelta]/2+2 \[CapitalDelta]\[Phi]},1];

\[Beta]0s0OverAfunc[\[CapitalDelta]\[Phi]_,\[CapitalDelta]_]:=-((2^(2-\[CapitalDelta]+2 \[CapitalDelta]\[Phi]) Gamma[2 \[CapitalDelta]] Gamma[\[CapitalDelta]\[Phi]]^4 Gamma[\[CapitalDelta]/2+\[CapitalDelta]\[Phi]] Gamma[-(1/2)+2 \[CapitalDelta]\[Phi]] HypergeometricPFQRegularized[{1/2 (-1+\[CapitalDelta]+2 \[CapitalDelta]\[Phi]),(1+\[CapitalDelta])/2,(1+\[CapitalDelta])/2,-(1/2)+2 \[CapitalDelta]\[Phi],-1+\[CapitalDelta]+2 \[CapitalDelta]\[Phi]},{1/2+\[CapitalDelta],1/2 (1+\[CapitalDelta]+2 \[CapitalDelta]\[Phi]),1/2 (-1+\[CapitalDelta]+4 \[CapitalDelta]\[Phi]),1/2 (-1+\[CapitalDelta]+4 \[CapitalDelta]\[Phi])},1])/(Sqrt[\[Pi]] Gamma[\[CapitalDelta]/2]^2 Gamma[-(\[CapitalDelta]/2)+\[CapitalDelta]\[Phi]]^2));

(*\[Beta]0s0OverAfunc[\[CapitalDelta]\[Phi]_,0]:=0;*)

bs0OverA[\[CapitalDelta]\[Phi]_,\[CapitalDelta]_][n_]:=(2 Gamma[\[CapitalDelta]] Gamma[1/2+\[CapitalDelta]] Gamma[n+\[CapitalDelta]\[Phi]]^4 Gamma[-(1/2)+n+2 \[CapitalDelta]\[Phi]]^2 Pochhammer[-(\[CapitalDelta]/2)+\[CapitalDelta]\[Phi],n] Pochhammer[-(1/2)+\[CapitalDelta]/2+\[CapitalDelta]\[Phi],n])/((n!)^2 Gamma[\[CapitalDelta]/2]^4 Gamma[-(\[CapitalDelta]/2)+\[CapitalDelta]\[Phi]] Gamma[1+n-\[CapitalDelta]/2+\[CapitalDelta]\[Phi]] Gamma[-(1/2)+\[CapitalDelta]/2+\[CapitalDelta]\[Phi]] Gamma[1/2+n+\[CapitalDelta]/2+\[CapitalDelta]\[Phi]] Gamma[-(1/2)+2 n+2 \[CapitalDelta]\[Phi]] Gamma[2 n+2 \[CapitalDelta]\[Phi]]);
as0OverA[\[CapitalDelta]\[Phi]_,\[CapitalDelta]_][n_]:=-((Gamma[\[CapitalDelta]] Gamma[1/2+\[CapitalDelta]] Gamma[n+\[CapitalDelta]\[Phi]]^4 Gamma[-(1/2)+n+2 \[CapitalDelta]\[Phi]]^2 Pochhammer[1/2 (-1+\[CapitalDelta])+\[CapitalDelta]\[Phi],n]^2 Pochhammer[-(\[CapitalDelta]/2)+\[CapitalDelta]\[Phi],n]^2 (-1+4 n+4 \[CapitalDelta]\[Phi]+(2 n-\[CapitalDelta]+2 \[CapitalDelta]\[Phi]) (-1+2 n+\[CapitalDelta]+2 \[CapitalDelta]\[Phi]) (HarmonicNumber[n]-2 HarmonicNumber[-1+n+\[CapitalDelta]\[Phi]]-HarmonicNumber[-(3/2)+n+2 \[CapitalDelta]\[Phi]]+2 HarmonicNumber[-2+4 n+4 \[CapitalDelta]\[Phi]]-Log[4])))/(2 (n!)^2 Gamma[\[CapitalDelta]/2]^4 Gamma[2 (n+\[CapitalDelta]\[Phi])] Gamma[1+n-\[CapitalDelta]/2+\[CapitalDelta]\[Phi]]^2 Gamma[1/2+n+\[CapitalDelta]/2+\[CapitalDelta]\[Phi]]^2 Gamma[-(1/2)+2 n+2 \[CapitalDelta]\[Phi]]));
A[\[CapitalDelta]\[Phi]_,\[CapitalDelta]_]:=(4^(-2+\[CapitalDelta])*Gamma[\[CapitalDelta]/2]^4*Gamma[(-1+\[CapitalDelta])/2+\[CapitalDelta]\[Phi]]^2*Gamma[-1/2*\[CapitalDelta]+\[CapitalDelta]\[Phi]]^2)/(Gamma[2*\[CapitalDelta]]*Gamma[\[CapitalDelta]\[Phi]]^4);

bCS[\[CapitalDelta]\[Phi]_,\[CapitalDelta]_,\[Beta]0s0OverA_][n_]:=(2^(-2 (n+\[CapitalDelta]\[Phi])) Sqrt[\[Pi]] Gamma[n+\[CapitalDelta]\[Phi]]^3 ((4^-\[CapitalDelta] Sqrt[\[Pi]] (\[CapitalDelta]-2 \[CapitalDelta]\[Phi]) (-1+\[CapitalDelta]+2 \[CapitalDelta]\[Phi]) Gamma[2 \[CapitalDelta]])/(Gamma[\[CapitalDelta]/2]^4 Gamma[1-\[CapitalDelta]/2+\[CapitalDelta]\[Phi]]^2 Gamma[(1+\[CapitalDelta])/2+\[CapitalDelta]\[Phi]]^2)-(\[Beta]0s0OverA Gamma[2 \[CapitalDelta]\[Phi]])/(Gamma[\[CapitalDelta]\[Phi]]^4 Gamma[-(1/2)+2 \[CapitalDelta]\[Phi]])) Gamma[-(1/2)+n+2 \[CapitalDelta]\[Phi]]^2)/(Gamma[1+n]^2 Gamma[1/2+n+\[CapitalDelta]\[Phi]] Gamma[-(1/2)+2 n+2 \[CapitalDelta]\[Phi]]);
aCS[\[CapitalDelta]\[Phi]_,\[CapitalDelta]_,\[Beta]0s0OverA_][n_]:=(2 HarmonicNumber[n]-3 HarmonicNumber[-1+n+\[CapitalDelta]\[Phi]]+HarmonicNumber[-(1/2)+n+\[CapitalDelta]\[Phi]]-2 HarmonicNumber[-(3/2)+n+2 \[CapitalDelta]\[Phi]]+2 HarmonicNumber[-(3/2)+2 n+2 \[CapitalDelta]\[Phi]]+Log[4])((4^(-1-n) \[Beta]0s0OverA (-1+4 \[CapitalDelta]\[Phi]) Gamma[1/2+\[CapitalDelta]\[Phi]] Gamma[n+\[CapitalDelta]\[Phi]]^3 Gamma[-(1/2)+n+2 \[CapitalDelta]\[Phi]]^2)/(Gamma[1+n]^2 Gamma[\[CapitalDelta]\[Phi]]^3 Gamma[1/2+n+\[CapitalDelta]\[Phi]] Gamma[1/2+2 \[CapitalDelta]\[Phi]] Gamma[-(1/2)+2 n+2 \[CapitalDelta]\[Phi]])-(4^(-1-n-\[CapitalDelta]) Sqrt[\[Pi]] (\[CapitalDelta]-2 \[CapitalDelta]\[Phi]) (-1+\[CapitalDelta]+2 \[CapitalDelta]\[Phi]) (-1+4 \[CapitalDelta]\[Phi]) Gamma[2 \[CapitalDelta]] Gamma[\[CapitalDelta]\[Phi]] Gamma[1/2+\[CapitalDelta]\[Phi]] Gamma[n+\[CapitalDelta]\[Phi]]^3 Gamma[-(1/2)+2 \[CapitalDelta]\[Phi]] Gamma[-(1/2)+n+2 \[CapitalDelta]\[Phi]]^2)/(Gamma[1+n]^2 Gamma[\[CapitalDelta]/2]^4 Gamma[2 \[CapitalDelta]\[Phi]] Gamma[1/2+n+\[CapitalDelta]\[Phi]] Gamma[1-\[CapitalDelta]/2+\[CapitalDelta]\[Phi]]^2 Gamma[1/2+\[CapitalDelta]/2+\[CapitalDelta]\[Phi]]^2 Gamma[1/2+2 \[CapitalDelta]\[Phi]] Gamma[-(1/2)+2 n+2 \[CapitalDelta]\[Phi]]));

bs1[\[CapitalDelta]\[Phi]_,\[CapitalDelta]_][n_]=-((2^(-2 n+\[CapitalDelta]-2 \[CapitalDelta]\[Phi]) (-1-2 n+\[CapitalDelta]-2 \[CapitalDelta]\[Phi]) (n+\[CapitalDelta]\[Phi]) Gamma[1+\[CapitalDelta]/2] Gamma[1/2+\[CapitalDelta]] Gamma[n+\[CapitalDelta]\[Phi]]^2 Gamma[1+n+\[CapitalDelta]\[Phi]] Gamma[1/2+n+2 \[CapitalDelta]\[Phi]]^2 Pochhammer[1/2-\[CapitalDelta]/2+\[CapitalDelta]\[Phi],n]^2)/((2 n+\[CapitalDelta]+2 \[CapitalDelta]\[Phi]) Gamma[1+n]^2 Gamma[1/2 (-1+\[CapitalDelta])] Gamma[(1+\[CapitalDelta])/2]^2 Gamma[3/2+n+\[CapitalDelta]\[Phi]] Gamma[3/2+n-\[CapitalDelta]/2+\[CapitalDelta]\[Phi]]^2 Gamma[\[CapitalDelta]/2+\[CapitalDelta]\[Phi]]^2 Gamma[1/2+2 n+2 \[CapitalDelta]\[Phi]]));as1[\[CapitalDelta]\[Phi]_,\[CapitalDelta]_][n_]=(2^(\[CapitalDelta]+4 (n+\[CapitalDelta]\[Phi])) Gamma[1+\[CapitalDelta]/2] Gamma[1/2+\[CapitalDelta]] Gamma[n+\[CapitalDelta]\[Phi]]^2 Gamma[1+n+\[CapitalDelta]\[Phi]]^2 Gamma[1/2+n+2 \[CapitalDelta]\[Phi]]^2 ((1-\[CapitalDelta]) \[CapitalDelta]-8 (n+\[CapitalDelta]\[Phi])^2 (1+2 n+2 \[CapitalDelta]\[Phi])+2 (-1-2 n+\[CapitalDelta]-2 \[CapitalDelta]\[Phi]) (n+\[CapitalDelta]\[Phi]) (1+2 n+2 \[CapitalDelta]\[Phi]) (2 n+\[CapitalDelta]+2 \[CapitalDelta]\[Phi]) (-(1/(n+\[CapitalDelta]\[Phi]))+HarmonicNumber[n]-2 HarmonicNumber[-1+n+\[CapitalDelta]\[Phi]]+2 HarmonicNumber[4 (n+\[CapitalDelta]\[Phi])]-HarmonicNumber[-(1/2)+n+2 \[CapitalDelta]\[Phi]]-Log[4])) Pochhammer[1/2-\[CapitalDelta]/2+\[CapitalDelta]\[Phi],n]^2)/(\[Pi] (1+2 n+2 \[CapitalDelta]\[Phi])^2 (2 n+\[CapitalDelta]+2 \[CapitalDelta]\[Phi])^2 Gamma[1+n]^2 Gamma[1/2 (-1+\[CapitalDelta])] Gamma[(1+\[CapitalDelta])/2]^2 Gamma[3/2+n-\[CapitalDelta]/2+\[CapitalDelta]\[Phi]]^2 Gamma[\[CapitalDelta]/2+\[CapitalDelta]\[Phi]]^2 Gamma[1+4 n+4 \[CapitalDelta]\[Phi]]);
\[Alpha]0s0prefunc[\[CapitalDelta]\[Phi]_,\[CapitalDelta]_][s_]:=(2^(-1+\[CapitalDelta]+2 \[CapitalDelta]\[Phi]) Gamma[s]^4 Gamma[1/2+\[CapitalDelta]] Gamma[s+\[CapitalDelta]/2-\[CapitalDelta]\[Phi]]^2 Gamma[\[CapitalDelta]/2+\[CapitalDelta]\[Phi]] HypergeometricPFQRegularized[{-(1/2)+\[CapitalDelta]/2+\[CapitalDelta]\[Phi],1/2-s+\[CapitalDelta]/2+\[CapitalDelta]\[Phi],1/2-s+\[CapitalDelta]/2+\[CapitalDelta]\[Phi],-(1/2)+2 \[CapitalDelta]\[Phi],-1+\[CapitalDelta]+2 \[CapitalDelta]\[Phi]},{1/2+\[CapitalDelta],1/2+\[CapitalDelta]/2+\[CapitalDelta]\[Phi],-(1/2)+s+\[CapitalDelta]/2+\[CapitalDelta]\[Phi],-(1/2)+s+\[CapitalDelta]/2+\[CapitalDelta]\[Phi]},1])/(Sqrt[\[Pi]] Gamma[1/2+2 s-2 \[CapitalDelta]\[Phi]]);

\[Beta]01func[\[CapitalDelta]\[Phi]_,0]:=
(4^(-\[CapitalDelta]-\[CapitalDelta]\[Phi]) \[Pi] (-1+\[CapitalDelta]) Gamma[1+2 \[CapitalDelta]] Gamma[\[CapitalDelta]\[Phi]]^3 Gamma[-(1/2)+2 \[CapitalDelta]\[Phi]])/(Gamma[1/2+\[CapitalDelta]/2]^4 Gamma[1/2+\[CapitalDelta]\[Phi]] Gamma[1/2-\[CapitalDelta]/2+\[CapitalDelta]\[Phi]]^2 Gamma[\[CapitalDelta]/2+\[CapitalDelta]\[Phi]]^2)-(2^(-1+\[CapitalDelta]+2 \[CapitalDelta]\[Phi]) (-1+\[CapitalDelta]) (-(1/2)-\[CapitalDelta]/2+\[CapitalDelta]\[Phi])^2 Gamma[1/2+\[CapitalDelta]] Gamma[1+\[CapitalDelta]] Gamma[\[CapitalDelta]\[Phi]]^4 Gamma[1/2+\[CapitalDelta]/2+\[CapitalDelta]\[Phi]] Gamma[-(1/2)+2 \[CapitalDelta]\[Phi]] HypergeometricPFQRegularized[{\[CapitalDelta]/2+\[CapitalDelta]\[Phi],1+\[CapitalDelta]/2,1+\[CapitalDelta]/2,-(1/2)+2 \[CapitalDelta]\[Phi],\[CapitalDelta]+2 \[CapitalDelta]\[Phi]},{3/2+\[CapitalDelta],1+\[CapitalDelta]/2+\[CapitalDelta]\[Phi],\[CapitalDelta]/2+2 \[CapitalDelta]\[Phi],\[CapitalDelta]/2+2 \[CapitalDelta]\[Phi]},1])/(\[Pi] Gamma[1/2+\[CapitalDelta]/2]^2 Gamma[1/2-\[CapitalDelta]/2+\[CapitalDelta]\[Phi]]^2)-(2 Gamma[1/2+\[CapitalDelta]] Gamma[1+\[CapitalDelta]] Gamma[3/4+\[CapitalDelta]/4+\[CapitalDelta]\[Phi]/2] Gamma[\[CapitalDelta]\[Phi]]^4 Gamma[1/2+2 \[CapitalDelta]\[Phi]] HypergeometricPFQRegularized[{1,\[CapitalDelta]\[Phi],\[CapitalDelta]\[Phi],-(1/2),1/2+\[CapitalDelta]/2-\[CapitalDelta]\[Phi],3/4+\[CapitalDelta]/4+\[CapitalDelta]\[Phi]/2},{2 \[CapitalDelta]\[Phi],1/2+\[CapitalDelta]/2,1/2+\[CapitalDelta]/2,-(1/4)+\[CapitalDelta]/4+\[CapitalDelta]\[Phi]/2,1+\[CapitalDelta]/2+\[CapitalDelta]\[Phi]},1])/(Gamma[1/2+\[CapitalDelta]/2]^2 Gamma[1/2-\[CapitalDelta]/2+\[CapitalDelta]\[Phi]]^2 Gamma[\[CapitalDelta]/2+\[CapitalDelta]\[Phi]])/.\[CapitalDelta]->0;

\[Beta]01func[\[CapitalDelta]\[Phi]_,\[CapitalDelta]_]:=(4^(-\[CapitalDelta]-\[CapitalDelta]\[Phi]) \[Pi] (-1+\[CapitalDelta]) Gamma[1+2 \[CapitalDelta]] Gamma[\[CapitalDelta]\[Phi]]^3 Gamma[-(1/2)+2 \[CapitalDelta]\[Phi]])/(Gamma[1/2+\[CapitalDelta]/2]^4 Gamma[1/2+\[CapitalDelta]\[Phi]] Gamma[1/2-\[CapitalDelta]/2+\[CapitalDelta]\[Phi]]^2 Gamma[\[CapitalDelta]/2+\[CapitalDelta]\[Phi]]^2)-(Gamma[\[CapitalDelta]] Gamma[1/2+\[CapitalDelta]] Gamma[1+\[CapitalDelta]] Gamma[\[CapitalDelta]\[Phi]]^4 HypergeometricPFQRegularized[{-(1/2),-1+\[CapitalDelta],1/2+\[CapitalDelta]/2-\[CapitalDelta]\[Phi],1/2+\[CapitalDelta]/2-\[CapitalDelta]\[Phi],-(1/2)+\[CapitalDelta]/2},{1/2+\[CapitalDelta],1/2+\[CapitalDelta]/2,-(1/2)+\[CapitalDelta]/2+\[CapitalDelta]\[Phi],-(1/2)+\[CapitalDelta]/2+\[CapitalDelta]\[Phi]},1])/(Gamma[1/2+\[CapitalDelta]/2]^3 Gamma[1/2-\[CapitalDelta]/2+\[CapitalDelta]\[Phi]]^2)-(2^(-1+\[CapitalDelta]+2 \[CapitalDelta]\[Phi]) (-1+\[CapitalDelta]) Gamma[1/2+\[CapitalDelta]] Gamma[1+\[CapitalDelta]] Gamma[\[CapitalDelta]\[Phi]]^4 Gamma[1/2+\[CapitalDelta]/2+\[CapitalDelta]\[Phi]] Gamma[-(1/2)+2 \[CapitalDelta]\[Phi]] HypergeometricPFQRegularized[{\[CapitalDelta]/2+\[CapitalDelta]\[Phi],1+\[CapitalDelta]/2,1+\[CapitalDelta]/2,-(1/2)+2 \[CapitalDelta]\[Phi],\[CapitalDelta]+2 \[CapitalDelta]\[Phi]},{3/2+\[CapitalDelta],1+\[CapitalDelta]/2+\[CapitalDelta]\[Phi],\[CapitalDelta]/2+2 \[CapitalDelta]\[Phi],\[CapitalDelta]/2+2 \[CapitalDelta]\[Phi]},1])/(\[Pi] Gamma[1/2+\[CapitalDelta]/2]^2 Gamma[-(1/2)-\[CapitalDelta]/2+\[CapitalDelta]\[Phi]]^2);
\[Alpha]01prefunc[\[CapitalDelta]\[Phi]_,\[CapitalDelta]_][s_]:=-(1/(2 (-1+2 s+\[CapitalDelta]-2 \[CapitalDelta]\[Phi])^2 Gamma[1/2+\[CapitalDelta]/2]^4 Gamma[1/2-\[CapitalDelta]/2+\[CapitalDelta]\[Phi]]^2))(-1+\[CapitalDelta]) (-1+4 s+\[CapitalDelta]-4 \[CapitalDelta]\[Phi]) Gamma[s]^4 Gamma[1/2+\[CapitalDelta]] Gamma[1+\[CapitalDelta]] Gamma[2 s+\[CapitalDelta]-2 \[CapitalDelta]\[Phi]] Gamma[1/2+s+\[CapitalDelta]/2-\[CapitalDelta]\[Phi]] HypergeometricPFQRegularized[{-(1/2)+2 s-2 \[CapitalDelta]\[Phi],-1+2 s+\[CapitalDelta]-2 \[CapitalDelta]\[Phi],1/2+\[CapitalDelta]/2-\[CapitalDelta]\[Phi],1/2+\[CapitalDelta]/2-\[CapitalDelta]\[Phi],-(1/2)+s+\[CapitalDelta]/2-\[CapitalDelta]\[Phi]},{1/2+\[CapitalDelta],1/2+s+\[CapitalDelta]/2-\[CapitalDelta]\[Phi],-(1/2)+2 s+\[CapitalDelta]/2-\[CapitalDelta]\[Phi],-(1/2)+2 s+\[CapitalDelta]/2-\[CapitalDelta]\[Phi]},1]-1/(Sqrt[\[Pi]] Gamma[1/2+\[CapitalDelta]/2]^4 Gamma[1/2+2 s-2 \[CapitalDelta]\[Phi]] Gamma[-(1/2)-\[CapitalDelta]/2+\[CapitalDelta]\[Phi]]^2) 2^(-2+\[CapitalDelta]+2 \[CapitalDelta]\[Phi]) (-1+\[CapitalDelta]) Gamma[s]^4 Gamma[1/2+\[CapitalDelta]] Gamma[1+\[CapitalDelta]] Gamma[1/2+s+\[CapitalDelta]/2-\[CapitalDelta]\[Phi]]^2 Gamma[1/2+\[CapitalDelta]/2+\[CapitalDelta]\[Phi]] Gamma[-(1/2)+2 \[CapitalDelta]\[Phi]] HypergeometricPFQRegularized[{\[CapitalDelta]/2+\[CapitalDelta]\[Phi],1-s+\[CapitalDelta]/2+\[CapitalDelta]\[Phi],1-s+\[CapitalDelta]/2+\[CapitalDelta]\[Phi],-(1/2)+2 \[CapitalDelta]\[Phi],\[CapitalDelta]+2 \[CapitalDelta]\[Phi]},{3/2+\[CapitalDelta],1+\[CapitalDelta]/2+\[CapitalDelta]\[Phi],s+\[CapitalDelta]/2+\[CapitalDelta]\[Phi],s+\[CapitalDelta]/2+\[CapitalDelta]\[Phi]},1];
\[Alpha]01funcAdd[\[CapitalDelta]\[Phi]_,\[CapitalDelta]_]:=-((2^(-1+\[CapitalDelta]-2 \[CapitalDelta]\[Phi]) (-1+\[CapitalDelta]) Gamma[1+\[CapitalDelta]/2] Gamma[1/2+\[CapitalDelta]] Gamma[\[CapitalDelta]\[Phi]]^3 Gamma[-(1/2)+2 \[CapitalDelta]\[Phi]] (Log[4]-3 PolyGamma[0,\[CapitalDelta]\[Phi]]+PolyGamma[0,1/2+\[CapitalDelta]\[Phi]]))/(Gamma[(1+\[CapitalDelta])/2]^3 Gamma[1/2+\[CapitalDelta]\[Phi]] Gamma[1/2-\[CapitalDelta]/2+\[CapitalDelta]\[Phi]]^2 Gamma[\[CapitalDelta]/2+\[CapitalDelta]\[Phi]]^2));


\[Alpha]01prefunc[\[CapitalDelta]\[Phi]_,0][s_]:=(2^(-3+2 \[CapitalDelta]\[Phi]) (-1+2 \[CapitalDelta]\[Phi]) Gamma[s]^4 Gamma[1/2+s-\[CapitalDelta]\[Phi]]^2 Gamma[-(1/2)+2 \[CapitalDelta]\[Phi]] HypergeometricPFQRegularized[{\[CapitalDelta]\[Phi],1-s+\[CapitalDelta]\[Phi],1-s+\[CapitalDelta]\[Phi],-(1/2)+2 \[CapitalDelta]\[Phi],2 \[CapitalDelta]\[Phi]},{3/2,1+\[CapitalDelta]\[Phi],s+\[CapitalDelta]\[Phi],s+\[CapitalDelta]\[Phi]},1])/(\[Pi]^2 Gamma[1/2+2 s-2 \[CapitalDelta]\[Phi]] Gamma[-(1/2)+\[CapitalDelta]\[Phi]])+(2^(-3+2 \[CapitalDelta]\[Phi]) (-1+4 s-4 \[CapitalDelta]\[Phi]) Gamma[s]^4 Gamma[-(1/2)+s-\[CapitalDelta]\[Phi]]^2 Gamma[3/4+s-\[CapitalDelta]\[Phi]/2] Gamma[1/2+2 \[CapitalDelta]\[Phi]] HypergeometricPFQRegularized[{1,s,s,-(1/2)+2 s-2 \[CapitalDelta]\[Phi],1/2-\[CapitalDelta]\[Phi],3/4+s-\[CapitalDelta]\[Phi]/2},{2 s,1/2+s-\[CapitalDelta]\[Phi],1/2+s-\[CapitalDelta]\[Phi],-(1/4)+s-\[CapitalDelta]\[Phi]/2,1+\[CapitalDelta]\[Phi]},1])/(\[Pi]^2 Gamma[2 \[CapitalDelta]\[Phi]] Gamma[1/2+\[CapitalDelta]\[Phi]]);
\[Mu][n_]:=-n^2;
\[Nu][\[CapitalDelta]\[Phi]_,\[CapitalDelta]_][n_]:=(-1+\[CapitalDelta]) \[CapitalDelta]+\[CapitalDelta]\[Phi]+1/2 n (-1+n+4 \[CapitalDelta]\[Phi]);
\[Rho][\[CapitalDelta]\[Phi]_,\[CapitalDelta]_][n_]:=-(((n+2 \[CapitalDelta]\[Phi])^2 (-1+n+4 \[CapitalDelta]\[Phi])^2)/(4 (-1+2 n+4 \[CapitalDelta]\[Phi]) (1+2 n+4 \[CapitalDelta]\[Phi])));
\[Mu]p[n_]:=-2 n;
\[Nu]p[\[CapitalDelta]\[Phi]_,\[CapitalDelta]_][n_]:=1/2 (-1+n+2 \[CapitalDelta]\[Phi])+1/2 (n+2 \[CapitalDelta]\[Phi]);
\[Rho]p[\[CapitalDelta]\[Phi]_,\[CapitalDelta]_][n_]:=-(((n+2 \[CapitalDelta]\[Phi]) (-1+n+4 \[CapitalDelta]\[Phi]) (1+4 n^3-6 \[CapitalDelta]\[Phi]+24 n^2 \[CapitalDelta]\[Phi]+32 \[CapitalDelta]\[Phi]^3+n (-2+48 \[CapitalDelta]\[Phi]^2)))/(2 (-1+2 n+4 \[CapitalDelta]\[Phi])^2 (1+2 n+4 \[CapitalDelta]\[Phi])^2));
(* The n=-1 coefficients always multiply b[-1]=0 / a[-1]=0 in the recursions, but their
   denominators contain (-3+4\[CapitalDelta]\[Phi]), which vanishes at \[CapitalDelta]\[Phi]=3/4 and turns 0*ComplexInfinity
   into Indeterminate. They must be exactly 0. *)
\[Rho][\[CapitalDelta]\[Phi]_,\[CapitalDelta]_][-1]=0;
\[Rho]p[\[CapitalDelta]\[Phi]_,\[CapitalDelta]_][-1]=0;

R[\[CapitalDelta]\[Phi]_,\[CapitalDelta]_][n_]:=If[Mod[2 n,2]==0,-(\[Pi]^(1/2)  Gamma[n+\[CapitalDelta]\[Phi]]^4 Gamma[-1/2+n+2\[CapitalDelta]\[Phi]]^2)/((n!)^2 Gamma[\[CapitalDelta]\[Phi]]^4 Gamma[2(n+\[CapitalDelta]\[Phi])]Gamma[2(n+\[CapitalDelta]\[Phi])-1/2]),0];
S[\[CapitalDelta]\[Phi]_,\[CapitalDelta]_][n_]:=If[Mod[2 n,2]==0,(\[Pi]^(1/2)    Gamma[n+\[CapitalDelta]\[Phi]]^4 Gamma[-1/2+n+2\[CapitalDelta]\[Phi]]^2)/((n!)^2 Gamma[\[CapitalDelta]\[Phi]]^4 Gamma[2(n+\[CapitalDelta]\[Phi])]Gamma[2(n+\[CapitalDelta]\[Phi])-1/2]) (-2HarmonicNumber[n+\[CapitalDelta]\[Phi]-1]-HarmonicNumber[n+2\[CapitalDelta]\[Phi]-3/2]+2HarmonicNumber[4n+4\[CapitalDelta]\[Phi]-2]+HarmonicNumber[n]-Log[4]),0];

Reven[\[CapitalDelta]\[Phi]_,\[CapitalDelta]_][n_]:=((-1)^(2 n) 2^(-2 (1+n+\[CapitalDelta]+\[CapitalDelta]\[Phi])) (-1+\[CapitalDelta]) Gamma[1+2 \[CapitalDelta]] Gamma[n+\[CapitalDelta]\[Phi]]^3 Gamma[-(1/2)+n+2 \[CapitalDelta]\[Phi]]^2 (-((64^(n+\[CapitalDelta]\[Phi]) n (14 n^2+2 (1-4 \[CapitalDelta]\[Phi])^2+n (-11+46 \[CapitalDelta]\[Phi])) Gamma[n+\[CapitalDelta]\[Phi]])/Gamma[4 (n+\[CapitalDelta]\[Phi])])+(2 \[Pi] (n+4 n^3-4 n \[CapitalDelta]\[Phi]-(1-4 \[CapitalDelta]\[Phi])^2 \[CapitalDelta]\[Phi]+4 n^2 (-1+5 \[CapitalDelta]\[Phi])))/(Gamma[1/2+n+\[CapitalDelta]\[Phi]] Gamma[1/2+2 n+2 \[CapitalDelta]\[Phi]])))/(Gamma[1+n]^2 Gamma[(1+\[CapitalDelta])/2]^4 Gamma[1/2-\[CapitalDelta]/2+\[CapitalDelta]\[Phi]]^2 Gamma[\[CapitalDelta]/2+\[CapitalDelta]\[Phi]]^2);
Rodd[\[CapitalDelta]\[Phi]_,\[CapitalDelta]_][n_]:=((-1)^(1+2 n) 2^(4 n+\[CapitalDelta]+4 \[CapitalDelta]\[Phi]) Gamma[1+\[CapitalDelta]/2] Gamma[1/2+\[CapitalDelta]] Gamma[n+\[CapitalDelta]\[Phi]]^2 Gamma[1+n+\[CapitalDelta]\[Phi]]^2 Gamma[1/2+n+2 \[CapitalDelta]\[Phi]]^2)/(\[Pi] (1+2 n+2 \[CapitalDelta]\[Phi]) Gamma[1+n]^2 Gamma[1/2 (-1+\[CapitalDelta])] Gamma[(1+\[CapitalDelta])/2]^2 Gamma[4 (n+\[CapitalDelta]\[Phi])] Gamma[1/2-\[CapitalDelta]/2+\[CapitalDelta]\[Phi]]^2 Gamma[\[CapitalDelta]/2+\[CapitalDelta]\[Phi]]^2);
Sodd[\[CapitalDelta]\[Phi]_,\[CapitalDelta]_][n_]:=(2^(1+4 n+\[CapitalDelta]+4 \[CapitalDelta]\[Phi]) Gamma[1+\[CapitalDelta]/2] Gamma[1/2+\[CapitalDelta]] Gamma[n+\[CapitalDelta]\[Phi]]^2 Gamma[1+n+\[CapitalDelta]\[Phi]]^2 Gamma[1/2+n+2 \[CapitalDelta]\[Phi]]^2 (-3-4 n-4 \[CapitalDelta]\[Phi]+2 (n+\[CapitalDelta]\[Phi]) (1+2 n+2 \[CapitalDelta]\[Phi]) (HarmonicNumber[n]-2 HarmonicNumber[-1+n+\[CapitalDelta]\[Phi]]+2 HarmonicNumber[4 (n+\[CapitalDelta]\[Phi])]-HarmonicNumber[-(1/2)+n+2 \[CapitalDelta]\[Phi]]-Log[4])))/(\[Pi] (1+2 n+2 \[CapitalDelta]\[Phi])^2 Gamma[1+n]^2 Gamma[1/2 (-1+\[CapitalDelta])] Gamma[(1+\[CapitalDelta])/2]^2 Gamma[1/2-\[CapitalDelta]/2+\[CapitalDelta]\[Phi]]^2 Gamma[\[CapitalDelta]/2+\[CapitalDelta]\[Phi]]^2 Gamma[1+4 n+4 \[CapitalDelta]\[Phi]]);
Seven[\[CapitalDelta]\[Phi]_,\[CapitalDelta]_][n_]:=If[n==0,-((2^(1-2 \[CapitalDelta]-2 \[CapitalDelta]\[Phi]) \[Pi] Gamma[1+2 \[CapitalDelta]] Gamma[\[CapitalDelta]\[Phi]]^3 Gamma[1/2+2 \[CapitalDelta]\[Phi]] (3+2 EulerGamma \[CapitalDelta]\[Phi]+4 \[CapitalDelta]\[Phi] PolyGamma[0,\[CapitalDelta]\[Phi]]-2 \[CapitalDelta]\[Phi] PolyGamma[0,2 \[CapitalDelta]\[Phi]]))/(Gamma[1/2 (-1+\[CapitalDelta])] Gamma[(1+\[CapitalDelta])/2]^3 Gamma[1/2+\[CapitalDelta]\[Phi]] Gamma[1/2-\[CapitalDelta]/2+\[CapitalDelta]\[Phi]]^2 Gamma[\[CapitalDelta]/2+\[CapitalDelta]\[Phi]]^2)),(2^(-6-2 n+\[CapitalDelta]-2 \[CapitalDelta]\[Phi]) Gamma[1/2+\[CapitalDelta]] Gamma[n+\[CapitalDelta]\[Phi]]^3 Gamma[-(1/2)+n+2 \[CapitalDelta]\[Phi]]^2 ((128\[Pi]  (6 n^2+3 n (-1+4 \[CapitalDelta]\[Phi])+\[CapitalDelta]\[Phi] (-1+4 \[CapitalDelta]\[Phi])) Gamma[\[CapitalDelta]/2+1] (-2 PolyGamma[0,n+\[CapitalDelta]\[Phi]]+PolyGamma[0,2 (n+\[CapitalDelta]\[Phi])]))/(Gamma[1/2+n+\[CapitalDelta]\[Phi]] Gamma[-(1/2)+2 n+2 \[CapitalDelta]\[Phi]])+(64^(n+\[CapitalDelta]\[Phi]) Gamma[1+\[CapitalDelta]/2] Gamma[-1+n+\[CapitalDelta]\[Phi]] (48 n^4+48 n^3 (-1+4 \[CapitalDelta]\[Phi])+2 \[CapitalDelta]\[Phi] (-1+4 \[CapitalDelta]\[Phi])^3+3 n (1-4 \[CapitalDelta]\[Phi])^2 (-1+8 \[CapitalDelta]\[Phi])+6 n^2 (3+4 \[CapitalDelta]\[Phi] (-7+16 \[CapitalDelta]\[Phi]))+2 n (-1+2 n+4 \[CapitalDelta]\[Phi]) (-1+4 n+4 \[CapitalDelta]\[Phi]) (6 n^2+3 n (-1+4 \[CapitalDelta]\[Phi])+\[CapitalDelta]\[Phi] (-1+4 \[CapitalDelta]\[Phi])) (PolyGamma[0,n]-PolyGamma[0,1/2+n+2 \[CapitalDelta]\[Phi]]+PolyGamma[0,1/2+2 n+2 \[CapitalDelta]\[Phi]])))/(n (-1+2 n+2 \[CapitalDelta]\[Phi]) (-1+2 n+4 \[CapitalDelta]\[Phi]) (-3+4 n+4 \[CapitalDelta]\[Phi]) (-1+4 n+4 \[CapitalDelta]\[Phi]) Gamma[4 (-1+n+\[CapitalDelta]\[Phi])])))/(\[Pi] (n!)^2 Gamma[1/2 (-1+\[CapitalDelta])] Gamma[(1+\[CapitalDelta])/2]^2 Gamma[1/2-\[CapitalDelta]/2+\[CapitalDelta]\[Phi]]^2 Gamma[\[CapitalDelta]/2+\[CapitalDelta]\[Phi]]^2)];
ROverA[\[CapitalDelta]\[Phi]_,\[CapitalDelta]_][n_]:=If[Mod[2 n,2]==0,-((2^(5-2 \[CapitalDelta]) Sqrt[\[Pi]] Gamma[2 \[CapitalDelta]] Gamma[n+\[CapitalDelta]\[Phi]]^4 Gamma[-(1/2)+n+2 \[CapitalDelta]\[Phi]]^2)/((n!)^2 Gamma[\[CapitalDelta]/2]^4 Gamma[2 (n+\[CapitalDelta]\[Phi])] Gamma[1/2 (-1+\[CapitalDelta])+\[CapitalDelta]\[Phi]]^2 Gamma[-(\[CapitalDelta]/2)+\[CapitalDelta]\[Phi]]^2 Gamma[-(1/2)+2 (n+\[CapitalDelta]\[Phi])])),0];
SOverA[\[CapitalDelta]\[Phi]_,\[CapitalDelta]_][n_]:=If[Mod[2 n,2]==0,(2^(5-2 \[CapitalDelta]) Sqrt[\[Pi]] Gamma[2 \[CapitalDelta]] Gamma[n+\[CapitalDelta]\[Phi]]^4 Gamma[-(1/2)+n+2 \[CapitalDelta]\[Phi]]^2 (HarmonicNumber[n]-2 HarmonicNumber[-1+n+\[CapitalDelta]\[Phi]]-HarmonicNumber[-(3/2)+n+2 \[CapitalDelta]\[Phi]]+2 HarmonicNumber[-2+4 n+4 \[CapitalDelta]\[Phi]]-Log[4]))/((n!)^2 Gamma[\[CapitalDelta]/2]^4 Gamma[2 (n+\[CapitalDelta]\[Phi])] Gamma[1/2 (-1+\[CapitalDelta])+\[CapitalDelta]\[Phi]]^2 Gamma[-(\[CapitalDelta]/2)+\[CapitalDelta]\[Phi]]^2 Gamma[-(1/2)+2 (n+\[CapitalDelta]\[Phi])]),0];

(*ROverA[\[CapitalDelta]\[Phi]_,0][n_]:=0;
SOverA[\[CapitalDelta]\[Phi]_,0][n_]:=0;*)

\[Beta]SS[\[CapitalDelta]\[Phi]_,\[CapitalDelta]_][n_]:=b0[2n]+2bCS[\[CapitalDelta]\[Phi],\[CapitalDelta], \[Beta]0s0OverA ][n]+bs0OverA[\[CapitalDelta]\[Phi],\[CapitalDelta]][n]; 
\[Alpha]SS[\[CapitalDelta]\[Phi]_,\[CapitalDelta]_][n_]:=a0[2n]+aCS[\[CapitalDelta]\[Phi],\[CapitalDelta], \[Beta]0s0OverA ][n]+as0OverA[\[CapitalDelta]\[Phi],\[CapitalDelta]][n];

(*\[Beta]SS[\[CapitalDelta]\[Phi]_,0][n_]:=0; 
\[Alpha]SS[\[CapitalDelta]\[Phi]_,0][n_]:=a0[2n];*)
\[Beta]pF[\[CapitalDelta]\[Phi]_,\[CapitalDelta]_][n_]:=bs1[\[CapitalDelta]\[Phi],\[CapitalDelta]][n]-2b1[2n+1];
\[Alpha]pF[\[CapitalDelta]\[Phi]_,\[CapitalDelta]_][n_]:=as1[\[CapitalDelta]\[Phi],\[CapitalDelta]][n]-2a1[2n+1];

agff[\[CapitalDelta]\[Phi]_,\[CapitalDelta]_]:=(2 Gamma[\[CapitalDelta]]^2 Gamma[\[CapitalDelta]+2 \[CapitalDelta]\[Phi]-1])/(Gamma[2 \[CapitalDelta]\[Phi]]^2 Gamma[2 \[CapitalDelta]-1] Gamma[\[CapitalDelta]-2 \[CapitalDelta]\[Phi]+1]);
defaultPrec=150;
Options[alphaminus]=Options[betaminus]=Options[alphaplus]=Options[betaplus]=Options[bosonicminusfuncList]=Options[fermionicplusfuncList]={PrecisionGoal->defaultPrec};

betaminus[ext_, order_, x_,OptionsPattern[]]:=
 Block[
{\[CapitalDelta]\[Phi] = SetPrecision[ext, OptionValue[PrecisionGoal]],
\[CapitalDelta] = SetPrecision[x, OptionValue[PrecisionGoal]], 
\[Beta]0s0OverA,
b0,
$MaxExtraPrecision
}, 
\[Beta]0s0OverA=\[Beta]0s0OverAfunc[\[CapitalDelta]\[Phi],\[CapitalDelta]];
If[Precision[x]===Infinity&&Precision[ext]===Infinity,
\[CapitalDelta]\[Phi] = ext;
\[CapitalDelta] = x;
$MaxExtraPrecision=Infinity;];
b0[-1] = 0; 
b0[0] = \[Beta]0s0OverA; 
b0[n_] := b0[n] = (ROverA[\[CapitalDelta]\[Phi],\[CapitalDelta]][(n - 1)/2] -(\[Rho][\[CapitalDelta]\[Phi],\[CapitalDelta]][n - 2]*b0[n - 2] + \[Nu][\[CapitalDelta]\[Phi],\[CapitalDelta]][n - 1]*b0[n - 1]))/\[Mu][n]; 
(*-b0[2order]-2bCS[\[CapitalDelta]\[Phi],\[CapitalDelta], \[Beta]0s0OverA ][order]-bs0OverA[\[CapitalDelta]\[Phi],\[CapitalDelta]][order]*)
-\[Beta]SS[\[CapitalDelta]\[Phi],\[CapitalDelta]][order]
];

betaminus[ext_, order_, 0, OptionsPattern[]]:=0;

betaminusSym[\[CapitalDelta]\[Phi]_, order_, \[CapitalDelta]_]:=
 Block[
{
\[Beta]0s0OverA,
b0,
$MaxExtraPrecision
}, 
\[Beta]0s0OverA=\[Beta]0s0OverAfunc[\[CapitalDelta]\[Phi],\[CapitalDelta]];
b0[-1] = 0; 
b0[0] = \[Beta]0s0OverA; 
b0[n_] := b0[n] = (ROverA[\[CapitalDelta]\[Phi],\[CapitalDelta]][(n - 1)/2] -(\[Rho][\[CapitalDelta]\[Phi],\[CapitalDelta]][n - 2]*b0[n - 2] + \[Nu][\[CapitalDelta]\[Phi],\[CapitalDelta]][n - 1]*b0[n - 1]))/\[Mu][n]; 
(*-b0[2order]-2bCS[\[CapitalDelta]\[Phi],\[CapitalDelta], \[Beta]0s0OverA ][order]-bs0OverA[\[CapitalDelta]\[Phi],\[CapitalDelta]][order]*)
-\[Beta]SS[\[CapitalDelta]\[Phi],\[CapitalDelta]][order]//FullSimplify
];

alphaminus[ext_, order_, x_, OptionsPattern[]] := Block[{
\[CapitalDelta]\[Phi] = SetPrecision[ext, 3/2*OptionValue[PrecisionGoal]],
\[CapitalDelta] = SetPrecision[x, 3/2*OptionValue[PrecisionGoal]], 
\[Beta]0s0OverA,
\[Alpha]0s0pre,
\[Alpha]0s0OverA,
prec=OptionValue[PrecisionGoal],
b0,
a0,
$MaxExtraPrecision
}, 
\[Alpha]0s0pre=If[ext==x&&Precision[x]===Infinity,
$MaxExtraPrecision=2prec;N[(\[Alpha]0s0prefunc[ext,x][ext  + 1/2*10^(-prec/2)] - \[Alpha]0s0prefunc[ext,x][ext - 1/2*10^(-prec/2)])/10^(-prec/2),prec],
(\[Alpha]0s0prefunc[\[CapitalDelta]\[Phi],\[CapitalDelta]][\[CapitalDelta]\[Phi]  + 1/2*10^(-prec/2)] - \[Alpha]0s0prefunc[\[CapitalDelta]\[Phi],\[CapitalDelta]][\[CapitalDelta]\[Phi]  - 1/2*10^(-prec/2)])/10^(-prec/2)];
\[CapitalDelta]\[Phi] = SetPrecision[ext, prec];
\[CapitalDelta] = SetPrecision[x, prec];
\[Beta]0s0OverA= \[Beta]0s0OverAfunc[\[CapitalDelta]\[Phi],\[CapitalDelta]];
\[Alpha]0s0OverA=-((2 Gamma[\[CapitalDelta]] Gamma[-(1/2)+2 \[CapitalDelta]\[Phi]])/(Gamma[\[CapitalDelta]/2]^4 Gamma[-(\[CapitalDelta]/2)+\[CapitalDelta]\[Phi]]^2))\[Alpha]0s0pre+EulerGamma \[Beta]0s0OverA;
If[Precision[x]===Infinity&&Precision[ext]===Infinity,
\[CapitalDelta]\[Phi] = ext;
\[CapitalDelta] = x;
$MaxExtraPrecision=Infinity;];
a0[0] =\[Alpha]0s0OverA; 
b0[-1] = 0; b0[0] = \[Beta]0s0OverA; 
b0[n_] := b0[n] = (ROverA[\[CapitalDelta]\[Phi],\[CapitalDelta]][(n - 1)/2] -(\[Rho][\[CapitalDelta]\[Phi],\[CapitalDelta]][n - 2]*b0[n - 2] + \[Nu][\[CapitalDelta]\[Phi],\[CapitalDelta]][n - 1]*b0[n - 1]))/\[Mu][n]; 
a0[-1] = 0; a0[0] =\[Alpha]0s0OverA; 
a0[n_] := a0[n] = (SOverA[\[CapitalDelta]\[Phi],\[CapitalDelta]][(n - 1)/2] - (\[Rho][\[CapitalDelta]\[Phi],\[CapitalDelta]][n-2] a0[n-2] +\[Nu][\[CapitalDelta]\[Phi],\[CapitalDelta]][n-1] a0[n-1]+\[Rho]p[\[CapitalDelta]\[Phi],\[CapitalDelta]][n-2] b0[n-2]+\[Nu]p[\[CapitalDelta]\[Phi],\[CapitalDelta]][n-1] b0[n-1]+\[Mu]p[n] b0[n]))/\[Mu][n]; 
-\[Alpha]SS[\[CapitalDelta]\[Phi],\[CapitalDelta]][order]//FullSimplify
];




alphaminus[ext_, order_, 0, OptionsPattern[]]:=Block[{
\[CapitalDelta]\[Phi] = SetPrecision[ext, OptionValue[PrecisionGoal]]
}, 
-agff[\[CapitalDelta]\[Phi],2\[CapitalDelta]\[Phi]+2#]&[order]
];

bosonicminusfuncList[ext_, order_, x_, OptionsPattern[]]:=Block[{
\[CapitalDelta]\[Phi] = SetPrecision[ext, 3/2*OptionValue[PrecisionGoal]],
\[CapitalDelta] = SetPrecision[x, 3/2*OptionValue[PrecisionGoal]], 
\[Beta]0s0OverA,
\[Alpha]0s0pre,
\[Alpha]0s0OverA,
prec=OptionValue[PrecisionGoal],
b0,
a0,
$MaxExtraPrecision
}, 
\[Alpha]0s0pre=If[ext==x&&Precision[x]===Infinity,
$MaxExtraPrecision=2prec;N[(\[Alpha]0s0prefunc[ext,x][ext  + 1/2*10^(-prec/2)] - \[Alpha]0s0prefunc[ext,x][ext - 1/2*10^(-prec/2)])/10^(-prec/2),prec],
(\[Alpha]0s0prefunc[\[CapitalDelta]\[Phi],\[CapitalDelta]][\[CapitalDelta]\[Phi]  + 1/2*10^(-prec/2)] - \[Alpha]0s0prefunc[\[CapitalDelta]\[Phi],\[CapitalDelta]][\[CapitalDelta]\[Phi]  - 1/2*10^(-prec/2)])/10^(-prec/2)];
\[CapitalDelta]\[Phi] = SetPrecision[ext, prec];
\[CapitalDelta] = SetPrecision[x, prec];
\[Beta]0s0OverA= \[Beta]0s0OverAfunc[\[CapitalDelta]\[Phi],\[CapitalDelta]];
\[Alpha]0s0OverA=-((2 Gamma[\[CapitalDelta]] Gamma[-(1/2)+2 \[CapitalDelta]\[Phi]])/(Gamma[\[CapitalDelta]/2]^4 Gamma[-(\[CapitalDelta]/2)+\[CapitalDelta]\[Phi]]^2))\[Alpha]0s0pre+EulerGamma \[Beta]0s0OverA;
If[Precision[x]===Infinity&&Precision[ext]===Infinity,
\[CapitalDelta]\[Phi] = ext;
\[CapitalDelta] = x;
$MaxExtraPrecision=Infinity;];
a0[0] =\[Alpha]0s0OverA; 
b0[-1] = 0; b0[0] = \[Beta]0s0OverA; 
b0[n_] := b0[n] = (ROverA[\[CapitalDelta]\[Phi],\[CapitalDelta]][(n - 1)/2] -(\[Rho][\[CapitalDelta]\[Phi],\[CapitalDelta]][n - 2]*b0[n - 2] + \[Nu][\[CapitalDelta]\[Phi],\[CapitalDelta]][n - 1]*b0[n - 1]))/\[Mu][n]; 
a0[-1] = 0; a0[0] =\[Alpha]0s0OverA; 
a0[n_] := a0[n] = (SOverA[\[CapitalDelta]\[Phi],\[CapitalDelta]][(n - 1)/2] - (\[Rho][\[CapitalDelta]\[Phi],\[CapitalDelta]][n-2] a0[n-2] +\[Nu][\[CapitalDelta]\[Phi],\[CapitalDelta]][n-1] a0[n-1]+\[Rho]p[\[CapitalDelta]\[Phi],\[CapitalDelta]][n-2] b0[n-2]+\[Nu]p[\[CapitalDelta]\[Phi],\[CapitalDelta]][n-1] b0[n-1]+\[Mu]p[n] b0[n]))/\[Mu][n]; 
{{0}~Join~(-\[Beta]SS[\[CapitalDelta]\[Phi],\[CapitalDelta]]/@Range[order]),-\[Alpha]SS[\[CapitalDelta]\[Phi],\[CapitalDelta]]/@Range[0,order]}//FullSimplify
];

bosonicminusfuncList[ext_, order_, 0, OptionsPattern[]]:=Block[{
\[CapitalDelta]\[Phi] = SetPrecision[ext, OptionValue[PrecisionGoal]]
}, 
{Table[0,{order+1}],-agff[\[CapitalDelta]\[Phi],2\[CapitalDelta]\[Phi]+2#]&/@Range[0,order]}
];

betatildezero[\[CapitalDelta]\[Phi]_,\[CapitalDelta]_]:=(2^(15-2 \[CapitalDelta]-2 \[CapitalDelta]\[Phi]) \[Pi] (-(((-1+\[CapitalDelta]) (1+\[CapitalDelta])^3 (-3+2 \[CapitalDelta]) (4+\[CapitalDelta]-2 \[CapitalDelta]\[Phi]) (-9+2 \[CapitalDelta] (3+\[CapitalDelta])+4 (5-2 \[CapitalDelta]\[Phi]) \[CapitalDelta]\[Phi]))/((-2+\[CapitalDelta])^3 \[CapitalDelta] (1+2 \[CapitalDelta]) (-5+\[CapitalDelta]+2 \[CapitalDelta]\[Phi])^2 (-3+\[CapitalDelta]+2 \[CapitalDelta]\[Phi])^2 (-1+\[CapitalDelta]+2 \[CapitalDelta]\[Phi])))+(-2+1/((\[CapitalDelta]-2 \[CapitalDelta]\[Phi]) (-5+\[CapitalDelta]+2 \[CapitalDelta]\[Phi])))/(2+\[CapitalDelta]-2 \[CapitalDelta]\[Phi])^2+(2 (-1+\[CapitalDelta]) (-1+2 \[CapitalDelta]) (-1+(-1+\[CapitalDelta]) \[CapitalDelta]) (2+1/(6+\[CapitalDelta]-\[CapitalDelta]^2-10 \[CapitalDelta]\[Phi]+4 \[CapitalDelta]\[Phi]^2)))/((-2+\[CapitalDelta])^3 (1+2 \[CapitalDelta]) (-5+\[CapitalDelta]+2 \[CapitalDelta]\[Phi])^2)) Gamma[-4+2 \[CapitalDelta]] Gamma[\[CapitalDelta]\[Phi]]^3 Gamma[-(1/2)+2 \[CapitalDelta]\[Phi]])/((-1+\[CapitalDelta]\[Phi])^3 (-5+4 \[CapitalDelta]\[Phi]) (-3+4 \[CapitalDelta]\[Phi]) Gamma[-1+\[CapitalDelta]/2]^4 Gamma[-(1/2)+\[CapitalDelta]\[Phi]] Gamma[1/2 (-5+\[CapitalDelta])+\[CapitalDelta]\[Phi]]^2 Gamma[-1-\[CapitalDelta]/2+\[CapitalDelta]\[Phi]]^2)+1/(\[Pi]^(3/2) Gamma[-1+\[CapitalDelta]/2] Gamma[-(\[CapitalDelta]/2)+\[CapitalDelta]\[Phi]]^2) 2^(-7+2 \[CapitalDelta]+2 \[CapitalDelta]\[Phi]) (9+2 (-5+\[CapitalDelta]) \[CapitalDelta]+2 \[CapitalDelta]\[Phi]) Gamma[-(3/2)+\[CapitalDelta]] Gamma[1/2 (-1+\[CapitalDelta])] Gamma[-1+\[CapitalDelta]\[Phi]]^4 Gamma[-2+\[CapitalDelta]/2+\[CapitalDelta]\[Phi]] Gamma[-(5/2)+2 \[CapitalDelta]\[Phi]] HypergeometricPFQRegularized[{-(1/2)+\[CapitalDelta]/2,-(1/2)+\[CapitalDelta]/2,-(5/2)+\[CapitalDelta]/2+\[CapitalDelta]\[Phi],-(5/2)+2 \[CapitalDelta]\[Phi],-5+\[CapitalDelta]+2 \[CapitalDelta]\[Phi]},{-(3/2)+\[CapitalDelta],-(3/2)+\[CapitalDelta]/2+\[CapitalDelta]\[Phi],-(7/2)+\[CapitalDelta]/2+2 \[CapitalDelta]\[Phi],-(7/2)+\[CapitalDelta]/2+2 \[CapitalDelta]\[Phi]},1]-1/(Sqrt[\[Pi]] (1+2 \[CapitalDelta]) Gamma[\[CapitalDelta]/2]^2 Gamma[-1-\[CapitalDelta]/2+\[CapitalDelta]\[Phi]]^2) 2^(-1-\[CapitalDelta]+2 \[CapitalDelta]\[Phi]) (-1+2 \[CapitalDelta]) (1+(-2+\[CapitalDelta]) \[CapitalDelta]^2) (-3+2 (-1+\[CapitalDelta]) \[CapitalDelta]+2 \[CapitalDelta]\[Phi]) Gamma[-3+2 \[CapitalDelta]] Gamma[-1+\[CapitalDelta]\[Phi]]^4 Gamma[-1+\[CapitalDelta]/2+\[CapitalDelta]\[Phi]] Gamma[-(5/2)+2 \[CapitalDelta]\[Phi]] HypergeometricPFQRegularized[{1/2+\[CapitalDelta]/2,1/2+\[CapitalDelta]/2,-(3/2)+\[CapitalDelta]/2+\[CapitalDelta]\[Phi],-(5/2)+2 \[CapitalDelta]\[Phi],-3+\[CapitalDelta]+2 \[CapitalDelta]\[Phi]},{1/2+\[CapitalDelta],-(1/2)+\[CapitalDelta]/2+\[CapitalDelta]\[Phi],-(5/2)+\[CapitalDelta]/2+2 \[CapitalDelta]\[Phi],-(5/2)+\[CapitalDelta]/2+2 \[CapitalDelta]\[Phi]},1]+1/(Sqrt[\[Pi]] (-1+4 \[CapitalDelta]^2) Gamma[\[CapitalDelta]/2]^2 Gamma[-2-\[CapitalDelta]/2+\[CapitalDelta]\[Phi]]^2) 2^(-4-\[CapitalDelta]+2 \[CapitalDelta]\[Phi]) (1+\[CapitalDelta])^3 (1+2 \[CapitalDelta] (3+\[CapitalDelta])+2 \[CapitalDelta]\[Phi]) Gamma[1+2 \[CapitalDelta]] Gamma[-1+\[CapitalDelta]\[Phi]]^4 Gamma[\[CapitalDelta]/2+\[CapitalDelta]\[Phi]] Gamma[-(5/2)+2 \[CapitalDelta]\[Phi]] HypergeometricPFQRegularized[{3/2+\[CapitalDelta]/2,3/2+\[CapitalDelta]/2,-(1/2)+\[CapitalDelta]/2+\[CapitalDelta]\[Phi],-(5/2)+2 \[CapitalDelta]\[Phi],-1+\[CapitalDelta]+2 \[CapitalDelta]\[Phi]},{5/2+\[CapitalDelta],1/2+\[CapitalDelta]/2+\[CapitalDelta]\[Phi],-(3/2)+\[CapitalDelta]/2+2 \[CapitalDelta]\[Phi],-(3/2)+\[CapitalDelta]/2+2 \[CapitalDelta]\[Phi]},1];

betaplus[ext_, order_, x_,OptionsPattern[]]:=Block[{\[CapitalDelta]\[Phi] = SetPrecision[ext, OptionValue[PrecisionGoal]],
\[CapitalDelta] = SetPrecision[x, OptionValue[PrecisionGoal]], 
\[Beta]01,
b1,
$MaxExtraPrecision
}, 
\[Beta]01=\[Beta]01func[\[CapitalDelta]\[Phi],\[CapitalDelta]];
If[Precision[x]===Infinity&&Precision[ext]===Infinity,
\[CapitalDelta]\[Phi] = ext;
\[CapitalDelta] = x;
$MaxExtraPrecision=Infinity;];
b1[-1] = 0; 
b1[0] =\[Beta]01; 
b1[n_?OddQ]:= b1[n] = (Reven[\[CapitalDelta]\[Phi],\[CapitalDelta]][(n - 1)/2] -(\[Rho][\[CapitalDelta]\[Phi],\[CapitalDelta]][n - 2]*b1[n - 2] + \[Nu][\[CapitalDelta]\[Phi],\[CapitalDelta]][n - 1]*b1[n - 1]))/\[Mu][n]; 
b1[n_?EvenQ]:= b1[n] = (Rodd[\[CapitalDelta]\[Phi],\[CapitalDelta]][(n - 2)/2] -(\[Rho][\[CapitalDelta]\[Phi],\[CapitalDelta]][n - 2]*b1[n - 2] + \[Nu][\[CapitalDelta]\[Phi],\[CapitalDelta]][n - 1]*b1[n - 1]))/\[Mu][n]; 
-\[Beta]pF[\[CapitalDelta]\[Phi],\[CapitalDelta]][order]//FullSimplify
];

alphaplus[ext_, order_, x_, OptionsPattern[]] := Block[{
\[CapitalDelta]\[Phi] = SetPrecision[ext, 3/2*OptionValue[PrecisionGoal]],
\[CapitalDelta] = SetPrecision[x, 3/2*OptionValue[PrecisionGoal]], 
\[Beta]01,
\[Alpha]01pre,
\[Alpha]01,
prec=OptionValue[PrecisionGoal],
b1,
a1,
$MaxExtraPrecision
}, 
\[Alpha]01pre=If[x==1,EulerGamma,If[(2ext+x==1)&&Precision[x]===Infinity,
$MaxExtraPrecision=2prec;N[(\[Alpha]01prefunc[ext,x][ext  + 1/2*10^(-prec/2)] - \[Alpha]01prefunc[ext,x][ext - 1/2*10^(-prec/2)])/10^(-prec/2),prec],
(\[Alpha]01prefunc[\[CapitalDelta]\[Phi],\[CapitalDelta]][\[CapitalDelta]\[Phi]  + 1/2*10^(-prec/2)] - \[Alpha]01prefunc[\[CapitalDelta]\[Phi],\[CapitalDelta]][\[CapitalDelta]\[Phi]  - 1/2*10^(-prec/2)])/10^(-prec/2)]];
(*\[Alpha]01pre=(\[Alpha]01prefunc[\[CapitalDelta]\[Phi],\[CapitalDelta]][\[CapitalDelta]\[Phi]  + 1/2*10^(-prec/2)] - \[Alpha]01prefunc[\[CapitalDelta]\[Phi],\[CapitalDelta]][\[CapitalDelta]\[Phi]  - 1/2*10^(-prec/2)])/10^(-prec/2);*)
\[CapitalDelta]\[Phi] = SetPrecision[ext, prec];
\[CapitalDelta] = SetPrecision[x, prec];
\[Beta]01=\[Beta]01func[\[CapitalDelta]\[Phi],\[CapitalDelta]];
\[Alpha]01=\[Alpha]01pre+EulerGamma \[Beta]01+\[Alpha]01funcAdd[\[CapitalDelta]\[Phi],\[CapitalDelta]];
If[Precision[x]===Infinity&&Precision[ext]===Infinity,
\[CapitalDelta]\[Phi] = ext;
\[CapitalDelta] = x;
$MaxExtraPrecision=Infinity;];
b1[-1] = 0; 
b1[0] =\[Beta]01; 
b1[n_?OddQ]:= b1[n] = (Reven[\[CapitalDelta]\[Phi],\[CapitalDelta]][(n - 1)/2] -(\[Rho][\[CapitalDelta]\[Phi],\[CapitalDelta]][n - 2]*b1[n - 2] + \[Nu][\[CapitalDelta]\[Phi],\[CapitalDelta]][n - 1]*b1[n - 1]))/\[Mu][n]; 
b1[n_?EvenQ]:= b1[n] = (Rodd[\[CapitalDelta]\[Phi],\[CapitalDelta]][(n - 2)/2] -(\[Rho][\[CapitalDelta]\[Phi],\[CapitalDelta]][n - 2]*b1[n - 2] + \[Nu][\[CapitalDelta]\[Phi],\[CapitalDelta]][n - 1]*b1[n - 1]))/\[Mu][n]; 
a1[-1] = 0; 
a1[0] =\[Alpha]01; 
a1[n_?OddQ]:= a1[n] = (Seven[\[CapitalDelta]\[Phi],\[CapitalDelta]][(n - 1)/2] -(\[Rho][\[CapitalDelta]\[Phi],\[CapitalDelta]][n - 2]*a1[n - 2] + \[Nu][\[CapitalDelta]\[Phi],\[CapitalDelta]][n - 1]*a1[n - 1]+\[Rho]p[\[CapitalDelta]\[Phi],\[CapitalDelta]][n-2] b1[n-2]+\[Nu]p[\[CapitalDelta]\[Phi],\[CapitalDelta]][n-1] b1[n-1]+\[Mu]p[n] b1[n]))/\[Mu][n]; 
a1[n_?EvenQ]:= a1[n] = (Sodd[\[CapitalDelta]\[Phi],\[CapitalDelta]][(n - 2)/2] -(\[Rho][\[CapitalDelta]\[Phi],\[CapitalDelta]][n - 2]*a1[n - 2] + \[Nu][\[CapitalDelta]\[Phi],\[CapitalDelta]][n - 1]*a1[n - 1]+\[Rho]p[\[CapitalDelta]\[Phi],\[CapitalDelta]][n-2] b1[n-2]+\[Nu]p[\[CapitalDelta]\[Phi],\[CapitalDelta]][n-1] b1[n-1]+\[Mu]p[n] b1[n]))/\[Mu][n]; 
-\[Alpha]pF[\[CapitalDelta]\[Phi],\[CapitalDelta]][order]//FullSimplify
];

fermionicplusfuncList[ext_, order_, x_, OptionsPattern[]]:= Block[{
\[CapitalDelta]\[Phi] = SetPrecision[ext, 3/2*OptionValue[PrecisionGoal]],
\[CapitalDelta] = SetPrecision[x, 3/2*OptionValue[PrecisionGoal]], 
\[Beta]01,
\[Alpha]01pre,
\[Alpha]01,
prec=OptionValue[PrecisionGoal],
b1,
a1,
$MaxExtraPrecision
}, 

\[Alpha]01pre=If[x==1,EulerGamma,If[(2ext+x==1)&&Precision[x]===Infinity,
$MaxExtraPrecision=2prec;N[(\[Alpha]01prefunc[ext,x][ext  + 1/2*10^(-prec/2)] - \[Alpha]01prefunc[ext,x][ext - 1/2*10^(-prec/2)])/10^(-prec/2),prec],
(\[Alpha]01prefunc[\[CapitalDelta]\[Phi],\[CapitalDelta]][\[CapitalDelta]\[Phi]  + 1/2*10^(-prec/2)] - \[Alpha]01prefunc[\[CapitalDelta]\[Phi],\[CapitalDelta]][\[CapitalDelta]\[Phi]  - 1/2*10^(-prec/2)])/10^(-prec/2)]];

\[CapitalDelta]\[Phi] = SetPrecision[ext, prec];
\[CapitalDelta] = SetPrecision[x, prec];
\[Beta]01=\[Beta]01func[\[CapitalDelta]\[Phi],\[CapitalDelta]];
\[Alpha]01=\[Alpha]01pre+EulerGamma \[Beta]01+\[Alpha]01funcAdd[\[CapitalDelta]\[Phi],\[CapitalDelta]];

If[Precision[x]===Infinity&&Precision[ext]===Infinity,
\[CapitalDelta]\[Phi] = ext;
\[CapitalDelta] = x;
$MaxExtraPrecision=Infinity;];
b1[-1] = 0; 
b1[0] =\[Beta]01; 
b1[n_?OddQ]:= b1[n] = (Reven[\[CapitalDelta]\[Phi],\[CapitalDelta]][(n - 1)/2] -(\[Rho][\[CapitalDelta]\[Phi],\[CapitalDelta]][n - 2]*b1[n - 2] + \[Nu][\[CapitalDelta]\[Phi],\[CapitalDelta]][n - 1]*b1[n - 1]))/\[Mu][n]; 
b1[n_?EvenQ]:= b1[n] = (Rodd[\[CapitalDelta]\[Phi],\[CapitalDelta]][(n - 2)/2] -(\[Rho][\[CapitalDelta]\[Phi],\[CapitalDelta]][n - 2]*b1[n - 2] + \[Nu][\[CapitalDelta]\[Phi],\[CapitalDelta]][n - 1]*b1[n - 1]))/\[Mu][n]; 
a1[-1] = 0; 
a1[0] =\[Alpha]01; 
a1[n_?OddQ]:= a1[n] = (Seven[\[CapitalDelta]\[Phi],\[CapitalDelta]][(n - 1)/2] -(\[Rho][\[CapitalDelta]\[Phi],\[CapitalDelta]][n - 2]*a1[n - 2] + \[Nu][\[CapitalDelta]\[Phi],\[CapitalDelta]][n - 1]*a1[n - 1]+\[Rho]p[\[CapitalDelta]\[Phi],\[CapitalDelta]][n-2] b1[n-2]+\[Nu]p[\[CapitalDelta]\[Phi],\[CapitalDelta]][n-1] b1[n-1]+\[Mu]p[n] b1[n]))/\[Mu][n]; 
a1[n_?EvenQ]:= a1[n] = (Sodd[\[CapitalDelta]\[Phi],\[CapitalDelta]][(n - 2)/2] -(\[Rho][\[CapitalDelta]\[Phi],\[CapitalDelta]][n - 2]*a1[n - 2] + \[Nu][\[CapitalDelta]\[Phi],\[CapitalDelta]][n - 1]*a1[n - 1]+\[Rho]p[\[CapitalDelta]\[Phi],\[CapitalDelta]][n-2] b1[n-2]+\[Nu]p[\[CapitalDelta]\[Phi],\[CapitalDelta]][n-1] b1[n-1]+\[Mu]p[n] b1[n]))/\[Mu][n]; 
{-\[Beta]pF[\[CapitalDelta]\[Phi],\[CapitalDelta]]/@Range[0,order],-\[Alpha]pF[\[CapitalDelta]\[Phi],\[CapitalDelta]]/@Range[0,order]}//FullSimplify
];

fermionicplusfuncList[ext_, order_, 0, OptionsPattern[]]:= Block[{
\[CapitalDelta]\[Phi] = SetPrecision[ext, 3/2*OptionValue[PrecisionGoal]],
\[CapitalDelta] =0, 
\[Beta]01,
\[Alpha]01pre,
\[Alpha]01,
prec=OptionValue[PrecisionGoal],
b1,
a1,
$MaxExtraPrecision,
x=0
}, 

\[Alpha]01pre=(\[Alpha]01prefunc[\[CapitalDelta]\[Phi],\[CapitalDelta]][\[CapitalDelta]\[Phi]  + 1/2*10^(-prec/2)] - \[Alpha]01prefunc[\[CapitalDelta]\[Phi],\[CapitalDelta]][\[CapitalDelta]\[Phi]  - 1/2*10^(-prec/2)])/10^(-prec/2);

\[CapitalDelta]\[Phi] = SetPrecision[ext, prec];
\[CapitalDelta] = SetPrecision[x, prec];

\[Beta]01=SetPrecision[\[Beta]01func[ext,\[CapitalDelta]], 3/2*OptionValue[PrecisionGoal]];

\[Alpha]01=\[Alpha]01pre+EulerGamma \[Beta]01+\[Alpha]01funcAdd[\[CapitalDelta]\[Phi],\[CapitalDelta]];

If[Precision[x]===Infinity&&Precision[ext]===Infinity,
\[CapitalDelta]\[Phi] = ext;
\[CapitalDelta] = x;
$MaxExtraPrecision=Infinity;];
b1[-1] = 0; 
b1[0] =\[Beta]01; 
b1[n_?OddQ]:= b1[n] = (Reven[\[CapitalDelta]\[Phi],\[CapitalDelta]][(n - 1)/2] -(\[Rho][\[CapitalDelta]\[Phi],\[CapitalDelta]][n - 2]*b1[n - 2] + \[Nu][\[CapitalDelta]\[Phi],\[CapitalDelta]][n - 1]*b1[n - 1]))/\[Mu][n]; 
b1[n_?EvenQ]:= b1[n] = (Rodd[\[CapitalDelta]\[Phi],\[CapitalDelta]][(n - 2)/2] -(\[Rho][\[CapitalDelta]\[Phi],\[CapitalDelta]][n - 2]*b1[n - 2] + \[Nu][\[CapitalDelta]\[Phi],\[CapitalDelta]][n - 1]*b1[n - 1]))/\[Mu][n]; 
a1[-1] = 0; 
a1[0] =\[Alpha]01; 
a1[n_?OddQ]:= a1[n] = (Seven[\[CapitalDelta]\[Phi],\[CapitalDelta]][(n - 1)/2] -(\[Rho][\[CapitalDelta]\[Phi],\[CapitalDelta]][n - 2]*a1[n - 2] + \[Nu][\[CapitalDelta]\[Phi],\[CapitalDelta]][n - 1]*a1[n - 1]+\[Rho]p[\[CapitalDelta]\[Phi],\[CapitalDelta]][n-2] b1[n-2]+\[Nu]p[\[CapitalDelta]\[Phi],\[CapitalDelta]][n-1] b1[n-1]+\[Mu]p[n] b1[n]))/\[Mu][n]; 
a1[n_?EvenQ]:= a1[n] = (Sodd[\[CapitalDelta]\[Phi],\[CapitalDelta]][(n - 2)/2] -(\[Rho][\[CapitalDelta]\[Phi],\[CapitalDelta]][n - 2]*a1[n - 2] + \[Nu][\[CapitalDelta]\[Phi],\[CapitalDelta]][n - 1]*a1[n - 1]+\[Rho]p[\[CapitalDelta]\[Phi],\[CapitalDelta]][n-2] b1[n-2]+\[Nu]p[\[CapitalDelta]\[Phi],\[CapitalDelta]][n-1] b1[n-1]+\[Mu]p[n] b1[n]))/\[Mu][n]; 
{-\[Beta]pF[\[CapitalDelta]\[Phi],\[CapitalDelta]]/@Range[0,order],-\[Alpha]pF[\[CapitalDelta]\[Phi],\[CapitalDelta]]/@Range[0,order]}//FullSimplify
];

functionalpair[ext_, order_, x_, opts___]:={bosonicminusfuncList[ext, order, x, opts],fermionicplusfuncList[ext, order, x, opts]};

dn[a_][n_]:=(Pochhammer[a,n]^4*Pochhammer[4a-1,2n])/(n!^2Pochhammer[2a,n]^2Pochhammer[4a+2n-1,2n]);
cn[a_][n_]:=-((2^(-3+6*a+4*n)*Gamma[1/2+a]*Gamma[a+n]^4*Gamma[-1/2+2*a+n]^2*(HarmonicNumber[n]-2*HarmonicNumber[-1+a+n]+HarmonicNumber[-1+2*a+n]-2*HarmonicNumber[-2+4*a+2*n]+2*HarmonicNumber[-2+4*a+4*n]))/(Pi*n!^2*Gamma[a]^3*Gamma[-1/2+2*a]*Gamma[-1+4*a+4*n]));




End[]
EndPackage[]
