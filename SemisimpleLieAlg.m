(* ::Package:: *)

BeginPackage["SemisimpleLieAlg`"];


Unprotect[ZeroMatrix,StandardBasis,CartanMatrix,CartanMatrixQ,ConnectedComponents,
		DynkinDiagram,SimpleRootLengths,ClearRootCaches,RootInnerProduct,PositiveRoots];


(* ::Subsubtitle::Closed:: *)
(*Sources*)


(* ::Text:: *)
(*Cahn*)
(*Wybourne*)
(*O'Raifeartaigh*)
(*Knapp*)
(*Inspiration for Cartan matrix labeling from Sage*)


(* ::Subsubtitle::Closed:: *)
(*Notes:*)


(* ::Text:: *)
(*The caching system for SimpleRootLengths and PositiveRoots is clumsy and maybe unnesc.  *)
(*PositiveRoots[ , ,"Method"->"ALG"] is really slow without caching SimpleRoots BUT maybe I should just always use the default "Method"->"SRS".*)


(* ::Subsubtitle:: *)
(*Usages*)


(* ::Text:: *)
(*All public functions must be given a usage \[LongDash] else they won't be public!*)


ZeroMatrix::usage = "ZeroMatrix[dim_] returns a dim\[Times]dim sparse array of zeros.";


StandardBasis::usage = "StandardBasis[{n1,n2,...},dim] returns a dim\[Times]dim sparse array with a 1 in the {n1,n2,...} position.";


CartanMatrix::usage = "CartanMatrix[type,rank] will return the Cartan matrix of the given type and rank.  If rank>10 then it returns a sparse array.";


CartanMatrixInverse::usage = "Gives the inverse of the CartanMatrix[type, rank]";


CartanType::usage = "CartanType[A].  If A is a simple Cartan matrix then CartanType returns its Cartan type.  If A is not simple, then CartanType returns a list of types along with their position in A.";


CartanMatrixQ::usage = "CartanMatrixQ[A] checks whether A is a valid Cartan matrix \[LongDash] ie: it has nonzero determinant; the diagonal elements are all 2; A_{ij}\[LessEqual]0 \[ForAll] i\[NotEqual]j; A_{ij}\[NotEqual]0 iff A_{ji}\[NotEqual]0; A_{ij}A_{ji}=0,1,2,3 for i\[NotEqual]j.";


ConnectedComponents::usage = "ConnectedComponents takes a graph as a list of edge rules or incidence matrix and returns the connected components.  In particular, we use this for breaking up Cartan matrices of semisimple Lie algebras.";


DynkinDiagram::usage = "DynkinDiagram is passed a Cartan matrix and optional labels and list of dotted edges.  The vertices are by default labeled from 1 to rnk.
An alternative usage is to pass it a Cartan label eg {\[OpenCurlyDoubleQuote]a\[CloseCurlyDoubleQuote],3} or \[OpenCurlyDoubleQuote]a3\[CloseCurlyDoubleQuote].  With the same optional extras. Labels can be turned off completely using the option \[OpenCurlyDoubleQuote]Labels\[CloseCurlyDoubleQuote]->False.
Extending the standard convention, Black dots correspond to the short roots and White to the long roots.  In the case where there is more than two root lengths (not a standard Dynkin diagram) each length is given a corresponding shade of grey.";


ClearRootCaches::usage = "ClearRootCaches[] clears all the stored SimpleRootLengths, RootInnerProduct and PositiveRoots calculations.";


SimpleRootLengths::usage = "Given a Cartan matrix of rank r, this returns a vector of length r containing the \!\(\*
StyleBox[\"squared\",\nFontWeight->\"Bold\"]\) lengths.
Options are \[OpenCurlyDoubleQuote]norm\[CloseCurlyDoubleQuote]->\[OpenCurlyDoubleQuote]Short\[CloseCurlyDoubleQuote],\[OpenCurlyDoubleQuote]Long\[CloseCurlyDoubleQuote] or None and \[OpenCurlyDoubleQuote]normTo\[CloseCurlyDoubleQuote]->norm.
SimpleRootLengths stores previously calculated results \[LongDash] they may be cleared using ClearSimpleRootLengths[].";


RootInnerProduct::usage = "RootInnerProduct[type_,simpleroots_][r1,r2] where type is a Cartan matrix or Cartan type, and r1 & r2 are linear combinations of simpleroots returns the innerproduct of the two roots. The normalisation options are the same as in SimpleRootLengths.";


PositiveRoots::usage = "PositiveRoots[type_,simpleroots_] where type is a Cartan matrix or Cartan type and the list of simpleroots is optional returns all the positive roots grouped according to their level.
PositiveRoots stores previously calculated results \[LongDash] they may be cleared using ClearPositiveRoots[].";


FundamentalWeights::usage = "FundamentalWeights[A\[Or]type] returns the fundamental weights in the SRS rep - ie it returns the inverse of the Cartan matrix.
FundamentalWeights[A\[Or]type,srs] expands out the fundamental weights in terms of the given symbols for the simple roots.";
PrimitiveWeights::usage = "An alias for FundamentalWeights.";


WeylReflection::usage = "WeylReflection[A\[Or]type,i_Integer] returns the reflection matrix that acts on the right of a weight in the DYN rep.";
SimpleReflections::usage = "SimpleReflections[A\[Or]type] = WeylReflection[A\[Or]type,#]&/@Range[Length[A]]";


WeylOrbit::usage = "WeylOrbit[A\[Or]type,\[Lambda]] returns the weights on the orbit of \[Lambda].";


WeylVector::usage = "WeylVector[A\[Or]type] returns the half-sum of all positive roots in the SRS rep. In the DYN rep it is {1,1,...,1}.";


WeightSystem::usage = "WeightSystem[A\[Or]type,M] returns the weights in the representation given by the Dynkin indices M. By default it has the option Multiplicity->False.";


WeightMults::usage = "WeightMults[A\[Or]type,WeightSystem] calculates the multiplicites for the supplied output of WeightSystem.";


SpindleForm::usage = "Displays weight systems in the traditional spindle form.  Apply to the output of of WeightSystem";


WeylDimension::usage = "WeylDimension[A\[Or]type,M] returns the dimension of the representation given by the Dynkin indices M.";


GraphRoots2D::usage = "Given a 2D Cartan matrix and optional simple root label and styles, this plots the root system.";
GraphRoots3D::usage = "Given a 3D Cartan matrix and optional simple root label, this plots the root system.";


GraphWeights2D::usage = "";
GraphWeights3D::usage = "";


LexPositive::usage = "";
LexOrderedQ::usage = "";
LexSort::usage = "";
SRSPositive::usage = "";
SRSOrderedQ::usage = "";
SRSSort::usage = "";
DynSort::usage = "";
DynOrderedQ::usage = "";


(* ::Subsubtitle::Closed:: *)
(*Begin Private*)


Begin["`Private`"]


(* ::Subsubtitle::Closed:: *)
(*Matrices*)


(* ::Subsection::Closed:: *)
(*Extend Position to SparseArrays*)


(* ::Text:: *)
(*Note that Position and Cases have the same optional arguments;  so*)


Unprotect[Position];
Position[A_SparseArray,pat_,others___]:=Cases[ArrayRules[A],HoldPattern[_->pat],others][[All,1]]
Protect[Position];


(* ::Subsection::Closed:: *)
(*Options, SparseSwitch, ZeroMatrix and StandardBasis*)


Options[ZeroMatrix]={"Sparse"->Automatic,"SparseAfter"->10};
Options[StandardBasis]={"Sparse"->Automatic,"SparseAfter"->10};
Options[CartanMatrix]={"Sparse"->Automatic,"SparseAfter"->10};


SparseSwitch::sparse = "not valid option for Sparse, must be either Automatic, True or False.";
SparseSwitch::after = "not valid option for SparseAfter, must be positive integer.";


SparseSwitch[rank_,sparse_,after_]:=(If[IntegerQ[after]&&after>0, Null, Message[SparseSwitch::after];Return[$Failed]];
	Switch[sparse,
		Automatic,If[rank<after,Normal@#,#],
		True,#,
		False,Normal@#,
		_,Message[SparseSwitch::sparse];Return[$Failed]
]&)


ZeroMatrix[n_Integer?Positive,opts:OptionsPattern[]]:=
	SparseSwitch[n,OptionValue[ZeroMatrix,"Sparse"],OptionValue[ZeroMatrix,"SparseAfter"]]@SparseArray[{},{n,n}]


StandardBasis[list:{_Integer?Positive..},dim_Integer?Positive,opts:OptionsPattern[]]/;dim>=Max[list]:=
	SparseSwitch[dim,OptionValue[StandardBasis,"Sparse"],OptionValue[StandardBasis,"SparseAfter"]]@SparseArray[list->1,dim]


(* ::Subsection::Closed:: *)
(*CartanMatrix and CartanMatrixStringSplit*)


(* ::Text:: *)
(*Low dimensional isomorphisms...*)
(*b1\[Congruent]c1\[Congruent]a1,	c2=b2\[Transpose]\[TildeFullEqual]b2,	d2\[Congruent]a1*a1,	d3\[TildeFullEqual]a3,		e4\[Congruent]a4,		e5\[Congruent]d5*)


CartanMatrix::type1 = "The type `1` is not supported";
CartanMatrix::type2 = "`1` in `2` is not a valid {type,rank} pair";
CartanMatrix::type3 = "The substring(s) `1` in `2` can not be parsed as valid {type,rank} pairs";
CartanMatrix::rank  = "The rank `1` is not valid for the type `2`";


CartanMatrix[type_String,rank_Integer?Positive,OptionsPattern[]]:=With[{spQ=OptionValue["Sparse"],
	aSeq=Sequence[Band[{1,1}]->2,Band[{2,1}]->-1,Band[{1,2}]->-1]},
	SparseSwitch[rank,OptionValue[CartanMatrix,"Sparse"],OptionValue[CartanMatrix,"SparseAfter"]][
	If[rank==1\[And]MemberQ[{"a","b","c"},ToLowerCase[type]],SparseArray[{{2}}],
	Switch[ToLowerCase[type],
		"a",SparseArray[{aSeq},{rank,rank}],
		"b",SparseArray[{{rank-1,rank}->-2,aSeq},{rank,rank}],
		"c",SparseArray[{{rank,rank-1}->-2,aSeq},{rank,rank}],
		"d",If[rank==2,{{2,0},{0,2}},
				If[rank>2,SparseArray[{{rank-1,rank}->0,{rank,rank-1}->0,{rank-2,rank}->-1,{rank,rank-2}->-1,aSeq},{rank,rank}],
			Message[CartanMatrix::rank,rank,type];Return[HoldForm[CartanMatrix[type,rank]]]]],
		"e",If[rank>3,SparseArray[{{3,rank}->-1,{rank,3}->-1,{rank-1,rank}->0,{rank,rank-1}->0,aSeq},{rank,rank}],
			Message[CartanMatrix::rank,rank,type];Return[HoldForm[CartanMatrix[type,rank]]]],
		"f",If[rank==4,SparseArray[{{2,-1,0,0},{-1,2,-2,0},{0,-1,2,-1},{0,0,-1,2}}],
			Message[CartanMatrix::rank,rank,type];Return[HoldForm[CartanMatrix[type,rank]]]],
		"g",If[rank==2,SparseArray[{{2,-3},{-1,2}}],
			Message[CartanMatrix::rank,rank,type];Return[HoldForm[CartanMatrix[type,rank]]]],
		_,Message[CartanMatrix::type1,type];Return[$Failed]]]]];


CartanMatrix[{type_String,rank_Integer?Positive},opts:OptionsPattern[]]:=CartanMatrix[type,rank,opts]


CartanMatrix[blocks:{{_String,_Integer?Positive}..},opts:OptionsPattern[]]:=Module[{x},
	ArrayFlatten[DiagonalMatrix[Array[x,Length@blocks]]/.x[i_]:>CartanMatrix[Sequence@@blocks[[i]],opts]]]


(* ::Text:: *)
(*(* Version 1 *)*)
(*CartanMatrixStringSplit[type_String] := Catch[Module[{str = StringTrim@ToLowerCase@type},*)
(*      	str = StringReplace[str, Whitespace -> ""];*)
(*      	str = StringSplit[str, {"*", "x"}];*)
(*      	str = Table[*)
(*          		If[StringMatchQ[term, CharacterRange["a", "g"] ~~ NumberString],*)
(*            			{StringTake[term, 1], ToExpression@StringDrop[term, 1]},*)
(*            			Message[CartanMatrix::type2, term, type]; Throw[$Failed]], {term, str}];*)
(*      	Return[str]*)
(*      ]]*)


(* ::Text:: *)
(*CartanMatrixStringSplit -- version 2:*)
(*No longer requires the separators x or * between factors.  Is this new method too flexible?*)
(*Maybe replace CharacterRange["a", "g"] with LetterCharacter to make more general.  The check on the range occurs in CartanMatrix.*)


CartanMatrixStringSplit[type_String]:=Catch[Module[{
  splits={Whitespace,"*","x"}, str=StringTrim@ToLowerCase@type, str2},
	str=DeleteCases[StringSplit[str,splits],""];(*Print@str;*)
	str2=Flatten[StringCases[#,CharacterRange["a","g"]~~NumberString]&/@str];(*Print@str2;*)
	str=DeleteCases[StringSplit[StringJoin[str],Reverse@SortBy[str2,StringLength]],""];(*Print@str;*)
	If[Length[Flatten[{str}]]>0,Message[CartanMatrix::type3,str,type];Throw[$Failed]];
	str=Table[{StringTake[term,1],ToExpression@StringDrop[term,1]},{term,str2}];
	Return[str]
]]


CartanMatrix[type_String,opts:OptionsPattern[]]/;StringLength[type]>1:=CartanMatrix[CartanMatrixStringSplit[type],opts]
CartanMatrix[type:{_String..},opts:OptionsPattern[]]:=CartanMatrix[StringJoin[type],opts]


CartanBlocksQ[expr_]:= MatchQ[expr,(_String|{_String,_Integer}|{{_String,_Integer}..})]


(* ::Subsection::Closed:: *)
(*CartanMatrixQ*)


Options[CartanMatrixQ] = {"Strict"->False};


CartanMatrixQ[A_, why:(True|False):False,opts:OptionsPattern[]]:=
	If[OptionValue[CartanMatrixQ,"Strict"],CartanMatrixQStrict[A,why],CartanMatrixQLax[A,why]]


(* ::Text:: *)
(*The "strict" method is a bit old... needs cleaning*)


CartanMatrixQ::notCartan = "The matrix is not Cartan because `1`";


CartanMatrixQStrict[A_?MatrixQ, why_] := Catch[Module[{x, dim = Length[A], PrintWhy}, 
   	PrintWhy[txt_] := If[why, Message[CartanMatrixQ::notCartan, txt]];
   	If[Det[A] === 0, PrintWhy["Det(A)=0"]; Throw[False]];
   	Do[If[A[[i, i]] === 2, Null, PrintWhy["A_{ii}!=2 \[ForAll] i"]; Throw[False]], {i, 1, dim}];
   	Do[x = A[[i, j]];
    		If[x > 0, PrintWhy["A_{ij}>0"]; Throw[False]];
    		If[x === 0, If[A[[j, i]] === 0, Null, PrintWhy["\[Exists] A_{ij}=0 with A_{ji}!=0, for some j>i"]; Throw[False]]];
    		If[x < 0, If[MemberQ[{1, 2, 3}, x A[[j, i]]], Null, PrintWhy["\[Exists] A_{ij}<0 with A_{ij}A_{ji}!=1,2,3"]; Throw[False]]]
    	, {i, 1, dim}, {j, i + 1, dim}];
   	True
   ]]


(* ::Text:: *)
(*Strictly speaking, the following duck-typing only checks if the CartanMatrix generates a Finite Lie algebra, Weyl group etc...  Do I want something more general than that?*)


CartanMatrixQLax[A_, why_]:=If[TrueQ[If[why,CartanType[A]==$Failed,Quiet[CartanType[A]==$Failed]]],False,True]


(* ::Subsection::Closed:: *)
(*CartanType*)


CartanType::notCartan = "Not a Cartan matrix because ``";
CartanType::unidentified = "``";


(* ::Text:: *)
(*Option to print out the positions of the components.*)


Options[CartanType]={"Positions"->False};


CartanType[A_,OptionsPattern[]] := Check[Which[
	TrueQ[!MatrixQ[A,MemberQ[{0,2,-1,-2,-3},#]&]], 
		Message[CartanType::notCartan,"it is not an matrix with elements in {0,2,-1,-2,-3}"],
	TrueQ[Det[A]==0],Message[CartanType::notCartan,"it has vanishing determinant (the roots are not lin indep)"],
	TrueQ[Or@@Thread[Diagonal[A]!=2]], Message[CartanType::notCartan,"its diagonal elements are not all 2"],
(* All ok so far.... *)
	True,Module[{cc = ConnectedComponents[A],pos=OptionValue["Positions"]},
   			If[Length@cc == 1,
      				If[TrueQ[pos],{First@cc,#},#]&@SimpleCartanType[A,Range@Length@A],
      				If[TrueQ[pos],#,Last@#]&/@({#,SimpleCartanType[A[[#, #]],#]}& /@ cc)]]],
$Failed,CartanType::notCartan]


(* ::Text:: *)
(*A lot of the final checks (eg in cases for F4 and En) are redundant - as they would have been picked up earlier (eg in Det!=0) - but they don't hurt.*)
(*Also don't really need the Throw[$Failed]'s below.*)


SimpleCartanType[A_,pos_]:=Catch[Module[{r=Length@A, edges, upedges, bedges,deg,ends,PrintDebug=Null},
	If[r==1,If[A=={{2}},Throw["A1"],
		Message[CartanType::notCartan,Row[{"\[Exists] a rank 1 component that's not A1 at position ",pos}]];Throw[$Failed]]];
	If[Length@Cases[Position[A,_?Positive],{i_,j_}/;i!=j]!=0,
		Message[CartanType::notCartan,Row[{"\[Exists] a positive element in component ",pos}]];Throw[$Failed]];
	edges=Position[A//Normal,_?Negative];    
	upedges=Cases[edges,{i_,j_}/;i>j];        PrintDebug["upper-edges: ",upedges,"\nlower-edges: ",Cases[edges,{j_,i_}/;i>j]];
	If[Sort[upedges]!=Sort[Reverse/@Cases[edges,{j_,i_}/;i>j]],
		Message[CartanType::notCartan,Row[{"\[Exists] A_{ij}<0 with A_{ji}=0 in component ",pos}]];Throw[$Failed]];
	If[Length[upedges]>=r, Message[CartanType::notCartan,Row[{"\[Exists] a cylce in component ",pos}]];Throw[$Failed]];
	bedges=Select[edges,A[[Sequence@@#]]<-1&]; PrintDebug["nonsimple upper-edges :",bedges];
	If[Length@bedges>1,
		Message[CartanType::notCartan,Row[{"\[Exists] two elements < -1 in component ",pos}]];Throw[$Failed]];
(* calculate degree of vertices, end-points etc... *)
	deg=Count[Flatten@upedges,#]&/@Range[r]; ends=Flatten@Position[deg,1]; PrintDebug["vertex degrees: ",deg,"\nend posns: ",ends];
	Switch[Length@ends,
		_?(#>3&),Message[CartanType::notCartan,"# of end points in Dynkin diagram \[Element] {2,3}"];Throw[$Failed],
		3,Switch[Length[Intersection[Flatten@Select[upedges,(MemberQ[#,Position[deg,3][[1,1]]]&)],ends]],
			1,Throw["E"<>ToString[r]], (* En,n>8 have det=0, En,n<6 are caught elswhere*)
			2,Throw["D"<>ToString[r]],
			3,If[r==4,Throw["D4"],Message[CartanType::notCartan,"\[Exists] trivalent vertex incident with 3 end points and not D4??"];Throw[$Failed]],
			_,Message[CartanType::notCartan,"\[Exists] trivalent vertex with neither 1, 2 or 3 endpoints"];Throw[$Failed]],
		2,Which[Length@bedges==0,Throw["A"<>ToString[r]],
				A[[Sequence@@bedges[[1]]]]==-3,If[r==2,Throw["G2"],
					Message[CartanType::notCartan,"\[Exists] triple edge but not G2"];Throw[$Failed]],
				A[[Sequence@@bedges[[1]]]]==-2,Switch[Length@Intersection[bedges[[1]],ends],
					0,If[r==4,Throw["F4"],
						Message[CartanType::notCartan,"\[Exists] a central double edges that's not F4, Bn or Cn"];Throw[$Failed]],
					1|2,If[OrderedQ[bedges[[1]]],Throw["B"<>ToString[r]],Throw["C"<>ToString[r]]] (*the 2 case is B2\[TildeEqual]C2*) ]
		]];
	Message[CartanType::unidentified,Short@A];Throw[$Failed]
]]


(* ::Subsection::Closed:: *)
(*CartanMatrixInverse*)


(* ::Text:: *)
(*The inverses are not sparse, so just use Arrays*)


CartanMatrixInverse["A"|"a",r_Integer?Positive]:=Array[If[#1>=#2,(r+1-#1)#2,(r+1-#2)#1]&,{r,r}]/(r+1)
CartanMatrixInverse["B"|"b",r_Integer?(#>1&)]:=Array[Switch[{#2,#1},{r,r-1},r-1,{_,r},#2/2,{_,r-1},#2,{i_,j_}/;i<j,#2,_,#1]&,{r,r}]
CartanMatrixInverse["C"|"c",r_Integer?(#>1&)]:=Array[Switch[{#1,#2},{r,r-1},r-1,{_,r},#1/2,{_,r-1},#1,{i_,j_}/;i<j,#1,_,#2]&,{r,r}]
CartanMatrixInverse["D"|"d",r_Integer?(#>2&)]:=Array[Switch[{#1,#2},{r-1,r-1}|{r,r},r/4,{r,r-1}|{r-1,r},(r-2)/4,
												{_,r|r-1},#1/2,{r|r-1,_},#2/2,{i_,j_}/;i<j,#1,_,#2]&,{r,r}]
CartanMatrixInverse["D"|"d",2]:={{2,0},{0,2}}


CartanMatrixInverse["E"|"e",r_Integer?(#>2&)]:=Inverse[CartanMatrix["E",r]]
CartanMatrixInverse["F"|"f",4]:=Inverse[CartanMatrix["F",4]]
CartanMatrixInverse["G"|"g",2]:=Inverse[CartanMatrix["G",2]]


CartanMatrixInverse[{type_String,rank_Integer?Positive},opts:OptionsPattern[]]:=CartanMatrixInverse[type,rank,opts]


CartanMatrixInverse[blocks:{{_String,_Integer?Positive}..},opts:OptionsPattern[]]:=Module[{x},
	ArrayFlatten[DiagonalMatrix[Array[x,Length@blocks]]/.x[i_]:>CartanMatrixInverse[Sequence@@blocks[[i]],opts]]]


CartanMatrixInverse[type_String,opts:OptionsPattern[]]/;StringLength[type]>1:=CartanMatrixInverse[CartanMatrixStringSplit[type],opts]
CartanMatrixInverse[type:{_String..},opts:OptionsPattern[]]:=CartanMatrixInverse[StringJoin[type],opts]


(* ::Text:: *)
(*Following code doesn't work because ConnectedComponents does not return each component in standard order...*)


(* ::Text:: *)
(*CartanMatrixInverse[A_?(ArrayQ[#, 2, IntegerQ] &)] := Module[{order, type},*)
(*    		{order, type} = CartanType[A, "Positions" -> True]\[Transpose];*)
(*    		order = SortBy[{Flatten[order], Range[Length@A]}\[Transpose], First][[All, 2]];*)
(*    		CartanMatrixInverse[type][[order, order]]]*)


(* ::Subsubtitle::Closed:: *)
(*Graphs and Plots -- Not finished*)


(* ::Subsection::Closed:: *)
(*ConnectedComponents*)


(* ::Text:: *)
(*Given a graph, we make it simple - ie remove self-loops and make all edges point from the lower numbered vertex to the higher *)
(*	(not really nesc. just include in the fold? - about 10% faster).*)
(*Fully contract all edges - the vertices that are left are the disconnected components of the original graph.*)


ConnectedComponents[graph:{(_->_)..}]:=Module[{edges,verts},
(* keep simple graph - take long way round for edges to handle lone self-loops *)
	edges=Union[Sort/@(graph/.Rule->List)];
	verts=Union@Flatten[edges]; 
	edges=DeleteCases[edges,{i_,i_}];
(* use every edge to identify vertices - equiv to fully contracting graph *)
	verts=Fold[#1/.Rule@@#1[[#2]]&,verts,edges];
(* number of connected components in contracted graph = that in original graph *)
	Flatten[Position[verts,#]]& /@ Union[verts]
]


ConnectedComponents[M_?MatrixQ]:=ConnectedComponents[Rule@@#&/@Position[M,_?(#!=0&)]]


(* ::Subsection::Closed:: *)
(*DynkinDiagram*)


(* ::Text:: *)
(*The default numbering follows Dynkin's conventions (eg http://www-math.mit.edu/~lesha/dynkin-diagrams.html).*)
(*In this program, that follows from the choice of Cartan matrices.*)


Options[DynkinDiagram]={"Labels"->True,"RowOrColumn"->Column};


DynkinDiagram::rootlengths = "more than two root lengths in component ``";
DynkinDiagram::neglength = "there is edge with negative length in component ``";


(* ::Text:: *)
(*Define VertexCoordinateRules for our built in Cartan matrices:  a,b,c,f and g are straight lines;  d and e are treated separately.*)


VCRWhich[A_,rnk_]:=Quiet[Which[
	MemberQ[CartanMatrix[#,rnk]&/@{"a","b","c","f","g"},A],Table[i->{i,0},{i,rnk}],
	A==CartanMatrix["d",rnk],Join[Table[i->{i,0},{i,rnk-2}],{rnk-1->{rnk-Sqrt[1-#^2],-#},rnk->{rnk-Sqrt[1-#^2],#}}&[.3]],
	A==CartanMatrix["e",rnk],Append[Table[i->{i,0},{i,rnk-1}],rnk->{3,.6}],
	_,Automatic
],CartanMatrix::"rank"]


(* ::Text:: *)
(*Define the VertexRenderingFunction and EdgeRenderingFunction:*)
(*The edges are of a fixed thickness and colour and are either dotted (with solid beginning and end) or not.  Dotting only works (and makes sense) for single lines.*)
(*The vertices are filled*)


DynkinEdgeRenderingFunction[dots:{{_Integer,_Integer}..},thick_:Thickness[.004],colour_:Black]:=Module[{lne,t,p=.25},
	lne[t_,coord_]:=t First[Sort@coord]+(1-t)Last[Sort@coord];
	({If[MemberQ[dots,Sort[#2]],
		{thick,Line[{lne[0,#],lne[p,#]}],Line[{lne[1,#],lne[1-p,#]}],Dotted,Line[{lne[p,#],lne[1-p,#]}]},
		{Dashing[{}],thick,Line[#]}]}&)]


(* ::Text:: *)
(*long is a greyscale: 1 = white, 0 = black*)


DynkinVertexRenderingFunction[long_,labels_,thick_:Thickness[.004],labelsQ_:True]:=({
	GrayLevel[long[[#2]]], 
	EdgeForm[Directive[thick,Black]],Disk[#,.1],
	Sequence@@If[labelsQ,{Black,Style[Text[labels[[#2]],#1+.15{1,1.5}],
	FontFamily->Times]},{}]}&)


(* ::Text:: *)
(*Diagrams*)


With[{NumMatQ=MatrixQ[#]\[And]NumericQ[Total[Flatten[#]]]&},
DynkinDiagram[A_?NumMatQ,opts:OptionsPattern[]]:=DynkinDiagram[A,Range[Length@A],opts];
DynkinDiagram[A_?NumMatQ,vLabels_List,opts:OptionsPattern[]]:=DynkinDiagram[A,vLabels,{{0,0}},opts];
DynkinDiagram[A_?NumMatQ,vLabels_List,dots:{_Integer,_Integer},opts:OptionsPattern[]]:=DynkinDiagram[A,vLabels,{dots},opts];
DynkinDiagram[A_?NumMatQ,vLabels_List,dots:{{_Integer,_Integer}..},opts:OptionsPattern[]]/;Length[A]==Length[vLabels]:=Module[
	{B,comps,num,AA,len,i,j,lengths,long,RowColumn=OptionValue[DynkinDiagram,"RowOrColumn"]},
	If[MemberQ[{Row,Column},RowColumn],RowColumn,Column];
	comps=ConnectedComponents[A]; (*Print[comps];*)
(* Calculate the lengths for each connected component AND generate edges *)
	Do[AA=(A[[cc,cc]]);len=Length@AA;B[cc]={};
		lengths=StandardBasis[{1},len,"Sparse"->False]; (* first length is init set to 1, the rest to zero *)
		For[i=1,i<=len,i++,
			For[j=i+1,j<=len,j++,
				If[A[[i,j]]A[[j,i]]<0,Message[DynkinDiagram::neglength,cc];Return[$Failed]];
				B[cc]=Nest[Append[#,i->j]&,B[cc],AA[[i,j]]AA[[j,i]]];
				If[AA[[i,j]]!=0,lengths[[j]]=lengths[[i]]AA[[j,i]]/AA[[i,j]]]]];
		long[cc]=If[Length@Union@lengths>2,
			Message[DynkinDiagram::rootlengths,cc];Flatten@((Position[Union@lengths,#]&/@lengths-1)/(Length@Union@lengths-1)),
			long[cc]=(Boole[#==Max[lengths]])&/@lengths];
			(*If[B[cc]=={},B[cc]={{0}}];*) (*Print[cc," : ",B[cc]," : ",long[cc]];*)
		,{cc,comps}];
(* if there is more than one component, then output a GraphicsRow/Column *)
	If[Length@comps>1,RowColumn[#(*,Spacings->{0,0},ImageMargins->None*)]&,First][
	MapIndexed[If[TrueQ[B[#]=={}],
		Graphics[{Text[vLabels[[#[[1]]]],{0,0},{-1,-.5}],Text[Style["\[EmptyCircle]",Medium],{0,-1}](*Circle[{0,0},.5]*)},ImageSize->30],
		GraphPlot[B[#],
			MultiedgeStyle->True,ImageSize->50len,VertexCoordinateRules->VCRWhich[A[[#,#]],len],
			EdgeRenderingFunction->DynkinEdgeRenderingFunction[Sort/@dots,Thickness[.003]],
			VertexRenderingFunction->DynkinVertexRenderingFunction[long[#],
				With[{x=Length[Flatten[comps[[;;#2[[1]]-1]]]]},vLabels[[1+x;;x+Length[#1]]]],
				Thickness[.003],OptionValue["Labels"]],
			ImageMargins->None]]&,comps]]
]]


(* ::Text:: *)
(*Maybe make the following more sensible - ie don't actually generate Cartan matrix etc... -- probably not worth the time!*)


DynkinDiagram[type_?CartanBlocksQ,etc___]:=DynkinDiagram[CartanMatrix[type],etc]


(* ::Subsection::Closed:: *)
(*GraphRoots2D and GraphRoots3D*)


(* ::Text:: *)
(*GraphRoots takes in roots in SRS rep.*)


(*angle from {1,0}*)
ang[x_]:=Sign[x[[2]]] ArcCos[x.{1,0}/Sqrt[x.x]]+\[Pi] KroneckerDelta[0,x[[2]]] KroneckerDelta[Sign[x[[1]]],-1]


GraphRoots2D::rank = "GraphRoots2D is only designed to work with rank two algebras";
GraphRoots3D::rank = "GraphRoots2D is only designed to work with rank three algebras";


GraphRoots2D[type_?CartanBlocksQ,symb_:"\[Alpha]",style_:{},rt1_:{1,0}] := GraphRoots2D[CartanMatrix[type],symb,style,rt1]
GraphRoots2D[A_?(ArrayQ[#,2,IntegerQ]&),symb_:"\[Alpha]",style_:{},rt1_:{1,0}] := Module[{Q,R,rts,sRts,r=Length@A},
	If[r!=2,Message[GraphRoots2D::rank];Return[$Failed]];
	rts=Join[#,-#]&@Flatten[PositiveRoots[A],1];
	Q=MatrixPower[A.DiagonalMatrix[SimpleRootLengths[A]]/2//N,1/2];(*Print@Q;*)
	R=RotationMatrix[{Q[[1]],rt1}];
	sRts=R.#&/@Q//Chop;rts=Sort[R.Q.#&/@rts,ang[#1]>ang[#2]&]//Chop;(*Print@rts;*)
	Graphics[{MapIndexed[Text[Style[Subscript[symb,#2[[1]]],style],#,-1.5Normalize[#]]&,sRts],
		Gray,Line[{{0,0},#}]&/@rts,
		Red,(Line[{{0,0},#}]&/@sRts),
		Green,Line[rts~Join~{rts[[1]]}]}]
]


GraphRoots3D[type_?CartanBlocksQ,symb_:"\[Alpha]",style_:{},rt1_:{1,0,0}] := GraphRoots3D[CartanMatrix[type],symb,style,rt1]
GraphRoots3D[A_?(ArrayQ[#,2,IntegerQ]&),symb_:"\[Alpha]",style_:{},rt1_:{1,0,0}] := Module[{Q,R,sRts,rts,r=Length@A},
	If[r!=3,Message[GraphRoots3D::rank];Return[$Failed]];
	rts=Join[#,-#]&@Flatten[PositiveRoots[A],1];
	Q=MatrixPower[A.DiagonalMatrix[SimpleRootLengths[A]]/2//N,1/2];(*Print@Q;*)
	R=RotationMatrix[{Q[[1]],rt1}];
	sRts=R.#&/@Q//Chop;rts=R.Q.#&/@rts//Chop;(*Print@rts;*)
	Graphics3D[{MapIndexed[Text[Style[Subscript[symb,#2[[1]]],style],1.1#]&,sRts],
	Gray,Line[{{0,0,0},#}]&/@rts[[4;;]],Gray,Red,Line[{{0,0,0},#}]&/@rts[[1;;3]]},Boxed->False]
]


(* ::Text:: *)
(*TODO*)


GraphRootsCoxeter[A_?(ArrayQ[#,2,IntegerQ]&)]:=Module[{},Null]


(* ::Subsection::Closed:: *)
(*GraphWeights2D and GraphWeights3D -- TODO*)


GraphWeights2D[A_?(ArrayQ[#,2,IntegerQ]&),weights_List,pw1_:{1,0}] := Module[{wts=Flatten[weights,1],r,
	Ai=Inverse[A](*,pw=PrimitiveWeights[A]*)},
		r=Length@wts;
		Graphics[]]


(* ::Subsection:: *)
(*Projection onto 2D*)


(* ::Subsection:: *)
(*Projection onto 3D*)


(* ::Subsubtitle::Closed:: *)
(*Roots*)


(* ::Text:: *)
(*SimpleRootLengths and SimpleRootInnerProduct can be effectively accessed through RootInnerProduct - the public function.*)
(*All of the above depend on the choice of root normalisation.*)


(* ::Text:: *)
(*PositiveRoots, on the other hand, does not depend on root normalisation and can be carried out using any repn*)


ClearRootCaches[] := (reloadSRLP[]; reloadRIP[];reloadPRALG[];reloadPRSRS[])


(* ::Subsection::Closed:: *)
(*SimpleRootLengths and SimpleRootInnerProduct*)


(* ::Text:: *)
(*Each simple component is normalised so that the shortest root has length squared 1 --- the rest then have length squared 2 or 3*)
(*The other common option is to normalise the long root to length squared 2 -- this simplifies some other formulas.*)
(*Of course neither of these lengths agree with a root innerproduct defined as the dual of the Killiing form.*)
(*Another option is to leave it unfixed by either normalizing the shortest or longest root to a variable OR by using the option "norm"->None.*)
(*In the latter case each component has an associated RootLengthSquared[root]*)


SimpleRootLengths::notCartan = "`` is not a valid Cartan matrix";
SimpleRootLengths::opt = "Bad option: ``";


Options[SimpleRootLengths]={"norm"->"Short","normTo"->1}; (*"norm"->"Long","normTo"->2*) (*"norm"->None *)
Options[SimpleRootLengthsPrivate]=Options[SimpleRootLengths];


(* ::Text:: *)
(*RootLengthSquared is a placeholder for unknown lengths when using the option "norm"->None|"None"*)


RootLengthSquared/: MakeBoxes[RootLengthSquared[root_Integer],form_] := InterpretationBox[SuperscriptBox[RowBox[{"\[LeftBracketingBar]",SubscriptBox["\[Alpha]", root],"\[RightBracketingBar]"}],2], RootLengthSquared[root]]
RootLengthSquared/: MakeBoxes[RootLengthSquared[root_],form_] := With[{rt=ToBoxes@root},InterpretationBox[SuperscriptBox[RowBox[{"\[LeftBracketingBar]",rt,"\[RightBracketingBar]"}],2], RootLengthSquared[root]]]


(* ::Text:: *)
(*Just returns a list of squared lengths corresponding to the simple roots in the order given by the Cartan matrix.*)
(*Use dynamic programming to speed things up...  minimal memory usage...  Has a slightly slower work around since SimpleRootLengths is Protected.*)


reloadSRLP[]:=(Clear[SimpleRootLengthsPrivate];SimpleRootLengthsPrivate[A_?MatrixQ,opts:OptionsPattern[]]:=
SimpleRootLengthsPrivate[A,opts]=Catch[
Module[{comps,cc,norm,normsQ,AA,i,j,len,lengths,result={}},
	comps=ConnectedComponents[A];  (*Print[comps]; *)
	normsQ=TrueQ[Length[OptionValue["normTo"]]>=Length[comps]];
(* Calculate the lengths for each connected component and generate edges *)
	Do[cc=comps[[i]];
		If[normsQ,norm=OptionValue["normTo"][[i]],norm=First[Flatten[{OptionValue["normTo"]}]]];
		AA=(A[[cc,cc]]);len=Length@AA;
		lengths=StandardBasis[{1},len,"Sparse"->False]; (* first length is init set to 1, the rest to zero *)
		For[i=1,i<=len,i++,For[j=i+1,j<=len,j++,
			Which[AA[[i,j]]<0&&AA[[j,i]]<0,lengths[[j]]=lengths[[i]] AA[[j,i]]/AA[[i,j]],
				AA[[i,j]]<0&&AA[[j,i]]>=0||AA[[i,j]]>=0&&AA[[j,i]]<0,
					Message[SimpleRootLengths::notCartan,Short[AA]];Throw[$Failed],
				_, Null]]];
		If[Length@Union@lengths[cc]>2,Print[{"!! More than two root lengths in component ",cc," !!"}]];
		AppendTo[result,Switch[OptionValue["norm"],
			None|"None",RootLengthSquared[First@cc]lengths,
			"Short",norm lengths/Min[DeleteCases[lengths,0]],
			"Long",norm lengths/Max[lengths,0],
			_,Message[SimpleRootLengths::opt,"norm"->OptionValue["norm"]];Throw[$Failed]]]
		,{i,1,Length@comps}]; 
	Flatten@result]])
reloadSRLP[]


SimpleRootLengths[A_?MatrixQ,opts:OptionsPattern[]]:=SimpleRootLengthsPrivate[A,opts]


SimpleRootLengths[type_?CartanBlocksQ,etc___]:=SimpleRootLengths[CartanMatrix[type],etc]


SimpleRootInnerProduct[A_?MatrixQ,opts:OptionsPattern[]][i_,j_]:=SimpleRootLengths[A,opts][[j]]A[[i,j]]/2


SimpleRootInnerProduct[type_?CartanBlocksQ,etc___][i_,j_]:=SimpleRootInnerProduct[CartanMatrix[type],etc][i,j]


(* ::Subsection::Closed:: *)
(*RootInnerProduct*)


(* ::Text:: *)
(*Here we need to supply a List of labels for the simple roots.*)


RootInnerProduct[A_?MatrixQ,simpRoots_List,opts:OptionsPattern[]][\[Alpha]_,\[Beta]_]/;And[Length@A==Length@simpRoots,And@@(MemberQ[simpRoots,#]&/@{\[Alpha],\[Beta]})]:=RootInnerProductPrivate[A,simpRoots,opts][\[Alpha],\[Beta]]


reloadRIP[]:=(Clear[RootInnerProductPrivate];
	RootInnerProductPrivate[A_,simpRoots_List,opts:OptionsPattern[]][\[Alpha]_,\[Beta]_]:=RootInnerProductPrivate[A,simpRoots,opts][\[Alpha],\[Beta]]=SimpleRootInnerProduct[A,opts][Position[simpRoots,\[Alpha]][[1,1]],Position[simpRoots,\[Beta]][[1,1]]]/.RootLengthSquared[n_Integer]:>RootLengthSquared[simpRoots[[n]]])
reloadRIP[]


RootInnerProduct[A_,simpRoots_,opts:OptionsPattern[]][\[Alpha]_+\[Beta]_,\[Gamma]_]:=RootInnerProduct[A,simpRoots,opts][\[Alpha],\[Gamma]]+RootInnerProduct[A,simpRoots,opts][\[Beta],\[Gamma]]
RootInnerProduct[A_,simpRoots_,opts:OptionsPattern[]][\[Alpha]_,\[Gamma]_+\[Beta]_]:=RootInnerProduct[A,simpRoots,opts][\[Alpha],\[Beta]]+RootInnerProduct[A,simpRoots,opts][\[Alpha],\[Gamma]]
RootInnerProduct[A_,simpRoots_,opts:OptionsPattern[]][a_ \[Alpha]_,\[Beta]_]:= a RootInnerProduct[A,simpRoots,opts][\[Alpha],\[Beta]]
RootInnerProduct[A_,simpRoots_,opts:OptionsPattern[]][\[Alpha]_,a_ \[Beta]_]:= a RootInnerProduct[A,simpRoots,opts][\[Alpha],\[Beta]]


RootInnerProduct[type_?CartanBlocksQ,etc___][i_,j_]:=RootInnerProduct[CartanMatrix[type],etc][i,j]


(* ::Subsection::Closed:: *)
(*PositiveRoots*)


(* ::Text:: *)
(*Calculate: 	algebraically 			"ALG", *)
(*		basis of simple roots 		"SRS" -- about 10 times faster*)
(*		Dynkin 				"DYN" -- same as SRS, just uses inverse Cartan matrix. No point really.*)


(* ::Text:: *)
(*PostiveRootsALG just uses the default normalization option in SimpleRootLengths, since it's normalization indep.*)


Options[PositiveRoots]={Method->"SRS","maxLevel"->150,"Debug"->False};
Options[PositiveRootsALG]={"maxLevel"->150,"Debug"->False};
Options[PositiveRootsSRS]={"maxLevel"->150,"Debug"->False};


PositiveRoots::maxLevel = "Reached the maximum root level of ``. Can increae using the option \[OpenCurlyDoubleQuote]maxLevel\[CloseCurlyDoubleQuote]->lev, lev=1,\[Ellipsis],\[Infinity]";
PositiveRoots::simpALG = "When using the algebraic method, you need to supply symbols for the simple roots.";
PositiveRoots::simp = "The number of supplied symbols for the simple roots does not match the rank of A.";
PositiveRoots::DYN = "DYN method is not implemented.";
PositiveRoots::Method = "`` is an invalid Method.";


PositiveRoots[A_?(ArrayQ[#,2,IntegerQ]&),simp_List:{},opts:OptionsPattern[]]:=Switch[OptionValue[Method],
	"ALG",If[Length[A]==Length[simp],
			PositiveRootsALG[A,simp,FilterRules[{opts},Except[Method]]],
			Message[PositiveRoots::simpALG];Return[$Failed]],
	"SRS"|Default,If[simp=={}||Length@simp==Length@A,
			PositiveRootsSRS[A,simp,FilterRules[{opts},Except[Method]]],
			Message[PositiveRoots::simpALG];Return[$Failed]],
	"DYN",Message[PositiveRoots::DYN];Return[$Failed],
	_,Message[PositiveRoots::Method,OptionValue[Method]];Return[$Failed]
]


PositiveRoots[type_?CartanBlocksQ,etc___]:=PositiveRoots[CartanMatrix[type],etc]


reloadPRALG[]:=(Clear[PositiveRootsALG];
PositiveRootsALG[A_?(ArrayQ[#,2,IntegerQ]&),simp_List,opts:OptionsPattern[]]:=PositiveRootsALG[A,simp,opts]=Module[
	{len=Length[A],maxLev=OptionValue["maxLevel"],pRts,lev=1,p,\[Alpha],\[Beta],PrintDebug},
	If[TrueQ[OptionValue["Debug"]],PrintDebug=Print,PrintDebug=Null&];
	pRts[0]={};pRts[1]=simp;
	While[lev<=maxLev,lev++;pRts[lev]={};                     PrintDebug["level: ",lev-1,", roots: ",pRts[lev-1]];
		Do[Do[p=1;While[MemberQ[pRts[lev-(p+1)],\[Beta]-p \[Alpha]],p++];  PrintDebug[{\[Beta],\[Alpha],p}];
			If[(p-1)-2 RootInnerProduct[A,simp][\[Alpha],\[Beta]]/RootInnerProduct[A,simp][\[Alpha],\[Alpha]]>0,
				pRts[lev]=Union[pRts[lev],{\[Beta]+\[Alpha]}]],
		{\[Alpha],simp}],{\[Beta],pRts[lev-1]}];
	If[pRts[lev]=={},lev--;Break[]]];
	If[lev>maxLev,lev--;Message[PositiveRoots::maxLevel,OptionValue["maxLevel"]]];
	Table[pRts[n],{n,lev}]])
reloadPRALG[]


reloadPRSRS[]:=(Clear[PositiveRootsSRS];
PositiveRootsSRS[A_?(ArrayQ[#,2,IntegerQ]&),simp_List,opts:OptionsPattern[]]:=Table[Sort[simp.#&/@lev],{lev,PositiveRootsSRS[A,{},opts]}];
PositiveRootsSRS[A_?(ArrayQ[#,2,IntegerQ]&),{},opts:OptionsPattern[]]:=PositiveRootsSRS[A,{},opts]=Module[
	{len=Length[A],maxLev=OptionValue["maxLevel"],pRts,lev=1,p,\[Alpha],\[Beta],\[Beta]A,PrintDebug},
	If[TrueQ[OptionValue["Debug"]],PrintDebug=Print,PrintDebug=Null&];
	pRts[0]={};pRts[1]=IdentityMatrix[len];
	While[lev<=maxLev,lev++;pRts[lev]={};                      PrintDebug["level: ",lev-1,", roots: ",pRts[lev-1]];
		Do[\[Beta]A=\[Beta].A;
			Do[p=1;While[MemberQ[pRts[lev-(p+1)],\[Beta]-p*\[Alpha]],p++];  PrintDebug[{\[Beta],\[Alpha],p}];
				If[(p-1)-\[Beta]A.\[Alpha]>0,pRts[lev]=Union[pRts[lev],{\[Beta]+\[Alpha]}]],
		{\[Alpha],pRts[1]}],{\[Beta],pRts[lev-1]}];
	If[pRts[lev]=={},lev--;Break[]];];
	If[lev>maxLev,lev--;Message[PositiveRoots::maxLevel,OptionValue["maxLevel"]]];
	Table[pRts[n],{n,lev}]])
reloadPRSRS[]


(* ::Subsubtitle::Closed:: *)
(*Weyl group*)


(* ::Text:: *)
(*Weyl reflections are all implemented as matrices acting on the right of Weights in the Dynkin representation.*)
(*The second argument is the simple root to be reflected about.*)


WeylReflection[A_?(ArrayQ[#,2,IntegerQ]&),i_Integer]/;0<i<=Length@A := Module[{R=IdentityMatrix[Length[A]]},
		R[[i]] = R[[i]] - A\[Transpose][[i]]; R]
WeylReflection[type_String,ord_Integer,i_Integer] := WeylReflection[CartanMatrix[{type,ord}],i]
WeylReflection[type:{{_String,_Integer}..},i_Integer] := WeylReflection[CartanMatrix[type],i]
WeylReflection[type_String,i_Integer] := WeylReflection[CartanMatrix[type],i]


WeylReflection[A_?(ArrayQ[#,2,IntegerQ]&),i_?(ArrayQ[#,1,IntegerQ]&)]:=WeylReflection[A,#]&/@i
WeylReflection[type_?CartanBlocksQ,i_?(ArrayQ[#,1,IntegerQ]&)]:=WeylReflection[type,#]&/@i


SimpleReflections[A_?(ArrayQ[#,2,IntegerQ]&)]:=WeylReflection[A,Range[Length[A]]]
SimpleReflections[type_?CartanBlocksQ]:=SimpleReflections[CartanMatrix[type]]


WeylOrbit[A_?(ArrayQ[#,2,IntegerQ]&),w_?(ArrayQ[#,1,IntegerQ]&),max_:100]:=With[
	{sr=SimpleReflections[A]},
		Reverse@FixedPoint[Union[Flatten[Append[Table[#.r,{r,sr}],#]&/@#,1]]&,{w},max]]
WeylOrbit[type_?CartanBlocksQ,w_?(ArrayQ[#,1,IntegerQ]&),max_:100]:=WeylOrbit[CartanMatrix[type],w,max]


(* ::Subsubtitle:: *)
(*Weights and reps*)


(* ::Subsection::Closed:: *)
(*FundamentalWeights (aka PrimitiveWeights)*)


(* ::Text:: *)
(*Returns the fundamental weights as components of their expansion in simple roots (ie the SRS rep)*)
(*If given a list of symbols for simple roots, then returns the simple root expansion*)


PrimitiveWeights = FundamentalWeights;


FundamentalWeights::length = "The rank of the algebra does not match the supplied number of symbols for the simple roots";


FundamentalWeights[A_?(ArrayQ[#,2,IntegerQ]&)] := Inverse[A]
FundamentalWeights[type_?CartanBlocksQ]:=CartanMatrixInverse[type]


FundamentalWeights[A_?(ArrayQ[#,2,IntegerQ]&),sr_List]:=(
		If[TrueQ[Length[sr]==Length[A]],
			FundamentalWeights[A].sr,
			Message[FundamentalWeights::length];Return[$Failed]])
FundamentalWeights[type_?CartanBlocksQ,sr_List]:=With[{fw=FundamentalWeights[type]},
		If[TrueQ[Length[sr]==Length[fw]],
			fw.sr,
			Message[FundamentalWeights::length];Return[$Failed]]]


(* ::Subsection::Closed:: *)
(*WeylVector*)


(* ::Text:: *)
(*Returns the Weyl vector, \[Rho], in the SRS rep.*)
(*In the DYN rep it has components (1,1,1,1...,1)*)


WeylVector[A_?(ArrayQ[#,2,IntegerQ]&)] := Total[Flatten[PositiveRoots[A], 1]]/2
WeylVector[type_?CartanBlocksQ] := Total[Flatten[PositiveRoots[type], 1]]/2


WeylVector[A_?(ArrayQ[#,2,IntegerQ]&),sr_?VectorQ] := sr.WeylVector[A]
WeylVector[type_?CartanBlocksQ,sr_?VectorQ] := sr.WeylVector[type]


WeylVector[A_?(ArrayQ[#,2,IntegerQ]&),"DYN"] := Array[1&,Length@A]
WeylVector[A_?CartanBlocksQ,"DYN"] := Array[1&,Length@CartanMatrix[A]]


(* ::Subsection::Closed:: *)
(*Weight Orderings ???*)


(* ::Text:: *)
(*A Lexographic ordering*)


LexPositive[w_VectorQ]:=Quiet[Check[DeleteCases[w,0][[1]]>0,False,{Part::"partw"}]]


LexOrderedQ[w:{_?VectorQ..}]:=Catch[Quiet[Check[
    		(Greater@@#[[All,1]]&)@NestWhile[#[[All, 2;;]]&, {w}, Equal@@#[[All,1]]&],
    		Throw[True], {Part::"partw"}]]]
LexOrderedQ[w:(_?VectorQ..)]:=LexOrderedQ[{w}]


LexSort[w:{_?VectorQ..}]:=Sort[w,SRSOrderedQ]


(* ::Text:: *)
(*Lexographic ordering in a Chevalley basis*)


SRSPositive[A_?(ArrayQ[#,2,IntegerQ]&),w_VectorQ]:=LexPositive[w.A]
SRSPositive[type_?CartanBlocksQ,w_VectorQ]:=LexPositive[w.CartanMatrix[type]]
SRSOrderedQ[A_?(ArrayQ[#,2,IntegerQ]&),w:{_?VectorQ..}]:=LexOrderedQ[#.A&/@w]
SRSOrderedQ[type_?CartanBlocksQ,w:{_?VectorQ..}]:=SRSOrderedQ[CartanMatrix[type],w]
SRSSort[A_?(ArrayQ[#,2,IntegerQ]&),w:{_?VectorQ..}]:=Sort[w,SRSOrderedQ[A,#]&]
SRSSort[type_?CartanBlocksQ,w:{_?VectorQ..}]:=Sort[w,SRSOrderedQ[CartanMatrix[type],#]&]


(* ::Text:: *)
(*Checks ordering of any number of weights in DYN rep (primitive weight expansion).*)
(*Actually, of course, works on any set of vectors.*)


(* ::Subsection:: *)
(*DominantWeights*)


(* ::Subsection::Closed:: *)
(*SpindleForm*)


SpindleForm[expr_List]:=Column[Row[Identity[#],"  ,  "]&/@expr,Alignment->Center]
SpindleForm[expr_]:=expr


(* ::Subsection::Closed:: *)
(*WeightSystem & *)


(* ::Text:: *)
(*Could speed up the "SRS" type results by using CartanMatrixInverse for simple algebras...  but would be minimal.*)


WeightSystem::nonDom = "You have supplied a non-dominant weight.";
WeightSystem::dims = "The dimensions of the Cartan matrix, the highest weight and (optional) the simple root labels must match.";
WeightSystem::maxLayer = "Reached the maximum weight layer of ``. Can increae using the option \[OpenCurlyDoubleQuote]maxLayer\[CloseCurlyDoubleQuote]->layer, layer=1,\[Ellipsis],\[Infinity]";
WeightSystem::Method = "`` is an invalid or not implemented method.";


Options[WeightSystem]={"maxLayer"->150,"Debug"->False,Multiplicity->False,Method->"DYN"};
Options[WeightSystemNoMultDYN]={"maxLayer"->150,"Debug"->False};
Options[WeightSystemDYN]={"maxLayer"->150,"Debug"->False};


WeightSystem[A_?(ArrayQ[#,2,IntegerQ]&),M_?(ArrayQ[#,1,IntegerQ]&),opts:OptionsPattern[]]:=Switch[OptionValue[Method],
	"DYN"|Default,If[TrueQ[Length[A]==Length[M]],
			If[OptionValue[Multiplicity]===False,
				WeightSystemNoMultDYN[A,M,FilterRules[{opts},Except[Method|Multiplicity]]],
				WeightSystemDYN[A,M,FilterRules[{opts},Except[Method|Multiplicity]]]],
			Message[WeightSystem::dims];Return[$Failed]],
	"SRS",WeightSystem[A,M,"SRS",FilterRules[{opts},Except[Method]]],
	_,Message[WeightSystem::Method,OptionValue[Method]];Return[$Failed]
]


WeightSystem[A_?(ArrayQ[#,2,IntegerQ]&),M_?(ArrayQ[#,1,IntegerQ]&),simps_List,opts:OptionsPattern[]]:=If[
	Length@A==Length@M==Length@simps,
	With[{Ai=Inverse[A],weights=WeightSystem[A,M,opts]},Map[#.Ai.simps&,weights,{2}]],
	Message[WeightSystem::dims]]


WeightSystem[A_?(ArrayQ[#,2,IntegerQ]&),M_?(ArrayQ[#,1,IntegerQ]&),"SRS",opts:OptionsPattern[]]:=If[
	Length@A==Length@M,
	With[{Ai=Inverse[A],weights=WeightSystem[A,M,opts]},Map[#.Ai&,weights,{2}]],
	Message[WeightSystem::dims]]


WeightSystem[type_?CartanBlocksQ,etc__]:=WeightSystem[CartanMatrix[type],etc]


(* ::Text:: *)
(*Naive construction of Weights WITHOUT multiplicities.  Takes Cartan matrix (or named blocks) + Highest weight.*)
(*Works in Dynkin rep where (\[Alpha]_i)_j=A_ {i j}*)


WeightSystemNoMultDYN[A_,M_,OptionsPattern[]]:= (
	If[Min[M]<0,Message[WeightSystem::nonDom]];
	Module[{r=Length@A,weights,layer=0,q,maxLayer=OptionValue["maxLayer"],PrintDebug},
		If[TrueQ[OptionValue["Debug"]],PrintDebug=Print,PrintDebug=Null&];
		weights[0]={M}; 
		While[layer<=maxLayer,layer++;weights[layer]={};           
		PrintDebug["layer: ",layer-1,", weights: ",weights[layer-1]];
			Do[Do[q=0;While[MemberQ[weights[layer-2-q],\[Lambda]+(q+1)*A[[j]]],q++];  
				PrintDebug["\[Lambda]=",\[Lambda]," \[Alpha]=",A[[j]]," q=",q," p=",\[Lambda][[j]]+q];
				If[\[Lambda][[j]]+q>0,weights[layer]=Union[weights[layer],{\[Lambda]-A[[j]]}]],
				{j,r}(*{\[Alpha],A}*)],{\[Lambda],weights[layer-1]}];
		If[weights[layer]=={},layer--;Break[]]];
		If[layer>maxLayer,layer--;Message[WeightSystem::maxLayer,OptionValue["maxLayer"]]];
		Table[weights[n],{n,0,layer}]
])


(* ::Text:: *)
(*http://en.wikipedia.org/wiki/Weyl_character_formula#Freudenthal.27s_formula*)


WeightSystemDYN[A_,M_,OptionsPattern[]]:= (
	If[Min[M]<0,Message[WeightSystem::nonDom]];
	Module[{r=Length@A,weights,mult,layer=0,q,maxLayer=OptionValue["maxLayer"],PrintDebug,
			\[Rho]=WeylVector[A,"DYN"],pRts=Map[#.A&,Flatten[PositiveRoots[A],1]],Ai=Inverse[A],G},
		If[TrueQ[OptionValue["Debug"]],PrintDebug=Print,PrintDebug=Null&];
		G=Ai.DiagonalMatrix[SimpleRootLengths[A,"norm"->None]]/2; (* metric for DYN rep *)
		weights[0]={M}; mult[M]=1;
		While[layer<=maxLayer,layer++;weights[layer]={};   
			PrintDebug["layer: ",layer-1,", weights: ",weights[layer-1]," mults: ",mult/@weights[layer-1]];        
			Do[Do[q=0;While[MemberQ[weights[layer-2-q],\[Lambda]+(q+1)*A[[j]]],q++];  
				If[\[Lambda][[j]]+q>0,weights[layer]=Union[weights[layer],{\[Lambda]-A[[j]]}]],
				{j,r}(*{\[Alpha],A}*)],{\[Lambda],weights[layer-1]}];
			If[weights[layer]=={},layer--;Break[]];
			With[{wts=Flatten[Table[weights[n],{n,0,layer-1}],1]},
				Do[mult[\[Lambda]]=2((M+\[Lambda]+2\[Rho]).G.(M-\[Lambda]))^(-1) Sum[Sum[mult[x](x.G.\[Alpha]),
				{x,Intersection[wts,Table[\[Lambda]+k*\[Alpha],{k,layer}]]}],{\[Alpha],pRts}],{\[Lambda],weights[layer]}]]
		];
		If[layer>maxLayer,layer--;Message[WeightSystem::maxLayer,OptionValue["maxLayer"]]];
		Map[Sequence@@Table@@{#,{mult[#]}}&,Table[weights[n],{n,0,layer}],{2}]
])


(* ::Text:: *)
(*(M + \[Rho]).G.(M + \[Rho]) - (\[Lambda] + \[Rho]).G.(\[Lambda] + \[Rho]) = (M + \[Lambda] + 2 \[Rho]).G.(M - \[Lambda])*)


(* ::Subsection::Closed:: *)
(*WeightMultiplicity*)


(* ::Text:: *)
(*Not designed to be optimal, but rather designed to be informative when using debug=True*)


WeightMults[A_?(ArrayQ[#,2,IntegerQ]&),wts_,debug_:False]:=
Module[{mult,M=wts[[1,1]],G,up,num,T=Length@wts,pRts=Map[#.A&,Flatten[PositiveRoots[A],1]],
	\[Rho]=WeylVector[A,"DYN"],Ai=Inverse[A],aa=Array[Subscript["\[Alpha]", #]&,Length@A],k,DebugPrint},
			DebugPrint=If[debug,Print[##]&,Null&];
			G=Ai.DiagonalMatrix[SimpleRootLengths[A,"norm"->None]]/2;
			mult[M]=1;
			DebugPrint["M=",M," : T=",T," : G=",G];
			Do[Do[mult[\[Lambda]]=((M+\[Lambda]+2\[Rho]).G.(M-\[Lambda]))^(-1);num=0;
				DebugPrint[Style[Row[{"\[Lambda]=",\[Lambda]}],Bold],": M+\[Lambda]+2\[Rho]=",M+\[Lambda]+2\[Rho],": M-\[Lambda]=",((M-\[Lambda]).Ai).aa,": denom = ",mult[\[Lambda]]];
				Table[(k=1;up={};While[MemberQ[Flatten[wts,1],\[Lambda]+k*\[Alpha]]&&k<layer,AppendTo[up,\[Lambda]+k*\[Alpha]];k++]);
					DebugPrint["\[Alpha]=",\[Alpha],"-chain=",up];
					num+=Sum[mult[x](x.G.\[Alpha]),{x,up}],{\[Alpha],pRts}];
				mult[\[Lambda]]=2*num mult[\[Lambda]];
				DebugPrint["\[Lambda]=",\[Lambda],": numer ",num," : mult = ",mult[\[Lambda]]],
			{\[Lambda],wts[[layer]]}],{layer,2,T}];
			Map[mult,wts,{2}]]
WeightMults[type_?CartanBlocksQ,wts_,debug_:False] := WeightMults[CartanMatrix[type],wts,debug]


(* ::Subsection::Closed:: *)
(*WeylDimension*)


(* ::Text:: *)
(*Highest weight must be given in Dynkin rep.  Positive roots are by default in SRS rep. So innerproduct is simple to calculate.*)


WeylDimension[A_?(ArrayQ[#,2,IntegerQ]&),hw_?VectorQ] := With[{srl=DiagonalMatrix[SimpleRootLengths[A]],\[Delta]=WeylVector[A,"DYN"]},
		With[{num=srl.(hw+\[Delta]),den=srl.\[Delta]},Product[\[Alpha].num/\[Alpha].den,{\[Alpha],Flatten[PositiveRoots[A],1]}]]]
WeylDimension[type_?CartanBlocksQ,hw_?VectorQ] := WeylDimension[CartanMatrix[type],hw]


(* ::Subsubtitle::Closed:: *)
(*Endings*)


(* ::Text:: *)
(*End Private :*)


End[];


(* ::Text:: *)
(*Protect all of my lovely public functions*)


Protect[ZeroMatrix,StandardBasis,CartanMatrix,CartanMatrixQ,ConnectedComponents,
		DynkinDiagram,SimpleRootLengths,ClearRootCaches,RootInnerProduct,PositiveRoots];


(* ::Text:: *)
(*End the package:*)


EndPackage[];
