#
#
# G is a permutation group acting on [1..n], usually Sym(n)
MSN_Group:=function(G, m, n)
   local SS, SS2, dim, Gens, GenNew, eGen, UG, iCol;
   SS:=Combinations([1..n], m+1);
   SS2:=[[1..n]];
   for iCol in [1..Length(SS)]
   do
     SS2[iCol+1]:=SS[iCol];
   od;
   return TranslateGroup(G, SS2, OnSets);
end;



#
#
# creation of the facets
MSN_Facets:=function(m,s,n)
   local SS, dim, SSN, ListFacet, eSSN, iCol, NEW, jCol, Pos;
   SS:=Combinations([1..n], m+1);
   dim:=Length(SS);
   SSN:=Combinations([1..n], m+2);
   ListFacet:=[];
   for eSSN in SSN
   do
     for iCol in [1..m+2]
     do
	NEW:=NullZero(1+Length(SS));
	for jCol in [1..m+2]
	do
	  Pos:=Position(SS, Difference(eSSN, [eSSN[jCol]]));
	  NEW[Pos+1]:=1;
	od;
	Pos:=Position(SS, Difference(eSSN, [eSSN[iCol]]));
	NEW[Pos+1]:=-s;
	ListFacet[Length(ListFacet)+1]:=NEW;
     od;
   od;
   if s<m then
     for iCol in [1..Length(SS)]
     do
	NEW:=NullZero(1+Length(SS));
	NEW[iCol+1]:=1;
	ListFacet[Length(ListFacet)+1]:=NEW;
     od;
   fi;
   return ListFacet;
end;


#
#
# n must be equal to m+3
ExtRay012:=function(Part, m, n)
  local SS, Tot, NEW, iCol, iPart, Pos, S, iS;
  SS:=Combinations([1..n], m+1);
  Tot:=[1..n];
  NEW:=[0];
  for iCol in [1..Length(SS)]
  do
    NEW[iCol+1]:=2;
  od;
  for iPart in [1..Length(Part)]
  do
    if Length(Part[iPart])=2 then
      Pos:=Position(SS, Difference(Tot, Part[iPart]));
      NEW[Pos+1]:=0;
    else
      S:=Combinations(Part[iPart], 2);
      for iS in [1..Length(S)]
      do
	Pos:=Position(SS, Difference(Tot, S[iS]));
	NEW[Pos+1]:=1;
      od;
    fi;
  od;
  return NEW;
end;
