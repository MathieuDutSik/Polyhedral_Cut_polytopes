OnPairMapping:=function(n, TheGRP)
  local ListSymb, i, j, ListPermGens, eGen, eList, eSymb;
  ListSymb:=[];
  for i in [1..n-1]
  do
    for j in [i+1..n]
    do
      Add(ListSymb, [i,j]);
    od;
  od;
  ListPermGens:=[];
  for eGen in GeneratorsOfGroup(TheGRP)
  do
    eList:=[1];
    for eSymb in ListSymb
    do
      Add(eList, 1+Position(ListSymb, OnSets(eSymb, eGen)));
    od;
    Add(ListPermGens, PermList(eList));
  od;
  return Group(ListPermGens);
end;


SymmetricGroupOnPairs:=function(n)
  return OnPairMapping(n, SymmetricGroup(n));
end;




__PermutedStabilizer:=function(GRP, eVect)
  local H, TheStab, eVal, Stab1, FC;
  H:=Set(eVect);
  TheStab:=Group(GeneratorsOfGroup(GRP));
  for eVal in H
  do
    FC:=Filtered([1..Length(eVect)], x->eVect[x]=eVal);
    Stab1:=Stabilizer(TheStab, FC, OnSets);
    TheStab:=PersoGroupPerm(SmallGeneratingSet(Stab1));
  od;
  return TheStab;
end;



__PermutedRepresentativeAction:=function(GRP, eVect1, eVect2)
  local H1, H2, eVectSecond1, H, gSearch, TheGrp, eVal, FC1, FC2, g, GrpStab;
  H1:=Collected(eVect1);
  H2:=Collected(eVect2);
  eVectSecond1:=ShallowCopy(eVect1);
  if H1<>H2 then
    return fail;
  fi;
  H:=Set(eVect1);
  gSearch:=();
  TheGrp:=Group(GeneratorsOfGroup(GRP));
  for eVal in H
  do
    FC1:=Filtered([1..Length(eVect1)], x->eVectSecond1[x]=eVal);
    FC2:=Filtered([1..Length(eVect2)], x->eVect2[x]=eVal);
    g:=RepresentativeAction(TheGrp, FC1, FC2, OnSets);
    if g=fail then
      return fail;
    fi;
    gSearch:=gSearch*g;
    eVectSecond1:=Permuted(eVectSecond1, g);
    if eVectSecond1=eVect2 then
      return gSearch;
    fi;
    GrpStab:=Stabilizer(TheGrp, eVectSecond1, OnSets);
    TheGrp:=Group(SmallGeneratingSet(GrpStab));
  od;
end;
