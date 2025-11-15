Read("Functions/MyPersoFunctions.g");
#
EXT56:=ReadAsFunction("Functions/TheGosset")();;
for i in [1..56]
do
  EXT56[i][1]:=1;
od;
DM56:=FuncDistMat(EXT56);
GRP56:=AutomorphismGroupEdgeColoredGraph(DM56);
DDA56:=DualDescriptionAdjacencies(EXT56);
GRA56:=DDA56.SkeletonGraph;
#
EXT27:=ReadAsFunction("Functions/TheSchlafli")();;
DM27:=FuncDistMat(EXT27);
GRP27:=AutomorphismGroupEdgeColoredGraph(DM27);
DDA27:=DualDescriptionAdjacencies(EXT27);
GRA27:=DDA27.SkeletonGraph;






MyFindOrbitSubGraphsNautyMethod:=function(nbv, TheGRA, FuncSelect, AskedSize)
  local GRP, O, Cand, eOrb, iP, CandSec, Stab, O2, eO2, eCand, FuncInsert;
  #
  GRP:=AutGroupGraph(TheGRA);
  O:=Orbits(GRP, [1..nbv], OnPoints);
  Cand:=[];
  for eOrb in O
  do
    AddSet(Cand, [Minimum(eOrb)]);
  od;
  Print("Find ",Length(O)," orbits at step 1\n");
  #
  for iP in [2..AskedSize]
  do
    CandSec:=[];
    for eCand in Cand
    do
      Stab:=Stabilizer(GRP, eCand, OnSets);
      O2:=Orbits(Stab, Difference([1..nbv], eCand), OnPoints);
      for eO2 in O2
      do
	AddSet(CandSec, Union(eCand, [eO2[1]]));
      od;
    od;
    Cand:=[];
    for eCand in CandSec
    do
      if FuncSelect(eCand)=true then
	AddSet(Cand, eCand);
      fi;
    od;
    CandSec:=[];
    FuncInsert:=function(eCand)
      local fCand, TheCanon;
      TheCanon:=CharacteristicGraphOfSubset(TheGRA, eCand);
      for fCand in CandSec
      do
        if fCand.TheCanon=TheCanon then
          return;
        fi;
      od;
      Add(CandSec, rec(eCand:=eCand, TheCanon:=TheCanon));
    end;
    for eCand in Cand
    do
      FuncInsert(eCand);
    od;
    Cand:=List(CandSec, x->x.eCand);
    Print("Find ", Length(Cand)," orbits at step ", iP, "\n");
  od;
  return Cand;
end;




MyFindOrbitSubGraphs:=function(nbv, TheGroup, FuncSelect, AskedSize)
  local O, Cand, eOrb, iP, CandSec, Stab, O2, eO2, eCand, FuncInsert;
  #
  O:=Orbits(TheGroup, [1..nbv], OnPoints);
  Cand:=[];
  for eOrb in O
  do
    AddSet(Cand, [Minimum(eOrb)]);
  od;
  Print("Find ",Length(O)," orbits at step 1\n");
  #
  for iP in [2..AskedSize]
  do
    CandSec:=[];
    for eCand in Cand
    do
      Stab:=Stabilizer(TheGroup, eCand, OnSets);
      O2:=Orbits(Stab, Difference([1..nbv], eCand), OnPoints);
      for eO2 in O2
      do
	AddSet(CandSec, Union(eCand, [eO2[1]]));
      od;
    od;
    Cand:=[];
    for eCand in CandSec
    do
      if FuncSelect(eCand)=true then
	AddSet(Cand, eCand);
      fi;
    od;
    CandSec:=[];
    FuncInsert:=function(eCand)
      local fCand;
      for fCand in CandSec
      do
        if RepresentativeAction(TheGroup, eCand, fCand, OnSets)<>fail then
          return;
        fi;
      od;
      Add(CandSec, eCand);
    end;
    for eCand in Cand
    do
      FuncInsert(eCand);
    od;
    Cand:=CandSec;
    Print("Find ",Length(Cand)," orbits at step ", iP, "\n");
  od;
  return Cand;
end;




MyFindAllOrbitAffineBasis:=function(Polytope, TheGroup, TheGraph)
  local Dimension, FuncLinear;
  Dimension:=RankMat(Polytope);
  FuncLinear:=function(eCand)
    local Matr, ePos, iVert, Solu, eV;
    Matr:=[];
    for ePos in eCand
    do
      Add(Matr, Polytope[ePos]);
    od;
    if RankMat(Matr)<Length(eCand) then
      return false;
    fi;
    for iVert in Difference([1..Length(Polytope)], eCand)
    do
      Solu:=SolutionMat(Matr, Polytope[iVert]);
      if Solu<>fail then
        for eV in Solu
        do
          if IsInt(eV)=false then
            return false;
          fi;
        od;
      fi;
    od;
    return true;
  end;
#  return MyFindOrbitSubGraphs(Length(Polytope), TheGroup, FuncLinear, Dimension);
  return MyFindOrbitSubGraphsNautyMethod(Length(Polytope), TheGraph, FuncLinear, Dimension);
end;

#
#
# find the simplices which can be separated from the rest of
# the polytope.
MyFindAllOrbitBaranovskiSimplex:=function(Polytope, TheGroup, TheGraph)
  local Dimension, GetValueVector, FuncSelection, FuncOK, ErdahlConeEXT, ErdahlConeFAC, ListORB, ListEmpty, ListNonEmpty, eCand, eVect, TheDim;
  Dimension:=RankMat(Polytope);
  GetValueVector:=function(eVect)
    return SymmetricMatrixToVector(TransposedMat([eVect])*[eVect]);
  end;
  ErdahlConeEXT:=ColumnReduction(List(Polytope, GetValueVector)).EXT;
  ErdahlConeFAC:=DualDescription(ErdahlConeEXT);
  Print("|FAC ErdahlCone|=", Length(ErdahlConeFAC), "\n");
  FuncSelection:=function(eCand)
    local Matr;
    Matr:=Polytope{eCand};
    return RankMat(Matr)=Length(eCand);
  end;
  FuncOK:=function(eCand)
    local TheCentralEXT, TheCentralFAC, RelevantCand;
    TheCentralEXT:=Sum(ErdahlConeEXT{eCand});
    TheCentralFAC:=Sum(Filtered(ErdahlConeFAC, x->TheCentralEXT*x=0));
    RelevantCand:=Filtered([1..Length(Polytope)], x->ErdahlConeEXT[x]*TheCentralFAC=0);
#    Print("RelevantCand=", RelevantCand, "  eCand=", eCand, "\n");
    return RelevantCand=eCand;
  end;
#  return MyFindOrbitSubGraphs(Length(Polytope), TheGroup, FuncSelection, Dimension);
  ListORB:=MyFindOrbitSubGraphsNautyMethod(Length(Polytope), TheGraph, FuncSelection, Dimension);
  ListNonEmpty:=[];
  ListEmpty:=[];
  for eCand in ListORB
  do
    if FuncOK(eCand)=true then
      Add(ListNonEmpty, eCand);
    else
      Add(ListEmpty, eCand);
    fi;
  od;
  return rec(ListEmpty:=ListEmpty, ListNonEmpty:=ListNonEmpty);
end;


#ListOrbSchlafli:=MyFindAllOrbitBaranovskiSimplex(EXT27, GRP27, GRA27);
ListOrbGosset:=MyFindAllOrbitBaranovskiSimplex(EXT56, GRP56, GRA56);
SaveDataToFile("ListOrbIndepGosset", ListOrbGosset);

#ListOrbGosset:=MyFindAllOrbitAffineBasis(EXT56, GRP56, GRA56);
#ListOrbSchlafli:=MyFindAllOrbitAffineBasis(EXT27, GRP27, GRA27);
