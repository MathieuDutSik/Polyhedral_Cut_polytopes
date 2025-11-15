Read("Functions/MyPersoFunctions.g");
#
GramMat_E6:=[[2,1,1,0,1,1],[1,4,1,1,1,3],[1,1,2,1,1,1],
[0,1,1,2,1,2],[1,1,1,1,2,2],[1,3,1,2,2,4]];

GRP_E6:=ArithmeticAutomorphismMatrixFamily_Souvignier("", [GramMat_E6], []);
DelCO_E6:=DelaunayDescriptionStandard(rec(GramMat:=GramMat_E6, MatrixGRP:=GRP_E6));
EXT_read:=DelCO_E6[1].EXT;
TheOrig:=ListWithIdenticalEntries(Length(EXT_read[1]),0);
TheOrig[1]:=1;
TheTrans:=TheOrig-EXT_read[1];
EXT27:=List(EXT_read, x->x+TheTrans);


GramMat_ER7:=ClassicalSporadicLattices("ER7");
GRP_ER7:=ArithmeticAutomorphismMatrixFamily_Souvignier("", [GramMat_ER7], []);
DelCO_ER7:=DelaunayDescriptionStandard(rec(GramMat:=GramMat_ER7, MatrixGRP:=GRP_ER7));
ListSizes:=List(DelCO_ER7, x->Length(x.EXT));
Pos:=Position(ListSizes, 35);
EXT_read:=DelCO_ER7[Pos].EXT;
TheOrig:=ListWithIdenticalEntries(Length(EXT_read[1]),0);
TheOrig[1]:=1;
TheTrans:=TheOrig-EXT_read[1];
EXT35:=List(EXT_read, x->x+TheTrans);

DM27:=FuncDistMat(EXT27);
GRP27:=AutomorphismGroupEdgeColoredGraph(DM27);
DDA27:=DualDescriptionAdjacencies(EXT27);
GRA27:=DDA27.SkeletonGraph;
#
DM35:=FuncDistMat(EXT35);
GRP35:=AutomorphismGroupEdgeColoredGraph(DM35);
DDA35:=DualDescriptionAdjacencies(EXT35);
GRA35:=DDA35.SkeletonGraph;




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




MyFindOrbitSubGraphsMinimumOrbitMethod:=function(nbv, TheGroup, FuncSelect, AskedSize)
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
      local fCand, TheCanon;
      TheCanon:=Minimum(Orbit(TheGroup, eCand, OnSets));
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




MyFindOrbitSubGraphsBacktrackMethod:=function(nbv, TheGroup, FuncSelect, AskedSize)
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





#
# find the simplices which can be separated from the rest of
# the polytope.
MyFindAllOrbitBaranovskiSimplex:=function(Polytope, TheGroup, TheGraph)
  local Dimension, GetValueVector, FuncSelection, FuncOK, ErdahlConeEXT, ErdahlConeFAC, ListORB, ListEmpty, ListNonEmpty, eCand;
  Dimension:=RankMat(Polytope);
  GetValueVector:=function(eVect)
    return SymmetricMatrixToVector(TransposedMat([eVect])*[eVect]);
  end;
  ErdahlConeEXT:=List(Polytope, GetValueVector);
  ErdahlConeFAC:=DualDescription(ErdahlConeEXT);
  Print("|FAC ErdahlCone|=", Length(ErdahlConeFAC), "\n");
  FuncSelection:=function(eCand)
    local Matr;
    Matr:=Polytope{eCand};
    return RankMat(Matr)=Length(eCand);
  end;
  FuncOK:=function(eCand)
    local TheCentralEXT, TheCentralFAC, RelevantCand;
    TheCentralEXT:=Sum(List(Polytope{eCand}, GetValueVector));
    TheCentralFAC:=Sum(Filtered(ErdahlConeFAC, x->TheCentralEXT*x=0));
    RelevantCand:=Filtered([1..Length(Polytope)], x->ErdahlConeEXT[x]*TheCentralFAC=0);
#    Print("RelevantCand=", RelevantCand, "  eCand=", eCand, "\n");
    return RelevantCand=eCand;
  end;
#  ListORB:=MyFindOrbitSubGraphsBacktrackMethod(Length(Polytope), TheGroup, FuncSelection, Dimension);
#  ListORB:=MyFindOrbitSubGraphsNautyMethod(Length(Polytope), TheGraph, FuncSelection, Dimension);
  ListORB:=MyFindOrbitSubGraphsMinimumOrbitMethod(Length(Polytope), TheGroup, FuncSelection, Dimension);
  Print("We are here 1\n");
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
  Print("We are here 2\n");
  return rec(ListEmpty:=ListEmpty, ListNonEmpty:=ListNonEmpty);
end;


ListOrbSchlafli:=MyFindAllOrbitBaranovskiSimplex(EXT27, GRP27, GRA27);
ListOrb35:=MyFindAllOrbitBaranovskiSimplex(EXT35, GRP35, GRA35);

#SaveDataToFile("BaraCone27", rec(GramMat_E6:=GramMat_E6, EXT27:=EXT27, ListOrbit:=ListOrbSchlafli));
SaveDataToFile("BaraCone35", rec(GramMat_ER7:=GramMat_ER7, EXT35:=EXT35, ListOrbit:=ListOrb35));
