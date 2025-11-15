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


GramMat_E7:=[[4,3,2,3,1,1,1],[3,4,3,3,1,1,1],[2,3,4,3,2,2,1],
[3,3,3,4,2,1,2],[1,1,2,2,2,1,1],[1,1,2,1,1,2,0],[1,1,1,2,1,0,2]];
GRP_E7:=ArithmeticAutomorphismMatrixFamily_Souvignier("", [GramMat_E7], []);
DelCO_E7:=DelaunayDescriptionStandard(rec(GramMat:=GramMat_E7, MatrixGRP:=GRP_E7));
ListSizes:=List(DelCO_E7, x->Length(x.EXT));
Pos:=Position(ListSizes, 56);
EXT_read:=DelCO_E7[Pos].EXT;
TheOrig:=ListWithIdenticalEntries(Length(EXT_read[1]),0);
TheOrig[1]:=1;
TheTrans:=TheOrig-EXT_read[1];
EXT56:=List(EXT_read, x->x+TheTrans);

DM27:=FuncDistMat(EXT27);
GRP27:=AutomorphismGroupEdgeColoredGraph(DM27);
DDA27:=DualDescriptionAdjacencies(EXT27);
GRA27:=DDA27.SkeletonGraph;
#
DM56:=FuncDistMat(EXT56);
GRP56:=AutomorphismGroupEdgeColoredGraph(DM56);
DDA56:=DualDescriptionAdjacencies(EXT56);
GRA56:=DDA56.SkeletonGraph;



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


TestAchillCriterion:=function(EXT, GramMat, ListOrbitSimplex)
  local eOrb, test, TheTest;
  TheTest:="relative interior";
#  TheTest:="closed";
  for eOrb in ListOrbitSimplex.ListNonEmpty
  do
    Print("eOrb=", eOrb, "\n");
    test:=TestSimplexPropertyInExtreme(EXT, GramMat, eOrb, TheTest);
    if test=false then
      return false;
    fi;
  od;
  return true;
end;




ListOrbSchlafli:=MyFindAllOrbitBaranovskiSimplex(EXT27, GRP27, GRA27);
SaveDataToFile("BaraCone27", rec(GramMat_E6:=GramMat_E6, EXT27:=EXT27, ListOrbit:=ListOrbSchlafli));
test27:=TestAchillCriterion(EXT27, GramMat_E6, ListOrbSchlafli);


ListOrbGosset:=MyFindAllOrbitBaranovskiSimplex(EXT56, GRP56, GRA56);
SaveDataToFile("BaraCone56", rec(GramMat_E7:=GramMat_E7, EXT56:=EXT56, ListOrbit:=ListOrbGosset));
test56:=TestAchillCriterion(EXT56, GramMat_E7, ListOrbGosset);









#ListOrbGosset:=MyFindAllOrbitAffineBasis(EXT56, GRP56, GRA56);
#ListOrbSchlafli:=MyFindAllOrbitAffineBasis(EXT27, GRP27, GRA27);
