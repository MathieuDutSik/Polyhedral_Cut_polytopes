Read("FunctionSet.g");

n:=6;
eRec:=GetRecordForOperations(n);
GRP:=LinPolytope_Automorphism(eRec.EXT);

RecHeuristic:=rec(MaxNumberUndone:=1000,
                  NbAllow:=100,
                  MaxDepth:=4);

eFacMet:=eRec.FAC_Met[1];
eSetInc:=Filtered([1..Length(eRec.EXT)], x->eRec.EXT[x]*eFacMet=0);
ListOrbitToTest:=[eSetInc];
TestConn:=AdvancedBalinskiConnectivity_Standard(eRec.EXT, GRP, ListOrbitToTest, RecHeuristic);


DirectTest:=function()
  local DDA, eDiff, eSetDiff, GRA, i, eAdj, pos;
  DDA:=DualDescriptionAdjacencies(eRec.EXT);
  eDiff:=Difference(Set(DDA.FAC), Set(eRec.FAC_Met));
  eSetDiff:=List(eDiff, x->Position(DDA.FAC, x));
  GRA:=NullGraph(Group(()), Length(eSetDiff));
  for i in [1..Length(eSetDiff)]
  do
    for eAdj in Adjacency(DDA.RidgeGraph, eSetDiff[i])
    do
      pos:=Position(eSetDiff, eAdj);
      if pos<>fail then
        AddEdgeOrbit(GRA, [i, pos]);
      fi;
    od;
  od;
  return GRA;
end;
