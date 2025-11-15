FCT__CutPolytope:=function(n, DataPolyhedral)
  local LSET, EXT, eSet, V, i, j, H, DistMat, eList, dist, FuncDistMat, FuncStabilizer, FuncIsomorphy, DD, BF, h, DataBank, GRP, FuncInvariant;
  LSET:=Combinations([1..n-1]);
  #
  EXT:=[];
  for eSet in LSET
  do
    V:=[1];
    for i in [1..n-1]
    do
      for j in [i+1..n]
      do
        H:=Intersection(eSet, [i,j]);
        if Length(H)=1 then
          Add(V, 1);
        else
          Add(V, 0);
        fi;
      od;
    od;
    Add(EXT, V);
  od;
  #
  DistMat:=[];
  for i in [1..Length(EXT)]
  do
    eList:=[];
    for j in [1..Length(EXT)]
    do
      dist:=0;
      for h in [2..Length(EXT[1])]
      do
        if EXT[i][h]<>EXT[j][h] then
          dist:=dist+1;
        fi;
      od;
      Add(eList, dist);
    od;
    Add(DistMat, eList);
  od;
  #
  FuncDistMat:=function(ListInc)
    return List(ListInc, x->DistMat[x]{ListInc});
  end;
  FuncStabilizer:=function(EXTask)
    local ListInc;
    ListInc:=List(EXTask, x->Position(EXT, x));
    if Length(ListInc)<=30 then
      return Group(());
    fi;
    return AutomorphismGroupEdgeColoredGraph(FuncDistMat(ListInc));
  end;
  FuncIsomorphy:=function(EXT1, EXT2)
    local test, ListInc1, ListInc2;
    ListInc1:=List(EXT1, x->Position(EXT, x));
    ListInc2:=List(EXT2, x->Position(EXT, x));
    test:=IsIsomorphicEdgeColoredGraph(FuncDistMat(ListInc1), FuncDistMat(ListInc2));
    if test=false then
      return false;
    else
      return PermList(test);
    fi;
  end;
  FuncInvariant:=function(EXT)
    return [Length(EXT), RankMat(EXT)];
  end;
  DataBank:=rec(BankPath:=DataPolyhedral.ThePath, Saving:=true);
  BF:=BankRecording(DataBank, FuncStabilizer, FuncIsomorphy, FuncInvariant, OnSetsGroupFormalism());
  #
  GRP:=AutomorphismGroupEdgeColoredGraph(FuncDistMat([1..Length(EXT)]));
  DD:=function()
    return __ListFacetByAdjacencyDecompositionMethod(EXT, GRP, DataPolyhedral, BF);
  end;
  Exec("mkdir -p ", DataPolyhedral.ThePath);
  return rec(EXT:=EXT, GRP:=GRP, DD:=DD, FuncDistMat:=FuncDistMat, FuncStabilizer:=FuncStabilizer, FuncIsomorphy:=FuncIsomorphy);
end;
