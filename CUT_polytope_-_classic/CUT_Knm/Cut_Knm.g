ListListSize:=[[1,3], [1,4], [1,5], 
               [2,2], [2,3], [2,4], [2,5], [2,6], [2,7], [2,8],
               [3,3], [3,4], [3,5], [3,6],
               [4,4], [4,5], [4,6], [4,7]];

ListUndone:=[];
for eListSize in ListListSize
do
  NameGraph:="K";
  for eVal in eListSize
  do
    NameGraph:=Concatenation(NameGraph, String(eVal));
  od;
  NameSave:=Concatenation("Data", NameGraph);
  if IsExistingFile(NameSave)=false then
    Add(ListUndone, eListSize);
  fi;
od;

for eListSize in ListUndone
do
  GRA:=GRAPH_GetMultiComplement(eListSize);
  NameGraph:="K";
  for eVal in eListSize
  do
    NameGraph:=Concatenation(NameGraph, String(eVal));
  od;
  #
  eRec:=CMC_GetCutPolytope(GRA);
  EXT:=eRec.ListVect;
  GRP:=LinPolytope_Automorphism(EXT);
  #
  ThePath:=Concatenation("SAVE_", NameGraph, "/");
  Exec("mkdir -p ", ThePath);
  #
  FuncStabilizer:=LinPolytope_Automorphism;
  FuncIsomorphy:=function(EXT1, EXT2)
    local ListInc1, ListInc2, ePerm;
    if Length(EXT1)<>Length(EXT2) then
      return false;
    fi;
    #
    ListInc1:=Set(List(EXT1, x->Position(EXT, x)));
    ListInc2:=Set(List(EXT2, x->Position(EXT, x)));
    ePerm:=RepresentativeAction(GRP, ListInc1, ListInc2, OnSets);
    if ePerm<>fail then
      return PermList(List(ListInc1, x->Position(ListInc2, OnPoints(x, ePerm))));
    fi;
    #
    return LinPolytope_Isomorphism(EXT1, EXT2);
  end;
  #
  FuncInvariant:=LinPolytope_Invariant;
  #
  IsRespawn:=function(OrdGRP, EXT, TheDepth)
    local nbInc;
    nbInc:=Length(EXT);
    if nbInc <= 150 then
      return false;
    fi;
    return true;
  end;
  #
  LimitNbVert:=100;
  #
  IsBankSave:=function(EllapsedTime, OrdGRP, EXT, TheDepth)
    if TheDepth=0 then
      return false;
    fi;
    if EllapsedTime>=600 then
      return true;
    fi;
    if Length(EXT)<=70 then
      return false;
    fi;
    return true;
  end;
  #
  TheFunc:=__DualDescriptionCDD_Reduction;
  #
  GetInitialRays_SamplingFramework:=function(EXT,nb)
    if Length(EXT)>200 then
      return SamplingStandard(EXT);
    else
      return GetInitialRays_LinProg(EXT,nb);
    fi;
  end;
  #
  TestNeedMoreSymmetry:=function(EXT)
    if Length(EXT)>350 then
      return true;
    else
      return false;
    fi;
  end;
  #
  Data:=rec(TheDepth:=0,
          ThePath:=ThePath,
          GetInitialRays:=GetInitialRays_SamplingFramework, 
          IsBankSave:=IsBankSave,
          GroupFormalism:=OnSetsGroupFormalism(LimitNbVert), 
          DualDescriptionFunction:=TheFunc, 
          TestNeedMoreSymmetry:=TestNeedMoreSymmetry,
          IsRespawn:=IsRespawn,
          Saving:=false,
          ThePathSave:="/irrelevant/");
  DataBank:=rec(BankPath:="/irrelevant/", Saving:=false);
  BF:=BankRecording(DataBank, FuncStabilizer, FuncIsomorphy, FuncInvariant, OnSetsGroupFormalism(LimitNbVert));
  LORB:=__ListFacetByAdjacencyDecompositionMethod(EXT, GRP, Data, BF);
  #
  fRec:=rec(EXT:=EXT, GRP:=GRP, LOrb:=LORB, eRec:=eRec);
  NameSave:=Concatenation("Data", NameGraph);
  SaveDataToFile(NameSave, fRec);
od;



