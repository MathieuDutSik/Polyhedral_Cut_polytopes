ListSize:=[4,7];
GRA:=GRAPH_GetMultiComplement(ListSize);

eRec:=CMC_GetCutPolytope(GRA);
EXT:=eRec.ListVect;
GRP:=LinPolytope_Automorphism(EXT);

FuncStabilizer:=LinPolytope_Automorphism;
FuncIsomorphy:=LinPolytope_Isomorphism;
FuncInvariant:=LinPolytope_Invariant;

IsRespawn:=function(OrdGRP, EXT, TheDepth)
  local nbInc;
  nbInc:=Length(EXT);
  if nbInc <= 85 then
    return false;
  fi;
  if nbInc <= 105 then
    if OrdGRP>10 then
      return true;
    else
      return false;
    fi;
  fi;
  return true;
end;

TestNeedMoreSymmetry:=function(EXT)
  if Length(EXT)<150 then
    return false;
  else
    return true;
  fi;
end;


IsBankSave:=function(EllapsedTime, OrdGRP, EXT, TheDepth)
  if TheDepth=0 then
    return false;
  fi;
  if EllapsedTime>=600 then
    return true;
  fi;
  if Length(EXT)<=90 then
    return false;
  fi;
  return true;
end;
ThePathWork:="./TheWork/";
Exec("mkdir -p ", ThePathWork);
ThePath:=Concatenation(ThePathWork, "tmp/");
Exec("mkdir -p ", ThePath);
PathSave:=Concatenation(ThePathWork, "PathSAVE/");
Exec("mkdir -p ", PathSave);

#TheFunc:=__DualDescriptionLRS_Reduction;
TheFunc:=__DualDescriptionCDD_Reduction;

Data:=rec(TheDepth:=0,
          ThePath:=ThePath,
          GetInitialRays:=GetInitialRays_LinProg, 
          IsBankSave:=IsBankSave,
          GroupFormalism:=OnSetsGroupFormalism(), 
          DualDescriptionFunction:=TheFunc, 
          IsRespawn:=IsRespawn,
          TestNeedMoreSymmetry:=TestNeedMoreSymmetry, 
          Saving:=true,
          ThePathSave:=PathSave);

BankPath:=Concatenation(ThePathWork, "TheBank/");
Exec("mkdir -p ", BankPath);
DataBank:=rec(BankPath:=BankPath, Saving:=true);
BF:=BankRecording(DataBank, FuncStabilizer, FuncIsomorphy, FuncInvariant, OnSetsGroupFormalism());

LORB:=__ListFacetByAdjacencyDecompositionMethod(EXT,
                                                GRP, 
                                                Data,
                                                BF);
SaveDataToFile("ListOrbitEXT", LORB);
