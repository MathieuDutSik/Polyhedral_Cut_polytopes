ListSize:=[1,1,1,1,1];
GRA:=GRAPH_GetMultiComplement(ListSize);

eRec:=CMC_GetCutPolytope(GRA);
EXT:=eRec.GetCUT_info().EXT;

FileSaveGRP:="TheData_GRP";
if IsExistingFile(FileSaveGRP) then
  GRP:=ReadAsFunction(FileSaveGRP)();
else
  GRP:=LinPolytope_Automorphism(EXT);
  SaveDataToFile(FileSaveGRP, GRP);
fi;
Print("We have GRP\n");



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


LinPolytope_Isomorphism;
FuncInvariant:=LinPolytope_InvariantMD5;

IsRespawn:=function(OrdGRP, EXT, TheDepth)
  local nbInc;
  nbInc:=Length(EXT);
  if nbInc <= 280 then
    return false;
  fi;
  if nbInc <= 335 then
    if OrdGRP >= 16 then
      return true;
    else
      return false;
    fi;
  fi;
  return true;
end;


IsBankSave:=function(EllapsedTime, OrdGRP, EXT, TheDepth)
  if TheDepth=0 then
    return false;
  fi;
  if EllapsedTime>=600 then
    return true;
  fi;
  if Length(EXT)<=240 then
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

GetInitialRays_SamplingFramework:=function(EXT,nb)
  if Length(EXT)>200 then
    return SamplingStandard(EXT);
  else
    return GetInitialRays_LinProg(EXT,nb);
  fi;
end;

TestNeedMoreSymmetry:=function(EXT)
  if Length(EXT)>350 then
    return true;
  else
    return false;
  fi;
end;


Data:=rec(TheDepth:=0,
          ThePath:=ThePath,
          GetInitialRays:=GetInitialRays_SamplingFramework, 
          IsBankSave:=IsBankSave,
          GroupFormalism:=OnSetsGroupFormalism(500), 
          DualDescriptionFunction:=TheFunc, 
          TestNeedMoreSymmetry:=TestNeedMoreSymmetry,
          IsRespawn:=IsRespawn,
          Saving:=true,
          ThePathSave:=PathSave);

BankPath:=Concatenation(ThePathWork, "TheBank/");
Exec("mkdir -p ", BankPath);
DataBank:=rec(BankPath:=BankPath, Saving:=true);
BF:=BankRecording(DataBank, FuncStabilizer, FuncIsomorphy, FuncInvariant, OnSetsGroupFormalism(500));

LORB:=__ListFacetByAdjacencyDecompositionMethod(EXT,
                                                GRP, 
                                                Data,
                                                BF);
SaveDataToFile("ListOrbitEXT", LORB);
