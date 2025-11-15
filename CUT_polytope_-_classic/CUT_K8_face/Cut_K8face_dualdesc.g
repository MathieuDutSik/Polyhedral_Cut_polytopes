#EXT:=ReadAsFunction("EXT64.gap")();
EXT:=ReadAsFunction("EXT60.gap")();
GRP:=LinPolytope_Automorphism(EXT);

FuncStabilizer:=LinPolytope_Automorphism;
FuncIsomorphy:=LinPolytope_Isomorphism;
FuncInvariant:=LinPolytope_Invariant;

IsRespawn:=function(OrdGRP, EXT, TheDepth)
  local nbInc;
  nbInc:=Length(EXT);
  if nbInc <= 59 then
    return false;
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
  if Length(EXT)<=45 then
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

Data:=rec(TheDepth:=0,
          ThePath:=ThePath,
          GetInitialRays:=GetInitialRays_SamplingFramework,
          IsBankSave:=IsBankSave,
          GroupFormalism:=OnSetsGroupFormalism(), 
          DualDescriptionFunction:=TheFunc, 
          IsRespawn:=IsRespawn,
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
