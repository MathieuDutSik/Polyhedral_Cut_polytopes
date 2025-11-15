ListSize:=[4,4];
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
  if nbInc<=42 then
    return false;
  fi;
  if nbInc <= 50 then
    return false;
  fi;
  if nbInc>90 and OrdGRP>50 then
    return true;
  fi;
  if nbInc >=96 and OrdGRP >= 32 then
    return true;
  fi;
  if nbInc >= 108 and OrdGRP >= 16 then
    return true;
  fi;
  if nbInc=144 and TheDepth=2 then
    return true;
  fi;
  if nbInc=120 and TheDepth=3 then
    return true;
  fi;
  if TheDepth=3 then
    return false;
  fi;
  if TheDepth>=2 then
    return false;
  fi;
  Print("We are in IsRespawn\n");
  if OrdGRP>=200 then
    return true;
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
  if Length(EXT)<=27 then
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
