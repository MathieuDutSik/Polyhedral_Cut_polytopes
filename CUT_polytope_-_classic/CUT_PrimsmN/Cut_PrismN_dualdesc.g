GetDualDescOfCutPolytopeGraph:=function(EXT, GRP)
  local FuncStabilizer, FuncIsomorphy, FuncInvariant, IsRespawn, IsBankSave, TheFunc, Data, DataBank, BF, LORB, ThePath, PathSave, BankPath;
  FuncStabilizer:=LinPolytope_Automorphism;
  FuncIsomorphy:=LinPolytope_Isomorphism;
  FuncInvariant:=LinPolytope_Invariant;

  IsRespawn:=function(OrdGRP, EXT, TheDepth)
    local nbInc;
    nbInc:=Length(EXT);
    if nbInc <= 70 then
      return false;
    fi;
    return true;
  end;
  IsBankSave:=function(EllapsedTime, OrdGRP, EXT, TheDepth)
    if TheDepth=0 then
      return false;
    fi;
    if EllapsedTime>=60 then
      return true;
    fi;
    return false;
  end;
  #TheFunc:=__DualDescriptionLRS_Reduction;
  TheFunc:=__DualDescriptionCDD_Reduction;

  ThePath:=Filename(POLYHEDRAL_tmpdir, "");
  PathSave:="./irrelevant/";
  BankPath:="./irrelevant/";
  Data:=rec(TheDepth:=0,
            ThePath:=ThePath,
            GetInitialRays:=GetInitialRays_LinProg, 
            IsBankSave:=IsBankSave,
            GroupFormalism:=OnSetsGroupFormalism(), 
            DualDescriptionFunction:=TheFunc, 
            IsRespawn:=IsRespawn,
            Saving:=false,
            ThePathSave:=PathSave);

  DataBank:=rec(BankPath:=BankPath, Saving:=false);
  BF:=BankRecording(DataBank, FuncStabilizer, FuncIsomorphy, FuncInvariant, OnSetsGroupFormalism());
  LORB:=__ListFacetByAdjacencyDecompositionMethod(EXT,
                                                  GRP, 
                                                  Data,
                                                  BF);
  return LORB;
end;

GetDualDescOfCutPrism:=function(N)
  local ePL, GRA, eRec, EXT, GRP, LORB, eFile, fRec;
  ePL:=Prism(N);
  GRA:=PlanGraphToGRAPE(ePL);
  eRec:=CMC_GetCutPolytope(GRA);
  EXT:=eRec.ListVect;
  GRP:=LinPolytope_Automorphism(EXT);
  LORB:=GetDualDescOfCutPolytopeGraph(EXT, GRP);
  fRec:=rec(eRec:=eRec, 
            EXT:=eRec.ListVect,
            LOrb:=LORB, 
            GRP:=GRP);
  eFile:=Concatenation("DataPrism", String(N));
  SaveDataToFile(eFile, fRec);
end;

GetDualDescOfCutAntiPrism:=function(N)
  local ePL, GRA, eRec, EXT, GRP, LORB, eFile, fRec;
  ePL:=AntiPrism(N);
  GRA:=PlanGraphToGRAPE(ePL);
  eRec:=CMC_GetCutPolytope(GRA);
  EXT:=eRec.ListVect;
  GRP:=LinPolytope_Automorphism(EXT);
  LORB:=GetDualDescOfCutPolytopeGraph(EXT, GRP);
  fRec:=rec(eRec:=eRec, 
            EXT:=eRec.ListVect,
            LOrb:=LORB, 
            GRP:=GRP);
  eFile:=Concatenation("DataAntiPrism", String(N));
  SaveDataToFile(eFile, fRec);
end;
