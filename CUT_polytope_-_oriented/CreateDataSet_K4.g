TreatPolyhedralCase:=function(CASE_info, FileSaveFinal, ScratchDir)
  local EXTrel, GRPrel, ThePathWork, ThePath, PathSave, FuncStabilizer, FuncIsomorphy, FuncInvariant, IsRespawn, TestNeedMoreSymmetry, IsBankSave, TheFunc, LimitNbVert, Data, BankPath, DataBank, BF, LORB;
  EXTrel:=CASE_info.EXT;
  GRPrel:=CASE_info.GRP;
  #
  ThePathWork:=ScratchDir;
  Exec("mkdir -p ", ThePathWork);
  ThePath:=Concatenation(ThePathWork, "tmp/");
  Exec("mkdir -p ", ThePath);
  PathSave:=Concatenation(ThePathWork, "PathSAVE/");
  Exec("mkdir -p ", PathSave);
  #
  FuncStabilizer:=LinPolytope_Automorphism;
  FuncIsomorphy:=LinPolytope_Isomorphism;
  FuncInvariant:=LinPolytope_Invariant;
  #
  IsRespawn:=function(OrdGRP, EXT, TheDepth)
    local nbInc;
    nbInc:=Length(EXT);
    if nbInc <= 90 then
      return false;
    fi;
    return true;
  end;
  #
  TestNeedMoreSymmetry:=function(EXT)
    if Length(EXT)<=50 then
      return false;
    else
      return true;
    fi;
  end;
  #
  IsBankSave:=function(EllapsedTime, OrdGRP, EXT, TheDepth)
    if TheDepth=0 then
      return false;
    fi;
    if EllapsedTime>=600 then
      return true;
    fi;
    if Length(EXT)<=50 then
      return false;
    fi;
    return true;
  end;
  #
#   TheFunc:=__DualDescriptionCDD_Reduction;
  TheFunc:=__DualDescriptionPPL_Reduction;
  #
  LimitNbVert:=100;
  Data:=rec(TheDepth:=0,
            ThePath:=ThePath,
            GetInitialRays:=GetInitialRays_LinProg,
            IsBankSave:=IsBankSave,
            GroupFormalism:=OnSetsGroupFormalism(LimitNbVert),
            TestNeedMoreSymmetry:=TestNeedMoreSymmetry,
            DualDescriptionFunction:=TheFunc,
            IsRespawn:=IsRespawn,
            Saving:=true,
            ThePathSave:=PathSave);
  #
  BankPath:=Concatenation(ThePathWork, "TheBank/");
  Exec("mkdir -p ", BankPath);
  DataBank:=rec(BankPath:=BankPath, Saving:=true);
  BF:=BankRecording(DataBank, FuncStabilizer, FuncIsomorphy, FuncInvariant, OnSetsGroupFormalism(LimitNbVert));
  #
  if IsExistingFilePlusTouch(FileSaveFinal)=false then
    LORB:=__ListFacetByAdjacencyDecompositionMethod(EXTrel, GRPrel, Data, BF);
    SaveDataToFilePlusTouch(FileSaveFinal, LORB);
  else
    LORB:=ReadAsFunction(FileSaveFinal)();
  fi;
  return LORB;
end;




ListRec:=ReadAsFunction("ListClassicalGraph")();

DirScratch:="SCRATCH/";

MaxNumberRelEdge:=24;

#ListTargetGraph:=["Cube", "Octahedron", "Petersen", "Prism3", "Prism5", "Antiprism4", "MobiusLadder4", "MobiusLadder5", 
#                  "K113", "K114", "K1113", "K1122", "K1123", "K6-K2", "K7-K3", "K7-C4", "K7-C7", "K7-K4", "K8-K5",
#                  "PyrPrism3", "PyrPyrCycl5"];
ListTargetGraph:=["Octahedron", "Prism3"];


ListRecSel:=Filtered(ListRec, x->Position(ListTargetGraph, x.name)<>fail);

for eRec in ListRecSel
do
  Print("name=", eRec.name, "\n");
  eRecOri:=CMC_GetCutPolytopeOriented(eRec.GRA);
  if Length(eRecOri.eSelect) <= MaxNumberRelEdge then
    #
    eScratchMET:=Concatenation(DirScratch, "OCUT_MET_info_", eRec.name, "/");
    Exec("mkdir -p ", eScratchMET);
    FileSaveMET:=Concatenation("OCUT_MetPolytopes/OCUT_MET_dualdesc_", eRec.name);
    FileIndicate:=Concatenation("IndicateWork/OCUT_MET_dualdesc_", eRec.name);
    if IsExistingFile(FileIndicate)=false and IsExistingFile(FileSaveMET)=false then
      TreatPolyhedralCase(eRecOri.METocut_info, FileSaveMET, eScratchMET);
    fi;
    #
    eScratchCUT:=Concatenation(DirScratch, "OCUT_CUT_info_", eRec.name, "/");
    Exec("mkdir -p ", eScratchCUT);
    FileSaveCUT:=Concatenation("OCUT_CutPolytopes/OCUT_CUT_dualdesc_", eRec.name);
    FileIndicate:=Concatenation("IndicateWork/OCUT_CUT_dualdesc_", eRec.name);
    if IsExistingFile(FileIndicate)=false and IsExistingFile(FileSaveCUT)=false then
      TreatPolyhedralCase(eRecOri.CUTocut_info, FileSaveCUT, eScratchCUT);
    fi;
    #
  fi;
od;

