for n in [3..7]
do
  eFileSave:=Concatenation("DATA/ListOrb", String(n));
  if IsExistingFilePlusTouch(eFileSave)=false then
    eRec:=FCT__CutPolytope(n);
    GRP:=LinPolytope_Automorphism(eRec.EXT);
    eRecGeneral:=GetStandardFCT_DualDescriptionStandard(eRec.EXT, GRP);
    WorkingDim:=RankMat(eRec.EXT);
    IsRespawn:=function(OrdGRP, EXT, TheDepth)
      if Length(EXT)<WorkingDim+7 then
        return false;
      fi;
      if OrdGRP>=100 and TheDepth<=2 then
        return true;
      fi;
      if Length(EXT)>45 and OrdGRP>20 then
        return true;
      fi;
      if Length(EXT)>=50 and OrdGRP>=8 then
        return true;
      fi;
      if Length(EXT) < WorkingDim+15 and OrdGRP < 20 then
        return false;
      fi;
      if OrdGRP<10 then
        return false;
      fi;
      return true;
    end;
    IsBankSave:=function(EllapsedTime, OrdGRP, EXT, TheDepth)
      if EllapsedTime>=600 then
        return true;
      fi;
      if Length(EXT)<=40 then
        return false;
      fi;
      return true;
    end;
    TmpPath:=Concatenation("./tmp", String(n), "/");
    Exec("mkdir -p ", TmpPath);
    BankPath:=Concatenation("./Bank", String(n), "/");
    Exec("mkdir -p ", BankPath);
    eRecGeneral.DataPolyhedral.IsRespawn:=IsRespawn;
    eRecGeneral.DataPolyhedral.Saving:=true;
    eRecGeneral.DataPolyhedral.IsBankSave:=IsBankSave;
    eRecGeneral.DataPolyhedral.ThePathSave:=TmpPath;
    BF:=BankRecording(rec(Saving:=true, BankPath:=BankPath), eRecGeneral.FuncStabilizer, eRecGeneral.FuncIsomorphy, eRecGeneral.FuncInvariant, OnSetsGroupFormalism());
    ListOrb:=__ListFacetByAdjacencyDecompositionMethod(eRec.EXT, GRP, eRecGeneral.DataPolyhedral, BF);
    SaveDataToFilePlusTouch(eFileSave, ListOrb);
  fi;
od;
