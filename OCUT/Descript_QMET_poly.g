for n in [3..5]
do
  eFileSave:=Concatenation("DATA/Poly_QMET_ListOrb", String(n));
  eRec:=Poly_QMET(n);
  GRP:=LinPolytope_Automorphism(eRec.FAC);
  Print("n=", n, " |GRP|=", Order(GRP), "\n");
  if IsExistingFilePlusTouch(eFileSave)=false then
    eRecGeneral:=GetStandardFCT_DualDescriptionStandard(eRec.FAC, GRP);
    WorkingDim:=RankMat(eRec.FAC);
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
    TmpPath:=Concatenation("./tmpPolyQMET", String(n), "/");
    Exec("mkdir -p ", TmpPath);
    BankPath:=Concatenation("./BankPolyQMET", String(n), "/");
    Exec("mkdir -p ", BankPath);
    eRecGeneral.DataPolyhedral.IsRespawn:=IsRespawn;
    eRecGeneral.DataPolyhedral.Saving:=true;
    eRecGeneral.DataPolyhedral.IsBankSave:=IsBankSave;
    eRecGeneral.DataPolyhedral.ThePathSave:=TmpPath;
    BF:=BankRecording(rec(Saving:=true, BankPath:=BankPath), eRecGeneral.FuncStabilizer, eRecGeneral.FuncIsomorphy, eRecGeneral.FuncInvariant, OnSetsGroupFormalism());
    ListOrb:=__ListFacetByAdjacencyDecompositionMethod(eRec.FAC, GRP, eRecGeneral.DataPolyhedral, BF);
    SaveDataToFilePlusTouch(eFileSave, ListOrb);
  else
    ListOrb:=ReadAsFunction(eFileSave)();
  fi;
  NewListOrb:=SortingOrbit(ListOrb, GRP);
  eFileOut:=Concatenation("OUT/Poly_QMET_ListOrbRay", String(n), ".txt");
  RemoveFileIfExist(eFileOut);
  output:=OutputTextFile(eFileOut, true);
  AppendTo(output, "Vertices of the polytope QMETP", n, "\n");
  nbOrb:=Length(NewListOrb);
  AppendTo(output, "|GRP|=", Order(GRP), " |FAC|=", Length(eRec.FAC), "\n");
  AppendTo(output, "nbOrb=", nbOrb, "\n");
  TotalNumberVertices:=0;
  for iOrb in [1..nbOrb]
  do
    Print("  iOrb=", iOrb, "/", nbOrb, "\n");
    eOrb:=NewListOrb[iOrb];
    eStab:=Stabilizer(GRP, eOrb, OnSets);
    OrbSize:=Order(GRP)/Order(eStab);
    TotalNumberVertices:=TotalNumberVertices + OrbSize;
    AppendTo(output, "iOrb=", iOrb, " |inc|=", Length(eOrb), " |stab|=", Order(eStab), " |O|=", OrbSize, "\n");
    eRay:=__FindFacetInequality(eRec.FAC, eOrb);
    eRayR:=eRay/eRay[1];
    O:=Orbit(eRec.GRPmatr, eRayR, OnPoints);
    TheSelection:=SelectionMinimalQMETelement(O);
    if Length(O)<>OrbSize then
      Print("Inconsistency between matrix group and\n");
      Print("Permutation group: Please panic\n");
      Print(NullMat(5));
    fi;
    eRayRed:=TheSelection{[2..Length(eRec.ListC)+1]};
    PrintQMETelement(output, eRec, eRayRed);
  od;
  AppendTo(output, "TotalNumberVertices=", TotalNumberVertices, "\n");
  AppendTo(output, "nbOrbit Vertices=", nbOrb, "\n");
  AppendTo(output, "|eRec.FAC|=", Length(eRec.FAC), "\n");
  AppendTo(output, "nbOrbit Facet=", Length(Orbits(GRP, [1..Length(eRec.FAC)], OnPoints)), "\n");
  AppendTo(output, "dim=", Length(eRec.FAC[1])-1, "\n");
  CloseStream(output);
od;
