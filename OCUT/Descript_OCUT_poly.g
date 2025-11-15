for n in [3..6]
do
  eFileSave:=Concatenation("DATA/Poly_OCUT_ListOrb", String(n));
  eRec:=Poly_OCUT(n);
  GRP:=LinPolytope_Automorphism(eRec.EXT);
  Print("|EXT|=", Length(eRec.EXT), " |GRP|=", Order(GRP), "\n");
  if IsExistingFilePlusTouch(eFileSave)=false then
    eRecGeneral:=GetStandardFCT_DualDescriptionStandard(eRec.EXT, GRP);
    WorkingDim:=RankMat(eRec.EXT);
    IsRespawn:=function(OrdGRP, EXT, TheDepth)
      if Length(EXT)<WorkingDim+50 then
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
    TmpPath:=Concatenation("./tmpPolyOCUT", String(n), "/");
    Exec("mkdir -p ", TmpPath);
    BankPath:=Concatenation("./BankPolyOCUT", String(n), "/");
    Exec("mkdir -p ", BankPath);
    eRecGeneral.DataPolyhedral.IsRespawn:=IsRespawn;
    eRecGeneral.DataPolyhedral.Saving:=true;
    eRecGeneral.DataPolyhedral.IsBankSave:=IsBankSave;
    eRecGeneral.DataPolyhedral.ThePathSave:=TmpPath;
    eRecGeneral.DataPolyhedral.DualDescriptionFunction:=__DualDescriptionCDD_Reduction;
    BF:=BankRecording(rec(Saving:=true, BankPath:=BankPath), eRecGeneral.FuncStabilizer, eRecGeneral.FuncIsomorphy, eRecGeneral.FuncInvariant, OnSetsGroupFormalism());
    ListOrb:=__ListFacetByAdjacencyDecompositionMethod(eRec.EXT, GRP, eRecGeneral.DataPolyhedral, BF);
    SaveDataToFilePlusTouch(eFileSave, ListOrb);
  else
    ListOrb:=ReadAsFunction(eFileSave)();
  fi;
  #
  NewListOrb:=SortingOrbit(ListOrb, GRP);
  eFileOut:=Concatenation("OUT/Poly_OCUT_ListOrb", String(n), ".txt");
  RemoveFileIfExist(eFileOut);
  output:=OutputTextFile(eFileOut, true);
  AppendTo(output, "Facets of the polytope OCUTP", n, "\n");
  AppendTo(output, "c");
  for i in [1..n-1]
  do
    for j in [i+1..n]
    do
      AppendTo(output, " d", i, j);
    od;
  od;
  for i in [2..n]
  do
    AppendTo(output, " w", i);
  od;
  AppendTo(output, "\n");
  nbOrb:=Length(NewListOrb);
  AppendTo(output, "nbOrb=", nbOrb, "\n");
  TotalNumberFacets:=0;
  for iOrb in [1..nbOrb]
  do
    Print("  iOrb=", iOrb, "/", nbOrb, "\n");
    eOrb:=NewListOrb[iOrb];
    eStab:=Stabilizer(GRP, eOrb, OnSets);
    OrbSize:=Order(GRP)/Order(eStab);
    TotalNumberFacets:=TotalNumberFacets + OrbSize;
    AppendTo(output, "iOrb=", iOrb, " |inc|=", Length(eOrb), " |stab|=", Order(eStab), " |O|=", OrbSize, "\n");
    eFac:=__FindFacetInequality(eRec.EXTproj, eOrb);
    O:=Orbit(eRec.GRPmatrResCongr, eFac, OnPoints);
    if Length(O)<>OrbSize then
      Print("Big error in group orbit. Please panic\n");
      Print(NullMat(5));
    fi;
    TheSelection:=SelectionMinimalOCUTfacet(eRec, O);
    WriteVector(output, TheSelection);
  od;
  AppendTo(output, "TotalNumberFacets=", TotalNumberFacets, "\n");
  AppendTo(output, "nbOrb Facet=", nbOrb, "\n");
  AppendTo(output, "|EXT|=", Length(eRec.EXT), "\n");
  AppendTo(output, "nbOrb Vertices=", Length(Orbits(GRP, [1..Length(eRec.EXT)], OnPoints)), "\n");
  AppendTo(output, "dim=", RankMat(eRec.EXT)-1, "\n");
  CloseStream(output);
od;
