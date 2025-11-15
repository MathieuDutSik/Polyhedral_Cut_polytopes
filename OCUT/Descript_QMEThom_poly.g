for n in [3..6]
do
  eFileSave:=Concatenation("DATA/Poly_QMEThom_ListOrb", String(n));
  eRec:=Poly_QMET(n);
  GRP:=LinPolytope_Automorphism(eRec.FAChom);
  Print("n=", n, " |GRP|=", Order(GRP), "\n");
  if IsExistingFilePlusTouch(eFileSave)=false then
    eRecGeneral:=GetStandardFCT_DualDescriptionStandard(eRec.FAChom, GRP);
    WorkingDim:=RankMat(eRec.FAChom);
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
    TmpPath:=Concatenation("./tmpPolyQMEThom", String(n), "/");
    Exec("mkdir -p ", TmpPath);
    BankPath:=Concatenation("./BankPolyQMEThom", String(n), "/");
    Exec("mkdir -p ", BankPath);
    eRecGeneral.DataPolyhedral.IsRespawn:=IsRespawn;
    eRecGeneral.DataPolyhedral.Saving:=true;
    eRecGeneral.DataPolyhedral.IsBankSave:=IsBankSave;
    eRecGeneral.DataPolyhedral.ThePathSave:=TmpPath;
    BF:=BankRecording(rec(Saving:=true, BankPath:=BankPath), eRecGeneral.FuncStabilizer, eRecGeneral.FuncIsomorphy, eRecGeneral.FuncInvariant, OnSetsGroupFormalism());
    ListOrb:=__ListFacetByAdjacencyDecompositionMethod(eRec.FAChom, GRP, eRecGeneral.DataPolyhedral, BF);
    SaveDataToFilePlusTouch(eFileSave, ListOrb);
  else
    ListOrb:=ReadAsFunction(eFileSave)();
  fi;
  NewListOrb:=SortingOrbit(ListOrb, GRP);
  nbFAChom:=Length(eRec.FAChom);
  eSetNew:=[];
  FACzero:=[];
  for iIdx in [1..nbFAChom]
  do
    if eRec.ListNatureIneq[eRec.ListIndexes[iIdx]]=0 then
      Add(FACzero, eRec.FAC[eRec.ListIndexes[iIdx]]);
      Add(eSetNew, iIdx);
    fi;
  od;
  eStabSymRev:=Stabilizer(GRP, eSetNew, OnSets);
  eFileOutPoly:=Concatenation("OUT/Poly_QMEThom_ListOrbRay", String(n), ".txt");
  RemoveFileIfExist(eFileOutPoly);
  eFileOutCone:=Concatenation("OUT/Poly_QMEThom_ConeInfo", String(n), ".txt");
  RemoveFileIfExist(eFileOutCone);
  outputPoly:=OutputTextFile(eFileOutPoly, true);
#  outputCone:=OutputTextFile(eFileOutCone, true);
  AppendTo(outputPoly, "Vertices of the polytope WQMETP", n, "\n");
  AppendTo(outputPoly, "|GRPsymRevSwitch|=", Order(GRP), " |FAC|=", nbFAChom, "\n");
  AppendTo(outputPoly, "|GRPsymRev|=", Order(eStabSymRev), "\n");
  nbOrb:=Length(NewListOrb);
  AppendTo(outputPoly, "nbOrb=", nbOrb, "\n");
  TotalNumberVertices:=0;
  for iOrb in [1..nbOrb]
  do
    Print("  iOrb=", iOrb, "/", nbOrb, "\n");
    eOrb:=NewListOrb[iOrb];
    eStab:=Stabilizer(GRP, eOrb, OnSets);
    OrbSiz:=Order(GRP)/Order(eStab);
    TotalNumberVertices:=TotalNumberVertices + OrbSiz;
    AppendTo(outputPoly, "iOrb=", iOrb, "/", nbOrb, " |inc|=", Length(eOrb), " |stab|=", Order(eStab), " |O|=", OrbSiz, "\n");
    eRay:=__FindFacetInequality(eRec.FAChom, eOrb)*eRec.NSP;
    eRayR:=eRay/eRay[1];
    O:=Orbit(eRec.GRPmatr, eRayR, OnPoints);
    TheSelection:=SelectionMinimalQMETelement(O);
    if Length(O)<>OrbSiz then
      Print("Inconsistency between matrix group and\n");
      Print("Permutation group: Please panic\n");
      Print(NullMat(5));
    fi;
    eRayRed:=TheSelection{[2..Length(eRec.ListC)+1]};
    PrintQMETelement(outputPoly, eRec, eRayRed);
    ListOrbSplit:=__IndividualLifting(eOrb, eStabSymRev, GRP);
    nbOrbSplit:=Length(ListOrbSplit);
    AppendTo(outputPoly, " under GRPsymRev, the orbit splits into ", nbOrbSplit, " orbits\n");
    for iOrbSplit in [1..nbOrbSplit]
    do
      AppendTo(outputPoly, "iOrbSplit=", iOrbSplit, "/", nbOrbSplit, "\n");
      eRay1:=__FindFacetInequality(eRec.FAChom, ListOrbSplit[iOrbSplit]);
      eRay:=eRay1*eRec.NSP;
      ListIncdTot:=Filtered(eRec.FAChom, x->x*eRay1=0);
      ListIncd:=Filtered(eRec.FAChom{eSetNew}, x->x*eRay1=0);
      rnkTot:=PersoRankMat(ListIncdTot);
      rnk:=PersoRankMat(ListIncd);
      eRayR:=eRay/eRay[1];
      eRayRed:=eRayR{[2..Length(eRec.ListC)+1]};
      Print("rnk=", rnk, " rnkTot=", rnkTot, " eRay1=", eRay1, " eRayRed=", eRayRed, "\n");
      PrintQMETelement(outputPoly, eRec, eRayRed);
    od;
    AppendTo(outputPoly, "\n");
#    Print(NullMat(5));
  od;
  AppendTo(outputPoly, "TotalNumberVertices=", TotalNumberVertices, "\n");
  AppendTo(outputPoly, "nbOrbitVertices=", nbOrb, "\n");
  AppendTo(outputPoly, "|eRec.FAChom|=", Length(eRec.FAChom), "\n");
  AppendTo(outputPoly, "nbOrb Facet=", Length(Orbits(GRP, [1..Length(eRec.FAChom)], OnPoints)), "\n");
  AppendTo(outputPoly, "dim=", Length(eRec.FAChom[1])-1, "\n");
  CloseStream(outputPoly);
#  CloseStream(outputCone);
od;
