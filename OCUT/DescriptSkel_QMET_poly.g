for n in [3..9]
do
  eFileSave:=Concatenation("DATA/Poly_QMET_skel", String(n));
  eRec:=Poly_QMET(n);
  GRP:=LinPolytope_Automorphism(eRec.FAC);
  Print("n=", n, " |GRP|=", Order(GRP), "\n");
  if IsExistingFilePlusTouch(eFileSave)=false then
    BoundingSet:=[];
    TheSkel:=SkeletonGraph(GRP, eRec.FAC, BoundingSet);
    SaveDataToFilePlusTouch(eFileSave, TheSkel);
  else
    TheSkel:=ReadAsFunction(eFileSave)();
  fi;
  TheDiam:=Diameter(TheSkel);
  eFileOut:=Concatenation("OUT/Poly_QMET_skel", String(n), ".txt");
  RemoveFileIfExist(eFileOut);
  output:=OutputTextFile(eFileOut, true);
  AppendTo(output, "Ridge graph of the polytope QMETP", n, "\n");
  nbFac:=Length(eRec.FAC);
  AppendTo(output, "|GRP|=", Order(GRP), " |FAC|=", nbFac, "\n");
  AppendTo(output, "diameter=", TheDiam, "\n");
  nbRep:=Length(eRec.ListRepresentatives);
  AppendTo(output, "nbRep facet=", nbRep, "\n");
  for iRep in [1..nbRep]
  do
    eRepr:=eRec.ListRepresentatives[iRep];
    eStab:=Stabilizer(GRP, eRepr, OnPoints);
    LAdj:=Adjacency(TheSkel, eRepr);
    eDiff:=Difference([1..nbFac], [eRepr]);
    LNonAdj:=Difference(eDiff, Set(LAdj));
    AppendTo(output, "iRep=", iRep, "  Representative ", eRec.LSymbIneq[eRepr], " |nonadj|=", Length(LNonAdj), " |stab|=", Order(eStab), "\n");
    for eNAdj in LNonAdj
    do
      AppendTo(output, "  ", eRec.LSymbIneq[eNAdj], "\n");
    od;
  od;
  CloseStream(output);
od;
