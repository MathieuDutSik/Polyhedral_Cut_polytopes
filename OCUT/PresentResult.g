for n in [3..6]
do
  Print("n=", n, "\n");
  eRec:=FCT__CutPolytope(n);
  eFileSave:=Concatenation("DATA/ListOrb", String(n));
  if IsExistingFile(eFileSave)=true then
    GRP:=LinPolytope_Automorphism(eRec.EXT);
    ListOrb:=ReadAsFunction(eFileSave)();
    eFileOut:=Concatenation("OUT/ListOrb", String(n), ".txt");
    RemoveFileIfExist(eFileOut);
    output:=OutputTextFile(eFileOut, true);
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
    nbOrb:=Length(ListOrb);
    AppendTo(output, "nbOrb=", nbOrb, "\n");
    for iOrb in [1..nbOrb]
    do
      Print("  iOrb=", iOrb, "/", nbOrb, "\n");
      eOrb:=ListOrb[iOrb];
      eStab:=Stabilizer(GRP, eOrb, OnSets);
      AppendTo(output, "iOrb=", iOrb, " |inc|=", Length(eOrb), " |stab|=", Order(eStab), "\n");
      eVect:=GetFacetDescription(eRec, eOrb);
      WriteVector(output, eVect);
    od;
    CloseStream(output);
  fi;
od;
