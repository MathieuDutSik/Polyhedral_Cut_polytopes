RetrieveSingularAutGraphs:=function(n)
  local ListGRA, ListSingular, nbGRA, eGRA, GRAcan, i, eAdj, eRec, GRP, iGRA, ord1, ord2, nRec, eFileB, output;
  ListGRA:=GetIsomorphismTypeGraph(n);
  ListSingular:=[];
  nbGRA:=Length(ListGRA);
  for iGRA in [1..nbGRA]
  do
    Print("iGRA=", iGRA, " / ", nbGRA, "\n");
    eGRA:=ListGRA[iGRA];
    GRAcan:=NullGraph(Group(()), n);
    for i in [1..n]
    do
      for eAdj in eGRA[i]
      do
        AddEdgeOrbit(GRAcan, [i, eAdj]);
      od;
    od;
    if IsConnectedGraph(GRAcan) then
      eRec:=CMC_GetCutPolytope(GRAcan);
      GRP:=LinPolytope_Automorphism(eRec.ListVect);
      if GRP<>eRec.GRPcanonic then
        ord1:=Order(GRP);
        ord2:=Order(eRec.GRPcanonic);
        nRec:=rec(eGRA:=eGRA, ord1:=ord1, ord2:=ord2);
        Add(ListSingular, nRec);
      fi;
    fi;
  od;
  eFileB:=Concatenation("TEX/SingularGraphs_", String(n), "_pre.tex");
  RemoveFileIfExist(eFileB);
  output:=OutputTextFile(eFileB, true);
  AppendTo(output, "Number of vertices \$n=", n, "\$.\n\n");
  AppendTo(output, "Number of singular graphs \$", Length(ListSingular), "\$.\n\n");
  AppendTo(output, "\\begin{itemize}\n");
  for nRec in ListSingular
  do
    AppendTo(output, "\\item Singular graph with \$|Aut(CUTP(G))| = ", nRec.ord1, "\$ and \$|ARes(G)| = ", nRec.ord2, "\$: \n");
    AppendTo(output, "\\begin{enumerate}\n");
    for i in [1..n]
    do
      AppendTo(output, "\\item vertex ", i, " adjacent to");
      for eAdj in nRec.eGRA[i]
      do
        AppendTo(output, " ", eAdj);
      od;
      AppendTo(output, "\n");
    od;
    AppendTo(output, "\\end{enumerate}\n");
  od;
  AppendTo(output, "\\end{itemize}\n");
  LATEX_Compilation(eFileB);
end;
