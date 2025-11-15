for m in [4..11]
do
  GRA:=GRAPH_Cycle(m);
  TheRec:=CMC_GetCutPolytope(GRA);
  TheCUT:=TheRec.GetCUT_info();
  GRP:=LinPolytope_Automorphism(TheCUT.EXT);
  TheAG:=Order(GRP) * 2^(1- GRA.order);
  TheQuot:=TheAG / Factorial(m);
  Print("m=", m, " A(G)=", TheAG, " TheQuot=", TheQuot, "\n");
od;
