for m in [4..6]
do
  GRA:=GRAPH_Path(m);
  TheRec:=CMC_GetCutPolytope(GRA);
  TheCUT:=TheRec.GetCUT_info();
  GRP:=LinPolytope_Automorphism(TheCUT.EXT);
  TheAG:=Order(GRP) * 2^(1- GRA.order);
  TheQuot:=TheAG / Factorial(m-1);
  Print("m=", m, " A(G)=", TheAG, " TheQuot=", TheQuot, "\n");
od;
