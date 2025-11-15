GRA_K2:=CompleteGraph(Group(()), 4);
GRA:=ComplementGraph(GRAPH_ZeroExtension(GRA_K2, 3));

eRec:=CMC_GetCutPolytope(GRA);
EXT:=eRec.ListVect;
GRP:=LinPolytope_Automorphism(EXT);
LOrb:=ReadAsFunction("ListOrbitEXT")();

fRec:=rec(EXT:=EXT, GRP:=GRP, LOrb:=LOrb, eRec:=eRec);
SaveDataToFile("DataK7-K4", fRec);
