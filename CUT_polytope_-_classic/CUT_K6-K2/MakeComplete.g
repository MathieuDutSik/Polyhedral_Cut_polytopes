GRA_K2:=CompleteGraph(Group(()), 2);
GRA:=ComplementGraph(GRAPH_ZeroExtension(GRA_K2, 4));

eRec:=CMC_GetCutPolytope(GRA);
EXT:=eRec.ListVect;
GRP:=LinPolytope_Automorphism(EXT);
LOrb:=ReadAsFunction("ListOrbitEXT")();

fRec:=rec(EXT:=EXT, GRP:=GRP, LOrb:=LOrb, eRec:=eRec);
SaveDataToFile("DataK6-K2", fRec);
