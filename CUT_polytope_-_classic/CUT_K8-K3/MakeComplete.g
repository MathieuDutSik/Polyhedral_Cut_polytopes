GRA_K2:=CompleteGraph(Group(()), 3);
GRA:=ComplementGraph(GRAPH_ZeroExtension(GRA_K2, 5));

eRec:=CMC_GetCutPolytope(GRA);
EXT:=eRec.ListVect;
GRP:=LinPolytope_Automorphism(EXT);
LOrb:=ReadAsFunction("ListOrbitEXT")();

fRec:=rec(EXT:=EXT, GRP:=GRP, LOrb:=LOrb, eRec:=eRec);
SaveDataToFile("DataK8-K3", fRec);
