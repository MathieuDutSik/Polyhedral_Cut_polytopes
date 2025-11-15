ListSize:=[4,5];
GRA:=GRAPH_GetMultiComplement(ListSize);

eRec:=CMC_GetCutPolytope(GRA);
EXT:=eRec.ListVect;
GRP:=LinPolytope_Automorphism(EXT);


LOrb:=ReadAsFunction("ListOrbitEXT")();

fRec:=rec(EXT:=EXT, GRP:=GRP, LOrb:=LOrb, eRec:=eRec);
SaveDataToFile("DataK45", fRec);
