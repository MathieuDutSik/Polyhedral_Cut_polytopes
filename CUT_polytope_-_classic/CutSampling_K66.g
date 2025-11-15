eListSize:=[6,6];

GRA:=GRAPH_GetMultiComplement(eListSize);
eRec:=CMC_GetCutPolytope(GRA);
eRecCUT:=eRec.GetCUT_info();
EXT:=eRecCUT.EXT;
#
TheLevel:=30;
LOrb:=SamplingStandardLevel(EXT, TheLevel);
ListLen:=List(LOrb, Length);
Print("ListLen=", ListLen, "\n");
