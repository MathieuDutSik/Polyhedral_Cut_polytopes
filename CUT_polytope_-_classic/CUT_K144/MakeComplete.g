ListSize:=[1,4,4];
GRA:=GRAPH_GetMultiComplement(ListSize);

eRec:=CMC_GetCutPolytope(GRA);
fREC:=eRec.GetCUT_info();
gREC:=eRec.GetMET_info();
