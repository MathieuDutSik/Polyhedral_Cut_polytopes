GRA:=CompleteGraph(Group(()), 7);

eRec:=CMC_GetCutPolytope(GRA);
fRec:=eRec.GetCUT_info();
gRec:=eRec.GetMET_info();

SaveDataToFile("Example_07_CUT7", fRec);
SaveDataToFile("Example_08_MET7", gRec);
