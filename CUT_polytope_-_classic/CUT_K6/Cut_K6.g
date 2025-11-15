GRA:=CompleteGraph(Group(()), 6);

eRec:=CMC_GetCutPolytope(GRA);
#EXT:=eRec.GetCUT_info().EXT;
EXT:=eRec.GetMET_info().EXT;
ListSet:=DualDescriptionSets(EXT);
