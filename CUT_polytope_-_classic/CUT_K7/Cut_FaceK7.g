GRA:=CompleteGraph(Group(()), 7);

eRec:=CMC_GetCutPolytope(GRA);
EXT:=eRec.GetCUT_info().EXT;
FAC:=eRec.GetMET_info().EXT;

eFAC:=FAC[1];
EXTfilt:=Filtered(EXT, x->x*eFAC=0);
Print("|EXTfilt|=", Length(EXTfilt), "\n");



ListSets:=DualDescriptionSets(EXTfilt);
Print("|ListSets|=", Length(ListSets), "\n");
