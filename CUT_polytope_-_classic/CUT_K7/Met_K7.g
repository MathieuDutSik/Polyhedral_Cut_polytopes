GRA:=CompleteGraph(Group(()), 7);

eRec:=CMC_GetCutPolytope(GRA);
EXT:=eRec.GetMET_info().EXT;
EXTfilt:=Filtered(EXT, x->x[1]=0);

ListSets:=DualDescriptionSets(EXTfilt);
