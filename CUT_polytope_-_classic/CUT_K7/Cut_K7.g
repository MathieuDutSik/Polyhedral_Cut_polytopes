GRA:=CompleteGraph(Group(()), 7);

eRec:=CMC_GetCutPolytope(GRA);
EXT:=eRec.GetCUT_info().EXT;
dim:=Length(EXT[1]);
EXTred:=List(EXT, x->x{[2..dim]});
EXTfilt:=Filtered(EXTred, x->x*x > 0);

ListSets:=DualDescriptionSets(EXTfilt);
