n:=6;
GRA:=CompleteGraph(Group(()), n);
eRecOri:=CMC_GetCutPolytopeOriented(GRA);



PrintPolytopeInfo:=function(string, EXT)
  local rnk, nbVert, GRP;
  rnk:=RankMat(EXT);
  nbVert:=Length(EXT);
  GRP:=LinPolytope_Automorphism(EXT);
  Print("INFO : ", string, " rnk/siz=", rnk, " / ", nbVert, " |G|=", Order(GRP), "\n");
end;


eRec1:=Poly_OCUT(n);
eRec2:=Poly_QMET(n);


PrintPolytopeInfo("EXTocut", eRecOri.CUTocut_info.EXT);
PrintPolytopeInfo("FACocut", eRecOri.METocut_info.EXT);


PrintPolytopeInfo("EXTocut (native)", eRec1.EXTproj);
PrintPolytopeInfo("METocut (native)", eRec2.FAChom);

test1:=LinPolytope_Isomorphism(eRec1.EXTproj, eRecOri.CUTocut_info.EXT);
Print("test1=", test1, "\n");
test2:=LinPolytope_Isomorphism(eRec2.FAChom, eRecOri.METocut_info.EXT);
Print("test2=", test2, "\n");

