#
# K_{1,4,4}
#

ListSize:=[1,4,4];
GRA:=GRAPH_GetMultiComplement(ListSize);

eRec:=CMC_GetCutPolytope(GRA);
EXT:=eRec.ListVect;
GRP:=LinPolytope_Automorphism(EXT);

SYMPOL_PrintMatrix("SYMPOL_K144.ext", EXT);
SYMPOL_PrintGroup("SYMPOL_K144.grp", Length(EXT), GRP);

#
# K_{8}
#

GRA:=ComplementGraph(NullGraph(Group(()),8));

eRec:=CMC_GetCutPolytope(GRA);
EXT:=eRec.ListVect;
GRP:=LinPolytope_Automorphism(EXT);

SYMPOL_PrintMatrix("SYMPOL_K8.ext", EXT);
SYMPOL_PrintGroup("SYMPOL_K8.grp", Length(EXT), GRP);

