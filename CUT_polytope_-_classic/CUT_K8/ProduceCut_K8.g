GRA:=CompleteGraph(Group(()), 8);


eRec:=CMC_GetCutPolytope(GRA);
eRecCUT:=eRec.GetCUT_info();


EXT:=eRecCUT.EXT;
GRP:=eRecCUT.GRP;

FileEXT:="CutK8.ext";
SYMPOL_PrintMatrix(FileEXT, EXT);

n:=Length(EXT);
FileGRP:="CutK8.grp";
SYMPOL_PrintGroup(FileGRP, n, GRP);
