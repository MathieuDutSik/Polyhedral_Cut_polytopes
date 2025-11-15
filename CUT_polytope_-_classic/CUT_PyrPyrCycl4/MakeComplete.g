GRA1:=GRAPH_Cycle(4);
GRA2:=GRAPH_Pyramid(GRA1);
GRA:=GRAPH_Pyramid(GRA2);

eRec:=CMC_GetCutPolytope(GRA);
EXT:=eRec.GetCUT_info().EXT;
FileSaveGRP:="TheData_GRP";
if IsExistingFile(FileSaveGRP) then
  GRP:=ReadAsFunction(FileSaveGRP)();
else
  GRP:=LinPolytope_Automorphism(EXT);
  SaveDataToFile(FileSaveGRP, GRP);
fi;

LOrb:=ReadAsFunction("ListOrbitEXT")();

fRec:=rec(EXT:=EXT, GRP:=GRP, LOrb:=LOrb, eRec:=rec(GRA:=eRec.GRA, ListVect:=EXT, ListSet:=eRec.GetListSet()));
SaveDataToFile("DataPyrPyrCycl4", fRec);
