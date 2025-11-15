ListSize:=[1,1,1,2,2];
GRA:=GRAPH_GetMultiComplement(ListSize);

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
SaveDataToFile("DataK11122", fRec);
