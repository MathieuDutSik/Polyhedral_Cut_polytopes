PL:=Prism(3);
GRA1:=PlanGraphToGRAPE(PL);
GRA:=GRAPH_Pyramid(GRA1);

eRec:=CMC_GetCutPolytope(GRA);
EXT:=eRec.ListVect;
FileSaveGRP:="TheData_GRP";
if IsExistingFile(FileSaveGRP) then
  GRP:=ReadAsFunction(FileSaveGRP)();
else
  GRP:=LinPolytope_Automorphism(EXT);
  SaveDataToFile(FileSaveGRP, GRP);
fi;

eRec:=CMC_GetCutPolytope(GRA);
EXT:=eRec.ListVect;

LOrb:=ReadAsFunction("ListOrbitEXT")();

fRec:=rec(EXT:=EXT, GRP:=GRP, LOrb:=LOrb, eRec:=eRec);
SaveDataToFile("DataPyrPrism3", fRec);
