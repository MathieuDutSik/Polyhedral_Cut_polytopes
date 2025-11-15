PL:=AntiPrism(4);
GRA:=PlanGraphToGRAPE(PL);

eRec:=CMC_GetCutPolytope(GRA);
EXT:=eRec.ListVect;
FileSaveGRP:="TheData_GRP";
if IsExistingFile(FileSaveGRP) then
  GRP:=ReadAsFunction(FileSaveGRP)();
else
  GRP:=LinPolytope_Automorphism(EXT);
  SaveDataToFile(FileSaveGRP, GRP);
fi;

LOrb:=ReadAsFunction("ListOrbitEXT")();

fRec:=rec(EXT:=EXT, GRP:=GRP, LOrb:=LOrb, eRec:=eRec);
SaveDataToFile("DataAntiPrism4", fRec);

