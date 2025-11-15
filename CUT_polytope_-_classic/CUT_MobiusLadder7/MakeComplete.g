GRA:=MobiusLadderGraph(7);

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
SaveDataToFile("DataMobiusLadder7", fRec);
