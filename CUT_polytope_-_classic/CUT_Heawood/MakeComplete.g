ListEdges:=[
[0,1],[0,5],[0,13],
[1,10],[1,2],
[2,7],[2,3],
[3,12],[3,4],
[4,9],[4,5],
[5,6],
[6,11],[6,7],
[7,8],
[8,13],[8,9],
[9,10],
[10,11],
[11,12],
[12,13]];
GRA:=NullGraph(Group(()),14);
for eEdge in ListEdges
do
  i:=eEdge[1]+1;
  j:=eEdge[2]+1;
  AddEdgeOrbit(GRA, [i,j]);
  AddEdgeOrbit(GRA, [j,i]);
od;

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
SaveDataToFile("DataHeawood", fRec);
