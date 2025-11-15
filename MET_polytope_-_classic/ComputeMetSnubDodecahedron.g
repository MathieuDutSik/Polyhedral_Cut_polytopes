PreListGraphName:=["SnubDodecahedron"];


PreListPL:=[];
PreListGRA:=[];
ListDisc:=[];
for eName in PreListGraphName
do
  ePL:=ArchimedeanPolyhedra(eName);
  eGRA:=PlanGraphToGRAPE(ePL);
  Add(PreListPL, ePL);
  Add(PreListGRA, eGRA);
  Add(ListDisc, [Length(ePL)]);
od;

ePerm:=SortingPerm(ListDisc);

ListGraphName:=Permuted(PreListGraphName, ePerm);
ListPL:=Permuted(PreListPL, ePerm);
ListGRA:=Permuted(PreListGRA, ePerm);
ListDiscSort:=Permuted(ListDisc, ePerm);
nbCase:=Length(ListGRA);

for iCase in [1..nbCase]
do
  eName:=ListGraphName[iCase];
  eDisc:=ListDiscSort[iCase];
  Print("iCase=", iCase, " / ", nbCase, " eName=", eName, " eDisc=", eDisc, "\n");
od;


ListOrbitInfo:=[];
for iCase in [1..nbCase]
do
  eName:=ListGraphName[iCase];
  Print("eName=", eName, "\n");
  ePL:=ListPL[iCase];
  eGRA:=ListGRA[iCase];
  FileSave:=Concatenation("DATA_MET/MET_info_", eName);
  if IsExistingFilePlusTouch(FileSave) then
    eOrbitInfo:=ReadAsFunction(FileSave)();
  else
    eRecFCT:=CMC_GetCutPolytope(eGRA);
    eOrbitInfo:=eRecFCT.GetMET_info_symbolic();
    Add(ListOrbitInfo, eOrbitInfo);
    SaveDataToFilePlusTouch(FileSave, eOrbitInfo);
  fi;
  Add(ListOrbitInfo, eOrbitInfo);
od;

