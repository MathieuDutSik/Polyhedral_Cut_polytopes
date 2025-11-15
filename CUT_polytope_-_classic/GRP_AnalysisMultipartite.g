FileSave:="DATA_GRP/ListEntries";
if IsExistingFilePlusTouch(FileSave) then
  ListCase:=ReadAsFunction(FileSave)();
else
  ListCase:=[];
fi;






FuncInsertQuantity:=function(ListSize)
  local GRA, TheRec, TheCUT, GRP, TheAG, TheQuot, eCase;
  GRA:=GRAPH_GetMultiComplement(ListSize);
  TheRec:=CMC_GetCutPolytope(GRA);
  TheCUT:=TheRec.GetCUT_info();
  GRP:=LinPolytope_Automorphism(TheCUT.EXT);
  TheAG:=Order(GRP) * 2^(1- GRA.order);
  TheQuot:=TheAG / Order(AutGroupGraph(GRA));
  eCase:=rec(ListSize:=ListSize, TheQuot:=TheQuot);
  Add(ListCase, eCase);
  Print("Now, |ListCase|=", Length(ListCase), "\n");
  SaveDataToFile(FileSave, ListCase);
end;

IsPresent:=function(ListSize)
  local eEnt;
  for eEnt in ListCase
  do
    if eEnt.ListSize=ListSize then
      return true;
    fi;
  od;
  return false;
end;


TheChoice:=1;
if TheChoice=1 then
  mMax:=9;
  nbIter:=10000;
  for iter in [1..nbIter]
  do
    m:=Random([3..mMax]);
    eSize:=Random([2..m]);
    eVect:=ListWithIdenticalEntries(m,0);
    for i in [1..m]
    do
      eVect[i] := Random([1..eSize]);
    od;
    PreListSize:=List(Collected(eVect), x->x[2]);
    ePerm:=SortingPerm(PreListSize);
    ListSize:=Permuted(PreListSize, ePerm);
    Print("iter=", iter, " eVect=", eVect, " ListSize=", ListSize, "\n");
    if Length(ListSize) > 1 then
      if IsPresent(ListSize)=false then
        FuncInsertQuantity(ListSize);
      fi;
    fi;
  od;
fi;

ListCaseSearch:=Filtered(ListCase, x->x.TheQuot > 1);

nbCaseSearch:=Length(ListCaseSearch);
for iCaseSearch in [1..nbCaseSearch]
do
  eCaseSearch:=ListCaseSearch[iCaseSearch];
  Print("iCaseSearch=", iCaseSearch, " ListSize=", eCaseSearch.ListSize, " TheQuot=", eCaseSearch.TheQuot, "\n");
od;
