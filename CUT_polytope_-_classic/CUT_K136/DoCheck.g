ListSize:=[1,3,6];
GRA:=GRAPH_GetMultiComplement(ListSize);

eRec:=CMC_GetCutPolytope(GRA);
EXT:=eRec.ListVect;
GRP:=LinPolytope_Automorphism(EXT);
LOrb:=ReadAsFunction("ListOrbitEXT_buggy")();

test:=RepresentativeAction(GRP, LOrb[2], LOrb[8], OnSets);

LOrbRed:=[];
FuncInsert:=function(eOrb)
  local fOrb;
  for fOrb in LOrbRed
  do
    if RepresentativeAction(GRP, eOrb, fOrb, OnSets)<>fail then
      return;
    fi;
  od;
  Add(LOrbRed, eOrb);
end;

for eOrb in LOrb
do
  FuncInsert(eOrb);
od;

SaveDataToFile("ListOrbitEXT_red", LOrbRed);

