n:=5;



GetOcut:=function(eSet)
  local eMat, i, j;
  eMat:=NullMat(n,n);
  for i in [1..n]
  do
    for j in [1..n]
    do
      if Position(eSet, i)<>fail and Position(eSet, j)=fail then
        eMat[i][j]:=1;
      fi;
    od;
  od;
  return eMat;
end;

DoOrientedSwitching:=function(eMat, eSet)
  local RetMat, i, j, eInt;
  RetMat:=NullMat(n,n);
  for i in [1..n]
  do
    for j in [1..n]
    do
      if i<>j then
        eInt:=Intersection(eSet, Set([i,j]));
        if Length(eInt)=1 then
          RetMat[i][j]:=1 - eMat[j][i];
        else
          RetMat[i][j]:=eMat[i][j];
        fi;
      fi;
    od;
  od;
  return RetMat;
end;


ListMat:=[];
ListPair:=[];

for eSet in Combinations([1..n])
do
  eMat:=GetOcut(eSet);
  for fSet in Combinations([1..n-1])
  do
    fMat:=DoOrientedSwitching(eMat, fSet);
    pos:=Position(ListMat, fMat);
    ePair:=rec(eSet:=eSet, fSet:=fSet);
    if pos=fail then
      Add(ListMat, fMat);
      Add(ListPair, [ePair]);
    else
      Add(ListPair[pos], ePair);
    fi;
  od;
od;

