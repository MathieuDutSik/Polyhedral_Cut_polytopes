GetHemimetricPolytope:=function(n, k)
  local ListSets, nbSet, FAC, eSet, i, eFAC, j, fSet, pos, eSign;
  ListSets:=Combinations([1..n],k);
  nbSet:=Length(ListSets);
  FAC:=[];
  for eSet in Combinations([1..n], k+1)
  do
    for i in [1..k+1]
    do
      eFAC:=ListWithIdenticalEntries(nbSet+1,0);
      for j in [1..k+1]
      do
        fSet:=Difference(eSet, [eSet[j]]);
        pos:=Position(ListSets, fSet);
        if i=j then
          eSign:=-1;
        else
          eSign:=1;
        fi;
        eFAC[pos+1]:=eSign;
      od;
      Add(FAC, eFAC);
    od;
    eFAC:=ListWithIdenticalEntries(nbSet+1,0);
    eFAC[1]:=k+1;
    for j in [1..k+1]
    do
      fSet:=Difference(eSet, [eSet[j]]);
      pos:=Position(ListSets, fSet);
      eFAC[pos+1]:=-1;
    od;
    Add(FAC, eFAC);
  od;
  return rec(ListSets:=ListSets, 
             FAC:=FAC);
end;
