#
#
#
#
# The group part
OSMET__ListElement:=function(n,m)
  local eSet, ListElement;
  ListElement:=[];
  for eSet in Combinations([1..n], m+1)
  do
    Add(ListElement, eSet);
    Add(ListElement, Permuted(eSet, (1,2) ));
  od;
  return ListElement;
end;

OSMET__Canonicalize:=function(eList)
  local eSet, ePerm;
  eSet:=Set(eList);
  ePerm:=PermList(List(eList, x->Position(eSet,x)));
  if SignPerm(ePerm)=1 then
    return eSet;
  else
    return Permuted(eSet, (1,2));
  fi;
end;


OSMET__InversionTransformation:=function(n,m)
  local ListElement, eList, iPos, eSet, fSet, jPos;
  ListElement:=OSMET__ListElement(n,m);
  eList:=[1];
  for iPos in [1..Length(ListElement)]
  do
    eSet:=ListElement[iPos];
    fSet:=Permuted(eSet, (1,2));
    jPos:=Position(ListElement, fSet);
    eList[iPos+1]:=jPos+1;
  od;
  return PermList(eList);
end;


OSMET__SymmetryGroup:=function(n,m)
  local ListElement, ListGen, eGen, iPos, eSet, fSet, jPos, eList;
  ListElement:=OSMET__ListElement(n,m);
  ListGen:=[];
  for eGen in GeneratorsOfGroup(SymmetricGroup(n))
  do
    eList:=[1];
    for iPos in [1..Length(ListElement)]
    do
      eSet:=ListElement[iPos];
      fSet:=OnTuples(eSet, eGen);
      jPos:=Position(ListElement, OSMET__Canonicalize(fSet));
      eList[iPos+1]:=jPos+1;
    od;
    Add(ListGen, PermList(eList));
  od;
  Add(ListGen, OSMET__InversionTransformation(n,m));
  return Group(ListGen);
end;







#
#
#
# Oriented partition generalities
SPAN:=function(OP, EltAdd)
  local k, size, iSet, TheList, U, kPos, idx;
  k:=Length(OP);
  TheList:=[];
  for iSet in [1..k]
  do
    U:=ShallowCopy(OP);
    U[iSet]:=Union(U[iSet], [EltAdd]);
    Add(TheList, U);
  od;
  for iSet in [0..k]
  do
    U:=ShallowCopy([]);
    idx:=1;
    for kPos in [0..k]
    do
      if kPos=iSet then
        U[kPos+1]:=[EltAdd];
      else
        U[kPos+1]:=OP[idx];
        idx:=idx+1;
      fi;
    od;
    Add(TheList, U);
  od;
  return TheList;
end;

# Increase level of a whole set of oriented partitions
IncreaseLevel:=function(ListO, EltAdd)
  local List, iO;
  List:=[];
  for iO in [1..Length(ListO)]
  do
    List:=Union(List, SPAN(ListO[iO], EltAdd));
  od;
  return List;
end;


# Span the list of oriented partitions on an arbitrary set
SpanOrientedPartGeneral:=function(n)
  local U, iElt, Resu, iRes, GRP, ListGenNew, eGen, eList, jPart, iPart, ePart, fPart;
  U:=[[[1]]];
  for iElt in [2..n]
  do
    U:=IncreaseLevel(U, iElt);
  od;
  Resu:=[];
  for iRes in [1..Length(U)]
  do
    if Length(U[iRes])>1 then
      Add(Resu, U[iRes]);
    fi;
  od;
  GRP:=SymmetricGroup(n);
  ListGenNew:=[];
  for eGen in GeneratorsOfGroup(GRP)
  do
    eList:=[];
    for iPart in [1..Length(Resu)]
    do
      ePart:=Resu[iPart];
      fPart:=OnTuplesSets(ePart, eGen);
      jPart:=Position(Resu, fPart);
      Add(eList, jPart);
    od;
    Add(ListGenNew, PermList(eList));
  od;
  eList:=[];
  for iPart in [1..Length(Resu)]
  do
    ePart:=Resu[iPart];
    fPart:=Reversed(ePart);
    jPart:=Position(Resu, fPart);
    Add(eList, jPart);
  od;
  Add(ListGenNew, PermList(eList));
  return rec(ListOrientedPartition:=Resu, GroupExt:=Group(ListGenNew));
end;




MergeOrientedPart:=function(OP1, OP2)
  local OP, i;
  OP:=ShallowCopy(OP1);
  for i in [1..Length(OP2)]
  do
    OP[Length(OP)+1]:=OP2[i];
  od;
  return(OP);
end;

Reversal:=function(OP)
  local k, U, iSet;
  k:=Length(OP);
  U:=[];
  for iSet in [1..k]
  do
    U[iSet]:=OP[k+1-iSet];
  od;
  return U;
end;

PartOutput:=function(output, OPset)
  local iOP, OP, iP, iS, i, VS;
  for iOP in [1..Length(OPset)]
  do
    OP:=OPset[iOP];
    for iP in [1..Length(OP)]
    do
      iS:=OP[iP];
      for i in [1..Length(iS)]
      do
        AppendTo(output, " ", iS[i]);
      od;
      if iP<Length(OP) then
        AppendTo(output, " |");
      fi;
    od;
  od;
end;


# do all cyclic shifts of an oriented partition
CyclicShifts:=function(OP)
  local k, List, iShift, S, U, iU;
  k:=Length(OP);
  List:=[];
  List[1]:=OP;
  for iShift in [2..k]
  do
    S:=List[iShift-1];
    U:=ShallowCopy([]);
    for iU in [2..k]
    do
      U[iU-1]:=S[iU];
    od;
    U[k]:=S[1];
    List[Length(List)+1]:=U;
  od;
  return List;
end;



# do all cyclic shifts of a set of oriented partitions
CyclicShiftsSet:=function(SETOP)
  local U, iU, List, iL; 
  List:=[];
  for iL in [1..Length(SETOP)]
  do
    U:=CyclicShifts(SETOP[iL]);
    for iU in [1..Length(U)]
    do
      List[Length(List)+1]:=U[iU];
    od;
  od;
  return List;
end;



PartToVector:=function(OP, ListEltM)
  local k, size, Mp1, TheVect, eSet, eList, eVal, iSub, eResList, ePerm;
  k:=Length(OP);
  size:=Sum(List(OP, x->Length(x)));
  Mp1:=Length(ListEltM[1]);
  TheVect:=[0];
  #
  for eSet in ListEltM
  do
    eList:=[];
    for eVal in eSet
    do
      for iSub in [1..Length(OP)]
      do
        if eVal in OP[iSub] then
          Add(eList, iSub);
        fi;
      od;
    od;
    eResList:=Set(eList);
    if Length(eResList)<Mp1 then
      Add(TheVect, 0);
    else
      ePerm:=PermList(List(eList, x->Position(eResList, x)));
      if SignPerm(ePerm)=1 then
        Add(TheVect, 1);
      else
        Add(TheVect, 0);
      fi;
    fi;
  od;
  return TheVect;
end;






OSMET__FuncInequality:=function(ListEltM, eList, iIdx)
  local TheIneq, i, fList;
  TheIneq:=ListWithIdenticalEntries(1+Length(ListEltM), 0);
  TheIneq[Position(ListEltM, OSMET__Canonicalize(eList))+1]:=-1;
  for i in [1..Length(eList)]
  do
    fList:=ShallowCopy(eList);
    fList[i]:=iIdx;
    TheIneq[Position(ListEltM, OSMET__Canonicalize(fList))+1]:=1;
  od;
  return TheIneq;
end;






OSMET__FacetInequalities:=function(n,m)
  local ListIneq, eSet, iIdx, RedList, PermRedList, ListElement, i, TheIneq, eList, eIneq, fIneq, iIneq, jIneq, GRP, ListGenNew, eGen, ListSymbol;
  ListIneq:=[];
  ListSymbol:=[];
  ListElement:=OSMET__ListElement(n,m);
  for eSet in Combinations([1..n],m+2)
  do
    for iIdx in eSet
    do
      RedList:=Difference(eSet, [iIdx]);
      Add(ListIneq, OSMET__FuncInequality(ListElement, RedList, iIdx));
      Add(ListSymbol, [RedList, iIdx]);
      PermRedList:=Permuted(RedList, (1,2));
      Add(ListIneq, OSMET__FuncInequality(ListElement, PermRedList, iIdx));
      Add(ListSymbol, [PermRedList, iIdx]);
    od;
  od;
  for i in [1..Length(ListElement)]
  do
    TheIneq:=ListWithIdenticalEntries(1+Length(ListElement), 0);
    TheIneq[i+1]:=1;
    Add(ListIneq, TheIneq);
  od;
  GRP:=OSMET__SymmetryGroup(n,m);
  ListGenNew:=[];
  for eGen in GeneratorsOfGroup(GRP)
  do
    eList:=[];
    for iIneq in [1..Length(ListIneq)]
    do
      eIneq:=ListIneq[iIneq];
      fIneq:=Permuted(eIneq, eGen);
      jIneq:=Position(ListIneq, fIneq);
      Add(eList, jIneq);
    od;
    Add(ListGenNew, PermList(eList));
  od;
  return rec(ListIneq:=ListIneq, GroupFac:=Group(ListGenNew));
end;
