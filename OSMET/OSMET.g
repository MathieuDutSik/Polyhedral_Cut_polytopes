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
  local U, iElt, GRP, ListGenNew, eGen, eList, jPart, iPart, ePart, fPart;
  U:=[[[1]]];
  for iElt in [2..n]
  do
    U:=IncreaseLevel(U, iElt);
  od;
  GRP:=SymmetricGroup(n);
  ListGenNew:=[];
  for eGen in GeneratorsOfGroup(GRP)
  do
    eList:=[];
    for iPart in [1..Length(U)]
    do
      ePart:=U[iPart];
      fPart:=OnTuplesSets(ePart, eGen);
      jPart:=Position(U, fPart);
      Add(eList, jPart);
    od;
    Add(ListGenNew, PermList(eList));
  od;
  eList:=[];
  for iPart in [1..Length(U)]
  do
    ePart:=U[iPart];
    fPart:=Reversed(ePart);
    jPart:=Position(U, fPart);
    Add(eList, jPart);
  od;
  Add(ListGenNew, PermList(eList));
  return rec(ListOrientedPartition:=U, GroupExt:=Group(ListGenNew));
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
      AppendTo(output, " ", OP[iP]);
    od;
    if iOP<Length(OPset) then
      AppendTo(output, " \\#");
    fi;
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



OSCUT__Description:=function(n,m)
  local HR, LELT, U, NewU, NewListGens, eGen, eList, iU, pos1, pos2, pos3, Lext;
  HR:=SpanOrientedPartGeneral(n);
  LELT:=OSMET__ListElement(n,m);
  U:=HR.ListOrientedPartition;
  NewU:=Filtered(U, x->Length(x)>m);
  NewListGens:=[];
  for eGen in GeneratorsOfGroup(HR.GroupExt)
  do
    eList:=[];
    for iU in [1..Length(NewU)]
    do
      pos1:=Position(U, NewU[iU]);
      pos2:=OnPoints(pos1, eGen);
      pos3:=Position(NewU, U[pos2]);
      Add(eList, pos3);
    od;
    Add(NewListGens, PermList(eList));
  od;
  Lext:=List(NewU, x->PartToVector(x, LELT));
  return rec(EXT:=Lext, Symbol:=NewU, GroupExt:=Group(NewListGens));
end;


FuncOrdPart:=function(eSol)
  local iPoint, eCase, iSize, Siz;
  iPoint:=0;
  eCase:=[];
  for iSize in [1..Length(eSol)]
  do
    Siz:=eSol[iSize];
    Add(eCase, [iPoint+1..iPoint+Siz]);
    iPoint:=iPoint+Siz;
  od;
  return eCase;
end;


ListTypeOrderedPartition:=function(n)
  local LE, Lreturn, eCase, eVect, i, j, k, GRP, eO, iPoint, Siz, eSol, LordPart, iSize;
  LE:=ListOfPartitions(n,n);
  Lreturn:=[];
  for eCase in LE
  do
    eVect:=[];
    for i in [1..Length(eCase)]
    do
      for j in [1..eCase[i]]
      do
        Add(eVect, i);
      od;
    od;
    k:=Length(eVect);
    GRP:=SymmetricGroup(k);
    for eO in Orbit(GRP, eVect, Permuted)
    do
      Add(Lreturn, eO);
    od;
  od;
  LordPart:=List(Lreturn, x->FuncOrdPart(x));
  return rec(Lreturn:=Lreturn, ListCases:=LordPart);
end;


TranslateGroupPlusMorphism:=function(GRP, LSET, Action)
  local eGen, ListGenNew, eList, ListGen, GRPnew;
  if Order(GRP)=1 then
    return Group(());
  fi;
  ListGenNew:=[];
  ListGen:=GeneratorsOfGroup(GRP);
  for eGen in ListGen
  do
    eList:=TranslateElement(eGen, LSET, Action);
    Add(ListGenNew, PermList(eList));
  od;
  GRPnew:=Group(ListGenNew);
  return rec(GRPnew:=GRPnew, phi:=GroupHomomorphismByImages(GRP, GRPnew, ListGen, ListGenNew));
end;



#
# return Group and list of extreme rays for
# a number of parts equal to m+1.
OSCUT_Family:=function(n,m)
  local GRP, ElT, ListCases, ListElement, ePart, ListEXT, Orb, Orb2, fPart, PermGRP, eVect, i, j, Lreturn, ListPart, LE, eCase, LordPart, Htrans;
  GRP:=SymmetricGroup(n);
  ListElement:=OSMET__ListElement(n,m);
  #
  LE:=ListOfPartitions(n,n);
  Lreturn:=[];
  for eCase in LE
  do
    eVect:=[];
    for i in [1..Length(eCase)]
    do
      for j in [1..eCase[i]]
      do
        Add(eVect, i);
      od;
    od;
    Add(Lreturn, eVect);
  od;
  LordPart:=List(Lreturn, x->FuncOrdPart(x));

  ListCases:=[];
  for ePart in LordPart
  do
    if Length(ePart)=m+1 then
      Add(ListCases, ePart);
    fi;
  od;
  ListPart:=[];
  for ePart in ListCases
  do
    Orb:=Orbit(GRP, ePart, OnSetsSets);
    Orb2:=List(Orb, x->Permuted(x, (1,2)));
    Append(ListPart, Orb);
    Append(ListPart, Orb2);
  od;
  ListEXT:=List(ListPart, x->PartToVector(x, ListElement));
  PermGRP:=OSMET__SymmetryGroup(n,m);
  Htrans:=TranslateGroupPlusMorphism(PermGRP, ListEXT, Permuted);
  return rec(ListPart:=ListPart, EXT:=ListEXT, GroupExt:=Htrans.GRPnew, phi:=Htrans.phi);
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


List2String:=function(eList)
  local eStr, i;
  eStr:="";
  for i in eList
  do
    eStr:=Concatenation(eStr, String(i));
  od;
  return eStr;
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
      Add(ListSymbol, Concatenation("OS(", List2String(RedList), ",", String(iIdx), ")"));
      PermRedList:=Permuted(RedList, (1,2));
      Add(ListIneq, OSMET__FuncInequality(ListElement, PermRedList, iIdx));
      Add(ListSymbol, Concatenation("OS(", List2String(PermRedList), ",", String(iIdx), ")"));
    od;
  od;
  for i in [1..Length(ListElement)]
  do
    TheIneq:=ListWithIdenticalEntries(1+Length(ListElement), 0);
    TheIneq[i+1]:=1;
    Add(ListIneq, TheIneq);
    Add(ListSymbol, Concatenation("NN(", List2String(ListElement[i]),")"));
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
#  eList:=[];
#  for iElt in [1..Length(ListElement)]
#  do
#    eSet:=ListElement[iElt];
#    ePermSet:=Permuted(eSet, (1,2));
#    Add(eList, Position(
#  od;
  return rec(ListIneq:=ListIneq, GroupFac:=Group(ListGenNew), ListSymbol:=ListSymbol);
end;



MatrixFromVector:=function(eVect, n, m)
  local ListElement, TheMat, len, iCol, eOrd, eDiff, U, SGN, W, i, j;
  if n<>m+3 then
    Print("Unconsistency between n and m found\n");
    Print(NullMat(5));
  fi;
  ListElement:=OSMET__ListElement(n,m);
  TheMat:=NullMat(n,n);
  len:=Length(eVect);
  if len<>Length(ListElement) then
    Print("Unconsistency\n");
    Print(NullMat(5));
  fi;
  for iCol in [1..len]
  do
    eOrd:=ListElement[iCol];
    eDiff:=Difference([1..n], Set(eOrd));
    U:=[];
    Append(U, eOrd);
    Append(U, eDiff);
    SGN:=SignPerm(PermList(U));
    if SGN=1 then
      W:=eDiff;
    else
      W:=Reversed(eDiff);
    fi;
    i:=W[1];
    j:=W[2];
    TheMat[i][j]:=eVect[iCol];
  od;
  return TheMat;
end;



SignVectorFromVector:=function(eVect, n, m)
  local ListElement, len, iCol, eOrd, eDiff, U, SGN, W, SignVect;
  if n<>m+2 then
    Print("Unconsistency between n and m found\n");
    Print(NullMat(5));
  fi;
  ListElement:=OSMET__ListElement(n,m);
  len:=Length(eVect);
  if len<>Length(ListElement) then
    Print("Unconsistency\n");
    Print(NullMat(5));
  fi;
  SignVect:=[];
  for iCol in [1..len]
  do
    eOrd:=ListElement[iCol];
    eDiff:=Difference([1..n], Set(eOrd));
    U:=[];
    Append(U, eOrd);
    Append(U, eDiff);
    SGN:=SignPerm(PermList(U));
    if SGN=1 then
      W:=[1, eDiff[1]];
    else
      W:=[-1, eDiff[1]];
    fi;
    Add(SignVect, [W, eVect[iCol]]);
  od;
  return SignVect;
end;





VertexSplitting:=function(n, m, eVect, i)
  local ListElement1, ListElement2, V, eElt, pos1, pos2, posElt, fElt, gElt;
  ListElement1:=OSMET__ListElement(n,m);
  ListElement2:=OSMET__ListElement(n+1,m);
  V:=[0];
  for eElt in ListElement2
  do
    pos1:=Position(n+1, eElt);
    pos2:=Position(i, eElt);
    if pos1<>fail and pos2<>fail then
      Add(V, 0);
    fi;
    if pos1=fail then
      posElt:=Position(ListElement1, eElt);
      Add(V, eVect[posElt+1]);
    fi;
    if pos1<>fail and pos2=fail then
      fElt:=ShallowCopy(eElt);
      fElt[pos1]:=i;
      gElt:=OSMET__Canonicalize(fElt);
      posElt:=Position(ListElement1, gElt);
      Add(V, eVect[posElt+1]);
    fi;
  od;
  return V;
end;

ZeroExtension:=function(n,m, eVect)
  local ListElement1, ListElement2, V, eElt, pos1, fElt, gElt, posElt;
  ListElement1:=OSMET__ListElement(n,m);
  ListElement2:=OSMET__ListElement(n+1,m+1);
  V:=[0];
  for eElt in ListElement2
  do
    pos1:=Position(n+1, eElt);
    if pos1=fail then
      Add(V, 0);
    else
      fElt:=OSMET__Canonicalize(eElt);
      gElt:=fElt{[1..m+1]};
      posElt:=Position(ListElement1, gElt);
      Add(V, eVect[posElt+1]);
    fi;
  od;
  return V;
end;





OSCUT__HintExpansionFormula:=function(n,m)
  local ListCases, W, ListElement, PLC, FAC, PermGRP, FuncTest, ePart, eVect, TheStab, TheImg, Lincd, Wcand, O, iOrb, eOrb, SVect, iElt, rnk, output;
  ListCases:=Filtered(ListTypeOrderedPartition(n).ListCases, x->Length(x)>m);
  W:=OSCUT_Family(n,m);
  ListElement:=OSMET__ListElement(n,m);
  PLC:=OSMET__FacetInequalities(n,m);
  FAC:=PLC.ListIneq;
  PermGRP:=OSMET__SymmetryGroup(n,m);
  FuncTest:=function(vList, eVect)
    local fVect;
    for fVect in vList
    do
      if fVect*eVect>0 then
        return false;
      fi;
    od;
    return true;
  end;
  output:=OutputTextFile(Concatenation("ExpansionFormula", String(n), "_", String(m)), true);
  for ePart in ListCases
  do
    eVect:=PartToVector(ePart, ListElement);
    TheStab:=Stabilizer(PermGRP, eVect, Permuted);
    TheImg:=Image(W.phi, TheStab);
    Print("|Stab|=", Order(TheStab), "  |Img|=", Order(TheImg), "\n");
    Lincd:=Filtered(FAC, x->x*eVect=0);
    rnk:=RankMat(Lincd);
    if rnk < Length(ListElement)-1 then
      Wcand:=Filtered([1..Length(W.EXT)], x->FuncTest(Lincd, W.EXT[x])=true);
      O:=Orbits(TheImg, Wcand, OnPoints);
      AppendTo(output, "ePart=", ePart, "\n");
      AppendTo(output, "eVect=", eVect, "\n");
      for iOrb in [1..Length(O)]
      do
        eOrb:=O[iOrb];
        AppendTo(output, "Orbit ", iOrb, "\n");
        SVect:=ListWithIdenticalEntries(Length(ListElement)+1, 0);
        for iElt in eOrb
        do
          SVect:=SVect+W.EXT[iElt];
          AppendTo(output, "incd=", W.ListPart[iElt], "\n");
        od;
        AppendTo(output, "SVect=", SVect, "\n");
      od;
      AppendTo(output, "\n");
    fi;
  od;
  CloseStream(output);
end;



OSCUT__Test:=function(n,m)
  local GRP, ElT, ePart, PLC, FAC, ListElement, k, ListStatus, eVect, ScalProd, TheStatus, i, Finc, rnk, ListCases;
  GRP:=SymmetricGroup(n);
  ElT:=ListTypeOrderedPartition(n);
  PLC:=OSMET__FacetInequalities(n,m);
  FAC:=PLC.ListIneq;
  ListElement:=OSMET__ListElement(n,m);
  ListCases:=[];
  for ePart in ElT.ListCases
  do
    k:=Length(ePart);
    if k>m then
      Add(ListCases, ePart);
    fi;
  od;
  ListStatus:=[];
  for ePart in ListCases
  do
    eVect:=PartToVector(ePart, ListElement);
    ScalProd:=List(FAC, x->x*eVect);
    TheStatus:="ray";
    for i in [1..Length(FAC)]
    do
      if ScalProd[i]<0 then
        TheStatus:="nonray";
      fi;
    od;
    if TheStatus="ray" then
      Finc:=Filtered([1..Length(FAC)], x->ScalProd[x] = 0);
      rnk:=RankMat(FAC{Finc});
      if rnk = Length(ListElement)-1 then
        TheStatus:="extreme ray";
      fi;
    fi;
    Add(ListStatus, TheStatus);
  od;
  return rec(ListStatus:=ListStatus, ListCases:=ListCases);
end;



PrintCutDimensionRayCases:=function(FileName, n)
  local m, output, W, nbCase, LRAY, LNRAY, LEXTRAY, iCase;
  output:=OutputTextFile(FileName, true);
  AppendTo(output, "\\section{ray status of oriented multicuts for dimension ", n, "}\n");
  for m in [1..n-2]
  do
    W:=OSCUT__Test(n,m);
    nbCase:=Length(W.ListStatus);
    AppendTo(output, "\\subsection{For m=", m, "}\n");
    LRAY:=Filtered([1..nbCase], x->W.ListStatus[x]="ray");
    LNRAY:=Filtered([1..nbCase], x->W.ListStatus[x]="nonray");
    LEXTRAY:=Filtered([1..nbCase], x->W.ListStatus[x]="extreme ray");
    if Length(LRAY)>0 then
      AppendTo(output, "Only ray\n");
      AppendTo(output, "\\begin{multicols}{3}\n");
      AppendTo(output, "\\noindent ");
      for iCase in LRAY
      do
        PartOutput(output, W.ListCases[iCase]);
        AppendTo(output, "\\\\\n");
      od;
      AppendTo(output, "\\end{multicols}\n");
    fi;
    if Length(LNRAY)>0 then
      AppendTo(output, "non ray\n");
      AppendTo(output, "\\begin{multicols}{3}\n");
      AppendTo(output, "\\noindent ");
      for iCase in LNRAY
      do
        PartOutput(output, W.ListCases[iCase]);
        AppendTo(output, "\\\\\n");
      od;
      AppendTo(output, "\\end{multicols}\n");
    fi;
    if Length(LEXTRAY)>0 then
      AppendTo(output, "Extreme ray\n");
      AppendTo(output, "\\begin{multicols}{3}\n");
      AppendTo(output, "\\noindent ");
      for iCase in LEXTRAY
      do
        PartOutput(output, W.ListCases[iCase]);
        AppendTo(output, "\\\\\n");
      od;
      AppendTo(output, "\\end{multicols}\n");
    fi;
  od;
  CloseStream(output);

end;






