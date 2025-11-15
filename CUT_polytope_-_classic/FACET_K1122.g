ListSize:=[1,1,2,2];
GRA:=GRAPH_GetMultiComplement(ListSize);


eRec:=CMC_GetCutPolytope(GRA);
RecCUT:=eRec.GetCUT_info();
EXTcut:=RecCUT.EXT;


GetHypIneq:=function(eList)
  local nbEdge, eIneq, len, i, j, ePair, pos;
  nbEdge:=Length(eRec.ListEdges);
  eIneq:=ListWithIdenticalEntries(1+nbEdge,0);
  len:=Length(eList);
  if len > GRA.order then
    Error("len is too big for us");
  fi;
  for i in [1..len]
  do
    for j in [i+1..len]
    do
      ePair:=[i,j];
      pos:=Position(eRec.ListEdges, ePair);
      if pos<>fail then
        eIneq[1+pos]:= - eList[i] * eList[j];
      fi;
    od;
  od;
  return eIneq;
end;


CheckIneq:=function(eIneq)
  local nbEdge, test, EXTface;
  nbEdge:=Length(eRec.ListEdges);
  test:=First(EXTcut, x->x*eIneq<0);
  if test<>fail then
    return false;
  fi;
  EXTface:=Filtered(EXTcut, x->x*eIneq=0);
  if RankMat(EXTface)<>nbEdge then
    return false;
  fi;
  return true;
end;


eIneq:=GetHypIneq([1,1,1,-1,-1,0]) + GetHypIneq([0,0,1,1,0,-1]);
test:=CheckIneq(eIneq);
Print("test=", test, "\n");
