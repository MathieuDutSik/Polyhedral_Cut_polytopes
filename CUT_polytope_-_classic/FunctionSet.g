SelectMinimalLengthElement:=function(ListCand)
  local eList, eMin, pos;
  eList:=List(ListCand, Length);
  eMin:=Minimum(eList);
  pos:=Position(eList, eMin);
  return ListCand[pos];
end;

GetRecordForOperations:=function(n)
  local EXT, FAC_Met, ListIncd, RnkEXT, FirstFac;
  EXT:=Set(CUT_Vertex(n, "vertex").EXT);
  FAC_Met:=MET_Facets(n, "polytope");
  RnkEXT:=RankMat(EXT);
  FirstFac:=Filtered([1..Length(EXT)], x->EXT[x]*FAC_Met[1]=0);
  return rec(EXT:=EXT,
             FAC_Met:=FAC_Met,
             RnkEXT:=RnkEXT,
             FirstFac:=FirstFac);
end;


IsAdjacentToMetFacet:=function(eRec, ListInc)
  local EXTinc, INTS, rnk, eFAC;
  EXTinc:=Set(eRec.EXT{ListInc});
  for eFAC in eRec.FAC_Met
  do
    INTS:=Filtered(EXTinc, x->x*eFAC=0);
    rnk:=RankMat(INTS);
    if rnk=eRec.RnkEXT-2 then
      return true;
    fi;
  od;
  return false;
end;

RandomWalkSearchCounterexample:=function(n)
  local eFileSave, eRec, ListCand, WorkFac, Lridge, RPLift, TheChoice, TheLift, test, ListFacet;
  eFileSave:=Concatenation("DATA/RecordDesc", String(n));
  eRec:=ReadAsFunction(eFileSave)();
  #
  eFileSave:=Concatenation("DATA/RecordFAC", String(n));
  ListCand:=ReadAsFunction(eFileSave)();
  WorkFac:=SelectMinimalLengthElement(ListCand);
  #
  ListFacet:=[];
  while(true)
  do
    Lridge:=GetInitialRays_LinProg(eRec.EXT{WorkFac}, 200);
    RPLift:=__ProjectionLiftingFramework(eRec.EXT, WorkFac);
    TheChoice:=SelectMinimalLengthElement(Lridge);
    TheLift:=RPLift.FuncLift(TheChoice);
    test:=IsAdjacentToMetFacet(eRec, TheLift);
    Add(ListFacet, TheLift);
    Print("Collected(ListLen)=", Collected(List(ListFacet, Length)), "\n");
    if Length(TheLift) < Length(WorkFac)+20 then
      if test=false then
        return WorkFac;
      else
        Print("We have a facet of incidence ", Length(TheLift), "; continuing\n");
        WorkFac:=TheLift;
      fi;
    fi;
  od;
end;


ExhaustiveSearchFromLowIncidence:=function(n, ListFacet)
  local eFileSave, eRec, TotalListFacet, TotalListStatus, FuncInsert, eFac, MinLen, i, TheLen, idxSel, TheFacet, Lridge, RPLift, eRidge, TheLift;
  eFileSave:=Concatenation("DATA/RecordDesc", String(n));
  eRec:=ReadAsFunction(eFileSave)();
  TotalListFacet:=[];
  TotalListStatus:=[];
  FuncInsert:=function(eFac)
    local test;
    test:=IsAdjacentToMetFacet(eRec, eFac);
    if test=false then
      Print("Victory!!!\n");
      Print(NullMat(5));
    fi;
    if Length(eFac)<70 then
      if Position(TotalListFacet, eFac)=fail then
        Add(TotalListFacet, eFac);
        Add(TotalListStatus, 1);
        Print("Now inserting |eFac|=", Length(eFac), " |FAC|=", Length(TotalListFacet), "\n");
      fi;
    fi;
  end;
  for eFac in ListFacet
  do
    FuncInsert(eFac);
  od;
  while(true)
  do
    MinLen:=2^n;
    for i in [1..Length(TotalListFacet)]
    do
      if TotalListStatus[i]=1 then
        TheLen:=Length(TotalListFacet[i]);
        if TheLen<MinLen then
          MinLen:=TheLen;
          idxSel:=i;
        fi;
      fi;
    od;
    TheFacet:=TotalListFacet[idxSel];
    TotalListStatus[idxSel]:=0;
    Print("Treating idxSel=", idxSel, " inc=", Length(TheFacet), "\n");
    RPLift:=__ProjectionLiftingFramework(eRec.EXT, TheFacet);
    Lridge:=DualDescriptionSets(eRec.EXT{TheFacet});
    for eRidge in Lridge
    do
      TheLift:=RPLift.FuncLift(eRidge);
      FuncInsert(TheLift);
    od;
  od;
end;
