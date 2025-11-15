Read("FunctionSet.g");
#
# there is a conjecture that every facet of CUTn is adjacent
# to a triangle inequality facet.
# we are searching a counterexample to this conjecture by doing
# random movement in the graph of facets

n:=10;
eFileSave:=Concatenation("DATA/RecordDesc", String(n));
eRec:=ReadAsFunction(eFileSave)();

eFileSave:=Concatenation("DATA/RecordFAC", String(n));
ListCand:=ReadAsFunction(eFileSave)();
WorkFac:=SelectMinimalLengthElement(ListCand);

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
      Print("We got one!!!, WorkFac=\n");
      Print(NullMat(5));
    else
      Print("We have a facet of incidence ", Length(TheLift), "; continuing\n");
      WorkFac:=TheLift;
    fi;
  fi;
od;

