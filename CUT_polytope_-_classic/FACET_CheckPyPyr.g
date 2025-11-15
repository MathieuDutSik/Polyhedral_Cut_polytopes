k:=2;
m:=2*k+1;
GRA1:=GRAPH_Cycle(m);
GRA2:=GRAPH_Pyramid(GRA1);
GRA:=GRAPH_Pyramid(GRA2);

eRec:=CMC_GetCutPolytope(GRA);

nbEdge:=Length(eRec.ListEdges);
eIneq:=ListWithIdenticalEntries(1+nbEdge,1);
eIneq[1]:=0;

ListPairs:=[[m-1,m],[m+1,m+2]];
for i in [1..k]
do
  Add(ListPairs, [2*i-1,m+1]);
  Add(ListPairs, [2*i-1,m+2]);
od;
for ePair in ListPairs
do
  pos:=Position(eRec.ListEdges, ePair);
  eIneq[pos+1]:=-1;
od;

RecCUT:=eRec.GetCUT_info();
EXTcut:=RecCUT.EXT;
test:=First(EXTcut, x->x*eIneq< 0);
ListScal:=List(EXTcut, x->x*eIneq);
EXTface:=Filtered(EXTcut, x->x*eIneq=0);
rnk:=RankMat(EXTface);
Print("m=", m, " test=", test, " len=", Length(EXTcut[1]), " rnk=", rnk, " |EXTface|=", Length(EXTface), "\n");
