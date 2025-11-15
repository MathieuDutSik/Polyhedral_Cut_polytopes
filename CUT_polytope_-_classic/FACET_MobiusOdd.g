k:=4;
m:=2*k+1;

GRA:=MobiusLadderGraph(m);

eRec:=CMC_GetCutPolytope(GRA);

nbEdge:=Length(eRec.ListEdges);
eIneq:=ListWithIdenticalEntries(1+nbEdge,1);
eIneq[1]:=0;


ListPair:=[rec(val:=1, pair:=[2*m-1,2*m]),
           rec(val:=1, pair:=[2,2*m-1]),
           rec(val:=-2, pair:=[1,2*m])];
for i in [1..k]
do
  Add(ListPair, rec(val:=2, pair:=[4*i-2,4*i]));
  Add(ListPair, rec(val:=2, pair:=[4*i-1,4*i+1]));
  Add(ListPair, rec(val:=1, pair:=[4*i-3,4*i-2]));
  Add(ListPair, rec(val:=-1, pair:=[4*i-1,4*i]));
  Add(ListPair, rec(val:=1, pair:=[4*i-3,4*i-1]));
  Add(ListPair, rec(val:=1, pair:=[4*i,4*i+2]));
od;

for eRecB in ListPair
do
  pos:=Position(eRec.ListEdges, eRecB.pair);
  eIneq[pos+1]:=eRecB.val;
od;

RecCUT:=eRec.GetCUT_info();
EXTcut:=RecCUT.EXT;
test:=First(EXTcut, x->x*eIneq< 0);
ListScal:=List(EXTcut, x->x*eIneq);
EXTface:=Filtered(EXTcut, x->x*eIneq=0);
rnk:=RankMat(EXTface);
Print("m=", m, " test=", test, " len=", Length(EXTcut[1]), " rnk=", rnk, " |EXTface|=", Length(EXTface), "\n");
