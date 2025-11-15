PL:=ArchimedeanPolyhedra("Octahedron");
LFac:=__FaceSet(PL);

ListSimpIneq:=[];
for eFac in LFac
do
  eSimp:=List(eFac, x->x[1]);
  Add(ListSimpIneq, Set(eSimp));
od;

nbVert:=Length(PL);

RecHMET:=CMC_GetCanonicalSimplicialInequality(nbVert, 2);

nbSet:=Length(RecHMET.ListSets);

eIneq:=ListWithIdenticalEntries(nbSet,0);
ListPos:=[];
for eSimp in ListSimpIneq
do
  pos:=Position(RecHMET.ListSets, eSimp);
  Add(ListPos, pos);
od;
eIneq{ListPos}:=ListWithIdenticalEntries(Length(ListPos),1);
eIneq[ListPos[1]]:=-1;

eSol:=SolutionMatNonnegative(RecHMET.ListIneq, eIneq);
