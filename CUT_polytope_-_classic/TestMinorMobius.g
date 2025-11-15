

for m in [4..10]
do
  GRA:=MobiusLadderGraph(m);
  K5gra:=CompleteGraph(Group(()), 5);
  IsMinor:=SAGE_TestIfHasMinor(GRA, K5gra);
  Print("m=", m, " IsMinor=", IsMinor, "\n");
od;