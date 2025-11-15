#

for n in [3..9]
do
  eRec:=ConePoly_PartialMetrics(n);
  FAC:=eRec.ListIneqPoly;
  GRP:=LinPolytope_Automorphism(FAC);
  BoundingSet:=[];
  eSkel:=SkeletonGraph(GRP, FAC, BoundingSet);
  eDiam:=Diameter(eSkel);
  Print("n=", n, " diam=", eDiam, "\n");


od;
