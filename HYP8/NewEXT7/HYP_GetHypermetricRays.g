R:=ReadAsFunction("DATA/NewExt7")();
GramMat:=FuncExtreme(R.EXT).GramMatrix;
ListOrbAff:=ReadAsFunction("DATA/OrbitAffineBasis")();

ListOrbIndep:=ReadAsFunction("DATA/OrbitIndependentSet")();
ListDet:=List(ListOrbIndep, x->AbsInt(DeterminantMat(R.EXT{x})));


n:=7;
nbV:=n+1;

ListExtremeHypermetric:=[];
TotalNumber:=0;
nbOrb:=Length(ListOrbAff);
for iOrb in [1..nbOrb]
do
  Print("iOrb=", iOrb, " / ", nbOrb, "\n");
  eOrb:=ListOrbAff[iOrb];
  eDist:=NullMat(nbV,nbV);
  EXTsel:=R.EXT{eOrb};
  for i in [1..nbV-1]
  do
    for j in [i+1..nbV]
    do
      eDiff:=EXTsel[i]{[2..n+1]} - EXTsel[j]{[2..n+1]};
      eNorm:=eDiff*GramMat*eDiff;
      eDist[i][j]:=eNorm;
      eDist[j][i]:=eNorm;
    od;
  od;
  GRP:=AutomorphismGroupEdgeColoredGraph(eDist);
  eSize:=Factorial(nbV)/Order(GRP);
  TotalNumber:=TotalNumber + eSize;
  Add(ListExtremeHypermetric, eDist);
od;
FileSave:="DATA/ListOrbAffER7_total";
eRecSave:=rec(ListExtremeHypermetric:=ListExtremeHypermetric,
              ListOrbAff:=ListOrbAff);
SaveDataToFile(FileSave, eRecSave);
