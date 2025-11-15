Read("Functions/MyPersoFunctions.g");
#
EXT56:=ReadAsFunction("Functions/TheGosset")();;
for i in [1..56]
do
  EXT56[i][1]:=1;
od;
DM56:=FuncDistMat(EXT56);
GRP56:=AutomorphismGroupEdgeColoredGraph(DM56);
GramMat:=FuncExtreme(EXT56).GramMatrix;


ListOrbGosset:=ReadAsFunction("DATA/ListOrbIndepGosset")();
ListDetEmpty:=List(ListOrbGosset.ListEmpty, x->AbsInt(DeterminantMat(EXT56{x})));
ListDetNonEmpty:=List(ListOrbGosset.ListNonEmpty, x->AbsInt(DeterminantMat(EXT56{x})));



ListHypermetric:=[];
nbOrb:=Length(ListOrbGosset.ListNonEmpty);
ListExtremeHypermetric:=[];
TotalNumber:=0;

n:=7;
nbV:=n+1;

ListOrbAff:=[];
for iOrb in [1..nbOrb]
do
  eOrb:=ListOrbGosset.ListNonEmpty[iOrb];
  eDet:=ListDetNonEmpty[iOrb];
  if eDet=1 then
    eDist:=NullMat(nbV,nbV);
    EXTsel:=EXT56{eOrb};
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
    Add(ListOrbAff, eOrb);
  fi;
od;

RecInfo:=rec(EXT:=EXT56,
             ListExtremeHypermetric:=ListExtremeHypermetric,
             TotalNumber:=TotalNumber, 
             ListOrbAff:=ListOrbAff);
SaveDataToFile("DATA/ListRecAffGosset", RecInfo);
