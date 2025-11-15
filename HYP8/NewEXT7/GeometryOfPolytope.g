Read("Functions/MyPersoFunctions.g");
R:=ReadAsFunction("DATA/NewExt7")();
EXT:=R.EXT;
DM:=FuncDistMat(EXT);
SymGrp:=AutomorphismGroupEdgeColoredGraph(DM);


O:=Orbits(SymGrp, [1..35], OnPoints);

EXTbyOrbit:=[];
for eO in O
do
  VC:=[];
  for iPos in eO
  do
    Add(VC, EXT[iPos]);
  od;
  Add(EXTbyOrbit, VC);
  Print(VC, "\n");
od;

FAC:=R.FAC;

OrbPref:=EXTbyOrbit[1];
Plane1:=NullspaceMat(TransposedMat(OrbPref))[1];


#for eFac in FAC
#do
#  SPC:=0;
#  for eExt in OrbPref
#  do
#    SPC:=SPC+eFac*eExt;
#  od;
#  Print("SPC=", SPC, "\n");
#od;

SymFac:=AutGroupGraph(R.GraphFac);

ListO:=Orbits(SymFac, [1..228], OnPoints);
for iOrb in [1..Length(ListO)]
do
  Print("iOrb=", iOrb, "\n");
  eOrb:=ListO[iOrb];
  eFac:=FAC[eOrb[1]];
  Print("eFac=", eFac, "\n");

  ExtSet:=[];
  for iExt in [1..Length(EXT)]
  do
    if EXT[iExt]*eFac=0 then
      Add(ExtSet, iExt);
    fi;
  od;
  IndGra:=InducedSubgraph(R.GraphExt, ExtSet);
  Stab:=Stabilizer(SymGrp, ExtSet, OnSets);
  Adja:=Adjacency(R.GraphFac, eOrb[1]);
  for jOrb in [1..Length(ListO)]
  do
    INTS:=Intersection(Set(Adja), ListO[jOrb]);
    if Length(INTS)>0 then
      Print("Intersection of size ", Length(INTS), " with orbit ", jOrb, "\n");
    fi;
  od;


  Print("OrbitSize=", Length(eOrb), "\n");
  Print("nbExt=", Length(ExtSet), "\n");
  Print("Grp =", Order(Stab), "\n");
  Print("Grp2=", Order(AutGroupGraph(IndGra)), "\n");
  Print("\n");
od;


for i in [1..4]
do
  SetScal:=[];
  for eExt in EXTbyOrbit[i]
  do
    AddSet(SetScal, eExt*Plane1);
  od;
  Print("i=", i, "  SetScal=", SetScal, "\n");
od;


OnPLANE1:=Union(O[1], O[3]);
VE:=InducedSubgraph(R.GraphExt, OnPLANE1);
GRVE:=AutGroupGraph(VE);
O1new:=[];
for eR in O[1]
do
  Add(O1new, Position(VE.names, eR)); 
od;

StabO1:=Stabilizer(GRVE, O1new, OnSets);

UnProj:=[];
for iPos in [1..Length(OnPLANE1)]
do
  jPos:=OnPLANE1[iPos];
  Add(UnProj, EXT[jPos]);
od;

Proj:=[];

for eV in UnProj
do
  V:=[];
  for i in [2..8]
  do
    Add(V, eV[i]);
  od;
  Add(Proj, V);
od;

DDA:=DualDescriptionAdjacencies(Proj);

ListFAC6:=[];
ListFAC6ByIncidence:=[];
for i in [1..Length(DDA.FAC)]
do
  if Length(Adjacency(DDA.RidgeGraph, i))=6 then
    Add(ListFAC6, DDA.FAC[i]);
    Linc:=[];
    for iExt in [1..Length(Proj)]
    do
      if Proj[iExt]*DDA.FAC[i]=0 then
        Add(Linc, iExt);
      fi;
    od;
    Add(ListFAC6ByIncidence, Linc);
  fi;
od;

ListNormalizer:=[];
ListStatus:=[];
DualityList:=[];
ListSet:=[];
for iFac in [1..Length(ListFAC6ByIncidence)]
do
  StabFac:=Stabilizer(GRVE, ListFAC6ByIncidence[iFac], OnSets);
  NorM:=Normalizer(GRVE, StabFac);
  Add(ListNormalizer, NorM);
  Add(ListStatus, 1);
  Add(DualityList, 0);
  Add(ListSet, []);
od;

for iFac in [1..Length(ListFAC6ByIncidence)]
do
  if ListStatus[iFac]=1 then
    GRP:=ListNormalizer[iFac];
    test:=1;
    for jFac in [iFac+1..Length(ListFAC6ByIncidence)]
    do
      if test=1 then
        if ListNormalizer[jFac]=GRP then
          test:=0;
          DualityList[iFac]:=jFac;
          DualityList[jFac]:=iFac;
          ListStatus[iFac]:=0;
          ListStatus[jFac]:=0;
          Hset:=Union(ListFAC6ByIncidence[iFac], ListFAC6ByIncidence[jFac]);
          ListSet[iFac]:=Hset;
          ListSet[jFac]:=Hset;
        fi;
      fi;
    od;
  fi;
od;


pos1:=Position(ListSet, Set(O1new));
pos2:=DualityList[pos1];

Set1:=[];
for eVert in ListFAC6ByIncidence[pos1]
do
  Add(Set1, VE.names[eVert]);
od;
Set2:=[];
for eVert in ListFAC6ByIncidence[pos2]
do
  Add(Set2, VE.names[eVert]);
od;

LAD:=[];
for i in [1..6]
do
  LA:=[];
  for j in [1..6]
  do
    if Distance(R.GraphExt, Set1[i], Set2[j])=1 then
      Add(LA, Set2[j]);
    fi;
  od;
  Add(LAD, LA);
od;

LAD2:=[];
for i in [1..6]
do
  LA:=[];
  for j in [1..6]
  do
    if Distance(R.GraphExt, Set1[i], O[2][j])=2 then
      Add(LA, O[2][j]);
    fi;
  od;
  if Length(LA)>2 then
    Print("there is a problem\n");
  fi;
  Add(LAD2, LA[1]);
od;


Scal:=0;
for i in [1..6]
do
  for j in [1..6]
  do
    VEC1:=EXT[Set1[i]]-EXT[Set1[j]];
    VEC2:=EXT[LAD2[i]]-EXT[LAD2[j]];
    RES:=VEC1+VEC2;
    Scal:=Scal+RES*RES;
  od;
od;
