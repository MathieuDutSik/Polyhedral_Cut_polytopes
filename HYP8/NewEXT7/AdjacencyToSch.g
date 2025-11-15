Read("Functions/MyPersoFunctions.g");
R:=ReadAsFunction("DATA/NewExt7")();
EXT:=R.EXT;
DM:=FuncDistMat(EXT);
SymGrp:=AutomorphismGroupEdgeColoredGraph(DM);

O:=Orbits(SymGrp, [1..35], OnPoints);

EXTsch:=Concatenation(EXT{O[1]}, EXT{O[3]});
AffBas:=CreateAffineBasis(EXTsch);
Add(AffBas, EXT[O[2][1]]);

EXT2:=HypermetricCoordinates(AffBas, Concatenation(EXT{O[1]}, 
EXT{O[3]}, EXT{O[2]}, EXT{O[4]}));


EXT3:=List(EXT2, x->Concatenation([1], x{[2..8]}));
EXT4:=EXT3{Difference([1..35], [28])};

EXTred:=EXT3{[1..27]};
TheMat:=MatrixHypermetricInequalities(EXTred);

LSEL:=Filtered(EXT4, x->x[8]<>0);
FindTrans:=function(eVect)
  local eW, eTrans, test, vVect, alpha;
  for eW in EXTred
  do
    eTrans:=eVect-eW;
    test:=true;
    for vVect in LSEL
    do
      alpha:=vVect[8];
      if Position(EXT3, vVect-alpha*eTrans)=fail then
        test:=false;
      fi;
    od;
    if test=true then
      return eTrans;
    fi;
  od;
end;
eTrans2:=ListWithIdenticalEntries(8,0);
eTrans2[8]:=1;

eTrans:=FindTrans([   1,   1,  -1,  -2,   1,   1,   0,   1 ]);
LPART:=List(LSEL, x->x-x[8]*eTrans);
LPos:=List(LPART, x->Position(EXT4, x));
EXT5:=Concatenation(EXT4{Difference([1..27], Set(LPos))}, EXT4{LPos});
for eVect in LSEL
do
  alpha:=eVect[8];
  Add(EXT5, eVect-alpha*eTrans+alpha*eTrans2);
od;

eVect:=EXT3[28];
alpha:=eVect[8];
OtherVect:=eVect-alpha*eTrans+alpha*eTrans2;
