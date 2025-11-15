Read("Functions/BaranovskiSearch");
Read("Functions/MyPersoFunctions.g");
Read("Functions/PermutedImprovedFunc.g");
Read("Functions/RecognitionDelaunay.g");



ListHypSymbol:=[];
ListHypIneq:=[];

ComputeHypermetricInformations:=function()
  local n, GRP, LORB, V, TheORB;
  n:=8;
  GRP:=SymmetricGroup(n);
  LORB:=GetConjecturalListFacetHYPn(n);
  for V in LORB
  do
    TheORB:=Orbit(GRP, V, Permuted);
    Append(ListHypSymbol, TheORB);
  od;
  ListHypIneq:=List(ListHypSymbol, FromHypermetricVectorToHypermetricFace);
end;

Print("Before determination of all facets\n");
ComputeHypermetricInformations();
Print("After determination of all facets\n");





EXT56:=
[ [ 1, 0, 0, 0, 0, 0, 0, 0 ], [ 1, 0, 0, 0, 0, 0, 1, 1 ],
  [ 1, 0, 0, 0, 0, 1, 0, 1 ], [ 1, 0, 0, 0, 0, 1, 1, 0 ],
  [ 1, 0, 0, 0, 1, 0, 0, 1 ], [ 1, 0, 0, 0, 1, 0, 1, 0 ],
  [ 1, 0, 0, 0, 1, 1, 0, 0 ], [ 1, 0, 0, 0, 1, 1, 1, 1 ],
  [ 1, 0, 0, 1, 0, 0, 0, 1 ], [ 1, 0, 0, 1, 0, 0, 1, 0 ],
  [ 1, 0, 0, 1, 0, 1, 0, 0 ], [ 1, 0, 0, 1, 0, 1, 1, 1 ],
  [ 1, 0, 0, 1, 1, 0, 0, 0 ], [ 1, 0, 0, 1, 1, 0, 1, 1 ],
  [ 1, 0, 0, 1, 1, 1, 0, 1 ], [ 1, 0, 0, 1, 1, 1, 1, 0 ],
  [ 1, 0, 1, 0, 0, 0, 0, 1 ], [ 1, 0, 1, 0, 0, 0, 1, 0 ],
  [ 1, 0, 1, 0, 0, 1, 0, 0 ], [ 1, 0, 1, 0, 0, 1, 1, 1 ],
  [ 1, 0, 1, 0, 1, 0, 0, 0 ], [ 1, 0, 1, 0, 1, 0, 1, 1 ],
  [ 1, 0, 1, 0, 1, 1, 0, 1 ], [ 1, 0, 1, 0, 1, 1, 1, 0 ],
  [ 1, 0, 1, 1, 0, 0, 0, 0 ], [ 1, 0, 1, 1, 0, 0, 1, 1 ],
  [ 1, 0, 1, 1, 0, 1, 0, 1 ], [ 1, 0, 1, 1, 0, 1, 1, 0 ],
  [ 1, 0, 1, 1, 1, 0, 0, 1 ], [ 1, 0, 1, 1, 1, 0, 1, 0 ],
  [ 1, 0, 1, 1, 1, 1, 0, 0 ], [ 1, 0, 1, 1, 1, 1, 1, 1 ],
  [ 1, 1, 1, 0, 0, 0, 0, 0 ], [ 1, -1, 0, 1, 1, 1, 1, 1 ],
  [ 1, 1, -1, 0, 0, 0, 0, 0 ], [ 1, -1, 2, 1, 1, 1, 1, 1 ],
  [ 1, 1, 0, 1, 0, 0, 0, 0 ], [ 1, -1, 1, 0, 1, 1, 1, 1 ],
  [ 1, 1, 0, -1, 0, 0, 0, 0 ], [ 1, -1, 1, 2, 1, 1, 1, 1 ],
  [ 1, 1, 0, 0, 1, 0, 0, 0 ], [ 1, -1, 1, 1, 0, 1, 1, 1 ],
  [ 1, 1, 0, 0, -1, 0, 0, 0 ], [ 1, -1, 1, 1, 2, 1, 1, 1 ],
  [ 1, 1, 0, 0, 0, 1, 0, 0 ], [ 1, -1, 1, 1, 1, 0, 1, 1 ],
  [ 1, 1, 0, 0, 0, -1, 0, 0 ], [ 1, -1, 1, 1, 1, 2, 1, 1 ],
  [ 1, 1, 0, 0, 0, 0, 1, 0 ], [ 1, -1, 1, 1, 1, 1, 0, 1 ],
  [ 1, 1, 0, 0, 0, 0, -1, 0 ], [ 1, -1, 1, 1, 1, 1, 2, 1 ],
  [ 1, 1, 0, 0, 0, 0, 0, 1 ], [ 1, -1, 1, 1, 1, 1, 1, 0 ],
  [ 1, 1, 0, 0, 0, 0, 0, -1 ], [ 1, -1, 1, 1, 1, 1, 1, 2 ] ];
DM56:=FuncDistMat(EXT56);
GRP56:=AutomorphismGroupEdgeColoredGraph(DM56);
DDA56:=DualDescriptionAdjacencies(EXT56);
GRA56:=DDA56.SkeletonGraph;


ListSubstructures:=[];
FuncInsertSubstructure:=function(TheSub, TheOrb)
  local TheGraph1, test, eStruct;
  TheGraph1:=CombinatorialModelGraphSubset(GRA56, TheSub);
  for eStruct in ListSubstructures
  do
    test:=EquivalenceVertexColoredGraph(TheGraph1.GRA, eStruct.GRA, TheGraph1.ThePartition);
    if test<>false then
      return;
    fi;
  od;
  if RankPolytope(EXT56{TheSub})<>2 then
    Print("We have a bug somewhere in the substructure search\n");
    Print(NullMat(5));
  fi;
  Print("We get a new one\n");
  Add(ListSubstructures, rec(TheSub:=TheSub, TheOrb:=TheOrb, GRA:=TheGraph1.GRA));
  if Length(ListSubstructures)=4 then
    Print("We are presumably finished\n");
    Print(NullMat(5));
  fi;
end;




FindAdjacentDelaunay:=function(TheOrb)
  local DM, DistVect, ListIncidentIndexes, ListIncidentIneq, GRP, TheStab, ListPermGens, eGen, eList, eFacet, PermGRP, BF, Ladj, IsRespawn, IsBankSave, TmpDir, DataPolyhedral, FuncStabilizer, FuncIsomorphy, WorkingDim, eAdj, iAdj, LS, Hcoord, ListIneqCoord, TheSub;
  DM:=List(DM56{TheOrb}, x->x{TheOrb});
  WorkingDim:=27;
  DistVect:=DistanceMatrixToDistanceVector(DM);
  ListIncidentIndexes:=Filtered([1..Length(ListHypIneq)], x->ListHypIneq[x]*DistVect=0);
  ListIncidentIneq:=ListHypIneq{ListIncidentIndexes};
  #
  GRP:=SymmetricGroupOnPairs(8);
  TheStab:=__PermutedStabilizer(GRP, DistVect);
  Print("|Incd|=", Length(ListIncidentIneq), "  |Stab|=", Order(TheStab), "\n");
  ListPermGens:=[];
  for eGen in TheStab
  do
    eList:=[];
    for eFacet in ListIncidentIneq
    do
      Add(eList, Position(ListIncidentIneq, Permuted(eFacet, eGen)));
    od;
    Add(ListPermGens, PermList(eList));
  od;
  PermGRP:=Group(ListPermGens);
  Print("We computed PermGRP\n");
  #
  #
  #
  FuncStabilizer:=function(EXTask)
    local GRP, ListGen, eGen, ListInc;
    ListInc:=List(EXTask, x->Position(ListIncidentIneq, x));
    GRP:=Stabilizer(PermGRP, ListInc, OnSets);
    ListGen:=[];
    for eGen in GeneratorsOfGroup(GRP)
    do
      Add(ListGen, PermList(List(ListInc, x->Position(ListInc, OnPoints(x, eGen)))));
    od;
    return PersoGroupPerm(ListGen);
  end;
  FuncIsomorphy:=function(EXT1, EXT2)
    local ePerm, ListInc1, ListInc2;
    if Length(EXT1)<>Length(EXT2) then
      return false;
    fi;
    ListInc1:=List(EXT1, x->Position(ListIncidentIneq, x));
    ListInc2:=List(EXT2, x->Position(ListIncidentIneq, x));
    ePerm:=RepresentativeAction(PermGRP, ListInc1, ListInc2, OnSets);
    if ePerm=fail then
      return false;
    else
      return PermList(List(ListInc1, x->Position(ListInc2, OnPoints(x, ePerm))));
    fi;
  end;
  BF:=BankRecording(rec(Saving:=false, BankPath:="/irrelevant/"), FuncStabilizer, FuncIsomorphy, OnSetsGroupFormalism());
  Print("We are here 1\n");
  #
  #
  #
  IsRespawn:=function(OrdGRP, EXT, TheDepth)
    if OrdGRP>=50000 and TheDepth<=2 then
      return true;
    fi;
    if OrdGRP<100 then
      return false;
    fi;
    if Length(EXT)<WorkingDim+7 then
      return false;
    fi;
    if TheDepth=2 then
      return false;
    fi;
    return true;
  end;
  IsBankSave:=function(EllapsedTime, OrdGRP, EXT, TheDepth)
    if TheDepth=0 then
      return false;
    fi;
    if EllapsedTime>=600 then
      return true;
    fi;
    if Length(EXT)<=WorkingDim+5 then
      return false;
    fi;
    return true;
  end;
  TmpDir:=Filename(POLYHEDRAL_tmpdir, "");
  DataPolyhedral:=rec(IsBankSave:=IsBankSave, 
        TheDepth:=0,
        IsRespawn:=IsRespawn, 
        Saving:=false,
        ThePathSave:="/irrelevant/",
        ThePath:=TmpDir,
        DualDescriptionFunction:=__DualDescriptionLRS_Reduction, 
        DualDescriptionFunction:=__DualDescriptionLRS_splitter,
        GroupFormalism:=OnSetsGroupFormalism());
  #
  #
  #
  Print("We are here 2\n");
#  Ladj:=__ListFacetByAdjacencyDecompositionMethod(ListIncidentIneq, PermGRP, DataPolyhedral, BF);
  Ladj:=__DualDescriptionCDD_Reduction(ListIncidentIneq, PermGRP, Filename(POLYHEDRAL_tmpdir, ""));
  Print("|Ladj|=", Length(Ladj), "\n");
  Print("We are here 3\n");
  for iAdj in [1..Length(Ladj)]
  do
    eAdj:=Ladj[iAdj];
    LS:=NullspaceMat(TransposedMat(ListIncidentIneq{eAdj}));
    Hcoord:=HypermetricCoordinates(EXT56{TheOrb}, EXT56);
    ListIneqCoord:=List(Hcoord, FromHypermetricVectorToHypermetricFace);
    TheSub:=Filtered([1..Length(ListIneqCoord)], x->ForAny(LS, y->y*ListIneqCoord[x]<>0)=false);
    FuncInsertSubstructure(TheSub, TheOrb);
    Print("iAdj=", iAdj, " Done\n");
  od;
end;


while(true)
do
  AffBas:=CreateAffineBasis(EXT56);
  eList:=List(AffBas, x->Position(EXT56, x));
  FindAdjacentDelaunay(eList);
od;

