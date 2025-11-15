Read("Functions/BaranovskiSearch");
Read("Functions/MyPersoFunctions.g");
Read("Functions/PermutedImprovedFunc.g");
Read("Functions/RecognitionDelaunay.g");

EXT:=ReadAsFunction("Functions/TheGosset")();
DistEXT:=FuncDistMat(EXT);

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

ComputeHypermetricInformations();



TheEXT35:=ReadAsFunction("Functions/TheEXT35")();;
TheGosset:=ReadAsFunction("Functions/TheGosset")();;
DM_EXT35:=FuncDistMat(TheEXT35);
DM_Gosset:=FuncDistMat(TheGosset);

IsInKnownList:=function(EXT2test)
  local DM;
  DM:=FuncDistMat(EXT2test);
  if IsIsomorphicEdgeColoredGraph(DM, DM_EXT35)<>false then
    Print("recognized as ErdahlRybnikov EXT35\n");
    return true;
  fi;
  if IsIsomorphicEdgeColoredGraph(DM, DM_Gosset)<>false then
    Print("recognized as Gosset G7 polytope\n");
    return true;
  fi;
  return false;
end;






FindAdjacentDelaunay:=function(TheOrb, VertSubset)
  local DM, DistVect, ListIncidentIndexes, ListIncidentIneq, GRP, TheStab, ListPermGens, eGen, eList, eFacet, PermGRP, RPLift, BF, Ladj, IsRespawn, IsBankSave, TmpDir, DataPolyhedral, FuncStabilizer, FuncIsomorphy, WorkingDim, eAdj, iAdj, ListAdj, ListTotalEXT, TheDistVect, TheLift, TheDistMat, TheGramMat, test, LS, Hcoord, ListIneqCoord, TheSub, RealListVertices;
  DM:=List(DistEXT{TheOrb}, x->x{TheOrb});
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
  RPLift:=__ProjectionLiftingFramework(ListHypIneq, ListIncidentIndexes);
  ListAdj:=[];
  for iAdj in [1..Length(Ladj)]
  do
    eAdj:=Ladj[iAdj];
    LS:=NullspaceMat(TransposedMat(ListIncidentIneq{eAdj}));
    Hcoord:=HypermetricCoordinates(EXT{TheOrb}, EXT);
    ListIneqCoord:=List(Hcoord, FromHypermetricVectorToHypermetricFace);
    TheSub:=Filtered([1..Length(ListIneqCoord)], x->ForAny(LS, y->y*ListIneqCoord[x]<>0)=false);
#    Print("Length(TheSub)=", Length(TheSub), "\n");
#    Print("TheSub=", TheSub, "\n");
    if TheSub=VertSubset then
      TheLift:=RPLift.FuncLift(eAdj);
      TheDistVect:=__FindFacetInequality(ListHypIneq, TheLift);
      TheDistMat:=DistanceVectorToDistanceMatrix(TheDistVect, 8);
      TheGramMat:=DistanceMatrixToGramMatrix(TheDistMat);
      if RankMat(TheGramMat)=7 then
        ListTotalEXT:=IdentityMat(8);
        Append(ListTotalEXT, ListHypSymbol{TheLift});
        RealListVertices:=ListOfVertices(TheDistVect);
        if Difference(Set(ListTotalEXT), Set(RealListVertices))<>[] then
          Print("Some facet ineq do not appear as vertex\n");
          Print(NullMat(5));
        fi;
        test:=IsInKnownList(RealListVertices);
        if test=false then
          Print("It is not in known list !!!!!!\n");
          Print("Maybe we missed a facet\n");
          Print("Or we find another extreme Delaunay in dimension 7\n");
          Print(NullMat(5));
        fi;
      fi;
      Add(ListAdj, TheDistVect);
    fi;
    Print("iAdj=", iAdj, " Done\n");
  od;
  return ListAdj;
end;


GetSpecialAffineBasis:=function(TheSub)
  local EXTred, AffBasis, IsOK, eVert, B;
  EXTred:=EXT{TheSub};
  while(true)
  do
    AffBasis:=CreateAffineBasis(EXTred);
    IsOK:=true;
    for eVert in EXT
    do
      B:=SolutionMat(AffBasis, eVert);
      if ForAny(B, x->IsInt(x)=false) then
        IsOK:=false;
      fi;
    od;
    if IsOK=true then
      return AffBasis;
    fi;
  od;
end;

ListSubstr:=ReadAsFunction("Rank2subsetGosset")();

ListResult:=[];
for MySub in ListSubstr
do
  AffBasis:=GetSpecialAffineBasis(MySub);
  eOrb:=List(AffBasis, x->Position(EXT, x));
  TheResult:=FindAdjacentDelaunay(eOrb, MySub);
  Add(ListResult, TheResult[1]);
od;

