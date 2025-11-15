Read("Functions/MyPersoFunctions.g");
#
# eCase:=rec(
#   Basis   the basis of the space considered
#   SuperMat   a positive definite matrix belonging to the space considered
#              such that EXT is a Delaunay for this matrix
#   TheGroup   the symmetry group of the case considered.
#              it has to satisfy 
#        U*M*TransposedMat(U)=M for any matrix M in the basis
#   EXT is a set of vertices invariant under GRP.
#     we want to find the linear inequality on elements
#     of the basis which characterize it being a Delaunay
# )
PermGroupToMatrixGroup:=function(EXT, PermGRP)
  return Group(List(GeneratorsOfGroup(PermGRP), x->__LemmaFindTransformation(EXT, EXT, x)));
end;


MatrixGroupToPermGroup:=function(EXT, MatrGRP)
  return Group(List(GeneratorsOfGroup(MatrGRP), x->PermList(List(EXT, y->Position(EXT, y*x)))));
end;



MyOrbitGeneration:=function(TheSet, GRPmat)
  local WorkSet, TheRet, eVert, O;
  WorkSet:=[];
  Append(WorkSet, TheSet);
  TheRet:=[];
  while(true)
  do
    eVert:=WorkSet[1];
    O:=Orbit(GRPmat, eVert);
    Append(TheRet, O);
    WorkSet:=Difference(WorkSet, Set(O));
    if Length(WorkSet)=0 then
      break;
    fi;
  od;
  return TheRet;
end;


SetReductionPermGroup:=function(EXT, GroupEXT)
  local EXTred, ListPermGens, ListSets, eGen, eList;
  EXTred:=Set(EXT);
  ListPermGens:=[];
  ListSets:=List(EXTred, x->Filtered([1..Length(EXT)], y->EXT[y]=x));
  for eGen in GeneratorsOfGroup(GroupEXT)
  do
    eList:=List(ListSets, x->Position(ListSets, OnSets(x, eGen)));
    Add(ListPermGens, PermList(eList));
  od;
  return rec(EXTred:=EXTred, GroupEXT:=Group(ListPermGens), ListSets:=ListSets);
end;



ComputeFacetInformation:=function(eCase)
  local n, CP, TheCVP, PermGRP, ListOrbFac, LVERT, eInc, EXTadj, eVert, O, LnonOrb, ListPermGens, ListMatGens, GRPmat, ListRelevantVert, RaysHYPcone, ListTangentMatrices, ExtractedBasis, ListIneq, LinearPartMatGroup, IsFinished, NewListRelevantVert, eGramMat, iSet, Lselect, ListGenerated, ListOrbRays, GRP, eSet, eOrb, TheRed, ListSizes;
  n:=Length(eCase.SuperMat);

  GRPmat:=PermGroupToMatrixGroup(eCase.EXT, eCase.GroupEXT);
  LinearPartMatGroup:=Group(List(GeneratorsOfGroup(GRPmat), x->List(x{[2..n+1]}, y->y{[2..n+1]})));
#  BasicCoherencyCheckups(rec(Basis:=eCase.Basis, TheGroup:=LinearPartMatGroup, SuperMat:=eCase.SuperMat));
  # check the Delaunay statement, a crucial thing
  CP:=CenterRadiusDelaunayPolytopeGeneral(eCase.SuperMat, eCase.EXT);
  TheCVP:=CVPVallentinProgram(eCase.SuperMat, CP.Center{[2..n+1]});
  if TheCVP.TheNorm<>CP.SquareRadius then
    Print("We find a shorter vector than it should have been\n");
    Print(NullMat(5));
  fi;
  #
  ListOrbFac:=DualDescriptionStandard(eCase.EXT, eCase.GroupEXT);
  LVERT:=[];
  for eInc in ListOrbFac
  do
    EXTadj:=FindAdjacentDelaunayPolytope(eCase.SuperMat, eCase.EXT, eInc);
    Append(LVERT, Difference(Set(EXTadj), Set(eCase.EXT)));
  od;
  LnonOrb:=Set(LVERT);
  Print("|Rel Vertices|=", Length(LnonOrb), "\n");
  ListRelevantVert:=MyOrbitGeneration(LnonOrb, GRPmat);
  ExtractedBasis:=RowReduction(eCase.EXT, n+1).EXT;
  while(true)
  do
    IsFinished:=true;
    ListIneq:=List(ListRelevantVert, x->VoronoiLinearInequality(ExtractedBasis, x, eCase.Basis));
    GRP:=MatrixGroupToPermGroup(ListRelevantVert, GRPmat);
    TheRed:=SetReductionPermGroup(ListIneq, GRP);
    ListOrbRays:=DualDescriptionStandard(TheRed.EXTred, TheRed.GroupEXT);
    RaysHYPcone:=[];
    for eOrb in ListOrbRays
    do
      for eSet in Orbit(TheRed.GroupEXT, eOrb, OnSets)
      do
        Add(RaysHYPcone, __FindFacetInequality(TheRed.EXTred, eSet));
      od;
    od;
    NewListRelevantVert:=[];
    for iSet in [1..Length(TheRed.ListSets)]
    do
      Lselect:=Filtered(RaysHYPcone, x->x*TheRed.EXTred[iSet]=0);
      if RankMat(Lselect)=Length(eCase.Basis)-1 then
        Append(NewListRelevantVert, ListRelevantVert{TheRed.ListSets[iSet]});
      fi;
    od;
    #
    ListTangentMatrices:=List(ListOrbRays, x->FuncComputeMat(eCase.Basis, __FindFacetInequality(TheRed.EXTred, x)));
    ListSizes:=List(ListOrbRays, x->Length(Orbit(TheRed.GroupEXT, x, OnSets)));
    for eGramMat in ListTangentMatrices
    do
      if RankMat(eGramMat)=7 then
        CP:=CenterRadiusDelaunayPolytopeGeneral(eGramMat, eCase.EXT);
        TheCVP:=CVPVallentinProgram(eGramMat, CP.Center{[2..n+1]});
        if CP.SquareRadius>TheCVP.TheNorm then
          ListGenerated:=List(TheCVP.ListVect, x->Concatenation([1], x));
          NewListRelevantVert:=Set(Concatenation(NewListRelevantVert, MyOrbitGeneration(ListGenerated, GRPmat)));
          IsFinished:=false;
        fi;
      fi;
    od;
    ListRelevantVert:=NewListRelevantVert;
    Print("|Relevant Vert|=", Length(ListRelevantVert), "\n");
    if IsFinished=true then
      break;
    fi;
  od;
  return rec(ListRelevantVert:=ListRelevantVert, ListTangentMatrices:=ListTangentMatrices, ListSizes:=ListSizes);
end;



SetFirstCoordTo1:=function(EXT)
  local n;
  n:=Length(EXT[1])-1;
  return List(EXT, x->Concatenation([1], x{[2..n+1]}));
end;


CreateLaminationFromExtremeDelaunay:=function(EXTwork)
  local n, DistMat, GroupEXT, EXT, V, H, GRPmat, FCR, SuperMat, i, j, eV, eCase, TheBasis, NewMat;
  n:=Length(EXTwork[1])-1;
  DistMat:=FuncDistMat(EXTwork);
  GroupEXT:=AutomorphismGroupEdgeColoredGraph(DistMat);
  EXT:=List(EXTwork, x->Concatenation(x, [0]));
  V:=ListWithIdenticalEntries(n+2, 0);
  V[1]:=1;
  V[n+2]:=1;
  Add(EXT, V);
  GRPmat:=PermGroupToMatrixGroup(EXT, GroupEXT);
  #
  FCR:=FuncExtreme(EXTwork);
  SuperMat:=NullMat(n+1,n+1);
  for i in [1..n]
  do
    for j in [1..n]
    do
      SuperMat[i][j]:=FCR.GramMatrix[i][j];
    od;
  od;
  eV:=FCR.TheCent*FCR.GramMatrix;
  for i in [1..n]
  do
    SuperMat[i][n+1]:=eV[i];
    SuperMat[n+1][i]:=eV[i];
  od;
  SuperMat[n+1][n+1]:=1+FCR.RhoSquare;
  NewMat:=NullMat(n+1,n+1);
  for i in [1..n]
  do
    for j in [1..n]
    do
      NewMat[i][j]:=FCR.GramMatrix[i][j];
    od;
  od;
  TheBasis:=[NewMat];
  for i in [1..n]
  do
    H:=NullMat(n+1,n+1);
    H[i][n+1]:=1;
    H[n+1][i]:=1;
    Add(TheBasis, H);
  od;
  H:=NullMat(n+1,n+1);
  H[n+1][n+1]:=1;
  Add(TheBasis, H);

  return rec(Basis:=TheBasis, 
             SuperMat:=RemoveFractionMatrix(SuperMat), 
             EXT:=EXT, 
             GroupEXT:=GroupEXT);
end;

CreateSchafliEnvironment:=function()
  local EXTsch;
  EXTsch:=SetFirstCoordTo1(ReadAsFunction("Functions/TheOtherSchafli")());
  return CreateLaminationFromExtremeDelaunay(EXTsch);
end;

eCase:=CreateSchafliEnvironment();
TheRet:=ComputeFacetInformation(eCase);
