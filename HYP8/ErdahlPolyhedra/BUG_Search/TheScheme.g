Read("Functions/PrintDelaunay");


MyDualDescriptionStandard:=function(EXT, PermGRP)
  local IsRespawn, IsBankSave, nbExt, WorkingDim, TmpDir, DataPolyhedral, FuncStabilizer, FuncIsomorphy, FuncInvariant, BF;
  nbExt:=Length(EXT);
  WorkingDim:=RankMat(EXT);
  #
  FuncStabilizer:=function(EXTask)
    local GRP, ListGen, eGen, ListInc;
    ListInc:=List(EXTask, x->Position(EXT, x));
    if Length(ListInc) < WorkingDim+5 then
      return Group(());
    fi;
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
    ListInc1:=List(EXT1, x->Position(EXT, x));
    ListInc2:=List(EXT2, x->Position(EXT, x));
    ePerm:=RepresentativeAction(PermGRP, ListInc1, ListInc2, OnSets);
    if ePerm=fail then
      return false;
    else
      return PermList(List(ListInc1, x->Position(ListInc2, OnPoints(x, ePerm))));
    fi;
  end;
  FuncInvariant:=function(EXT)
    return [Length(EXT), RankMat(EXT)];
  end;
  BF:=BankRecording(rec(Saving:=false, BankPath:="/irrelevant/"), FuncStabilizer, FuncIsomorphy, FuncInvariant, OnSetsGroupFormalism());
  #
  IsRespawn:=function(OrdGRP, EXT, TheDepth)
    if OrdGRP>=100 and TheDepth<=2 then
      return true;
    fi;
    if OrdGRP<100 then
      return false;
    fi;
    if Length(EXT)<WorkingDim+7 then
      return false;
    fi;
    if TheDepth=3 then
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
        GetInitialRays:=GetRAYS, 
        Saving:=false,
        ThePathSave:="/irrelevant/",
        ThePath:=TmpDir,
#        DualDescriptionFunction:=__DualDescriptionLRS_Reduction, 
#        DualDescriptionFunction:=__DualDescriptionLRS_splitter,
        DualDescriptionFunction:=__DualDescriptionCDD_Reduction,
        GroupFormalism:=OnSetsGroupFormalism());
  return __ListFacetByAdjacencyDecompositionMethod(EXT, PermGRP, DataPolyhedral, BF);
end;







FindCoordinateBasis:=function(NSP)
  local TheDim, k, eSubset, NSPselect;
  TheDim:=Length(NSP[1]);
  k:=Length(NSP);
  eSubset:=InitialSub([1..TheDim], k);
  while(true)
  do
    NSPselect:=List(NSP, x->x{eSubset});
    if AbsInt(DeterminantMat(NSPselect))=1 then
      return eSubset;
    fi;
    eSubset:=NextSubset([1..TheDim], eSubset);
    if eSubset=false then
      break;
    fi;
  od;
  Print("We should never be there\n");
  Print(NullMat(5));
end;


GetInvariantOfExtremeDelaunayPolytope:=function(EXT)
  local n, GramMatrix, ListValue, ListOcc, iVert, jVert, EXT1red, EXT2red, Vvector, TheDist, pos, DistDistrib, TheDet, SHV, TheMin, CP, eLcm, eEnt;
  n:=Length(EXT[1])-1;
  GramMatrix:=FuncExtreme(EXT).GramMatrix;
  ListValue:=[];
  ListOcc:=[];
  for iVert in [1..Length(EXT)-1]
  do
    for jVert in [iVert+1..Length(EXT)]
    do
      EXT1red:=EXT[iVert]{[2..n+1]};
      EXT2red:=EXT[jVert]{[2..n+1]};
      Vvector:=EXT1red-EXT2red;
      TheDist:=Vvector*GramMatrix*Vvector;
      pos:=Position(ListValue, TheDist);
      if pos=fail then
        Add(ListValue, TheDist);
        Add(ListOcc, 1);
      else
        ListOcc[pos]:=ListOcc[pos]+1;
      fi;
    od;
  od;
  DistDistrib:=Set(List([1..Length(ListValue)], x->[ListValue[x], ListOcc[x]]));
  TheDet:=DeterminantMat(GramMatrix);
  SHV:=ShortestVectorDutourVersion(GramMatrix);
  TheMin:=rec(min:=SHV[1]*GramMatrix*SHV[1], MinSiz:=Length(SHV));
  CP:=CenterRadiusDelaunayPolytopeGeneral(GramMatrix, EXT);
  eLcm:=1;
  for eEnt in CP.Center
  do
    eLcm:=LcmInt(eLcm, DenominatorRat(eEnt));
  od;
  return rec(DistDistrib:=DistDistrib, TheDet:=TheDet, TheMin:=TheMin, eLcm:=eLcm);
end;


FindAdjacentPerfectDelaunayPolytopes:=function(EXT, GRP)
  local n, GramMat, ListOrbFac, eInc, LOrbVert, ListValues, eLine, ListOrbHinge, ListIncidentValues, ListSelect, iAdj, eAdj, TheLift, ListAdjacent, FCR, GramMatNew, TheCentNew, TheSquareRadius, CP, TheCVP, EXTtotalNew, RPLift, EXTnew, GetValueVector, NSP, NewAffineBasis, NewBasis, eVect, eSubset, TheDiff, iCol, EXTret, nRet, EXTfiniteRet, GramMatRet, GramMatFiniteRet;
  n:=Length(EXT[1])-1;
  if Set(List(EXT, x->x[1]))<>[1] then
    Print("The first coordinate should be 1\n");
    Print(NullMat(5));
  fi;
  GramMat:=FuncExtreme(EXT).GramMatrix;
  LOrbVert:=GetVerticesAdjacentDelaunay(GramMat, EXT, GRP);
  Print("Find ", Length(LOrbVert), " possibly relevant vertices\n");
  GetValueVector:=function(eVect)
    return SymmetricMatrixToVector(TransposedMat([eVect])*[eVect]);
  end;
  Print("|Relevant vertices|=", Length(LOrbVert), "\n");
  ListIncidentValues:=List(EXT, GetValueVector);
  ListValues:=List(LOrbVert, GetValueVector);
  ListSelect:=List(ListIncidentValues, x->Position(ListValues, x));
  ListOrbHinge:=MyDualDescriptionStandard(ListIncidentValues, GRP);
  RPLift:=__ProjectionLiftingFramework(ListValues, ListSelect);
  ListAdjacent:=[];
  Print("|Cases|=", Length(ListOrbHinge), "\n");
  for iAdj in [1..Length(ListOrbHinge)]
  do
    Print("iAdj=", iAdj, "\n");
    eAdj:=ListOrbHinge[iAdj];
    TheLift:=RPLift.FuncLift(eAdj);
    EXTnew:=LOrbVert{TheLift};
    FCR:=FuncExtreme(EXTnew);
    GramMatNew:=FCR.GramMatrix;
    if RankMat(GramMatNew) = n then
      TheCentNew:=FCR.TheCent;
      TheSquareRadius:=FCR.RhoSquare;
      CP:=CenterRadiusDelaunayPolytopeGeneral(GramMatNew, EXTnew);
      if CP.SquareRadius<>FCR.RhoSquare or TheCentNew<>CP.Center{[2..n+1]} then
        Print("There is something illogic or we missed some maths\n");
        Print(NullMat(5));
      fi;
      TheCVP:=CVPVallentinProgram(GramMatNew, CP.Center{[2..n+1]});
      EXTtotalNew:=List(TheCVP.ListVect, x->Concatenation([1], x));
      Add(ListAdjacent, rec(EXTfinite:=EXTtotalNew, AddDim:=0));
    else
      NSP:=NullspaceIntMat(GramMatNew);
      NewAffineBasis:=[];
      eVect:=ListWithIdenticalEntries(Length(NSP[1])+1, 0);
      eVect[1]:=1;
      Add(NewAffineBasis, eVect);
      NewBasis:=[];
      eSubset:=FindCoordinateBasis(NSP);
      TheDiff:=Difference([1..Length(NSP[1])], eSubset);
      for iCol in TheDiff
      do
        eVect:=ListWithIdenticalEntries(Length(NSP[1])+1, 0);
        eVect[iCol+1]:=1;
        Add(NewAffineBasis, eVect);
        eVect:=ListWithIdenticalEntries(Length(NSP[1]), 0);
        eVect[iCol]:=1;
        Add(NewBasis, eVect);
      od;
      Append(NewAffineBasis, List(NSP, x->Concatenation([0], x)));
      Append(NewBasis, NSP);
      EXTret:=EXTnew*Inverse(NewAffineBasis);
      nRet:=Length(NSP[1])-Length(NSP);
      EXTfiniteRet:=List(EXTret, x->x{[1..1+nRet]});
      GramMatRet:=NewBasis*GramMatNew*TransposedMat(NewBasis);
      GramMatFiniteRet:=List(GramMatRet{[1..nRet]}, x->x{[1..nRet]});
      CP:=CenterRadiusDelaunayPolytopeGeneral(GramMatFiniteRet, EXTfiniteRet);
      TheCVP:=CVPVallentinProgram(GramMatFiniteRet, CP.Center{[2..nRet+1]});
      EXTtotalNew:=List(TheCVP.ListVect, x->Concatenation([1], x));
      Add(ListAdjacent, rec(EXTfinite:=EXTtotalNew, AddDim:=Length(NSP)));
    fi;
  od;
  return ListAdjacent;
end;


IsExtremeDelaunayPolytope:=function(EXT)
  local n, GramMat, CP, TheCVP, NewListEXT;
  n:=Length(EXT[1])-1;
  GramMat:=FuncExtreme(EXT).GramMatrix;
  CP:=CenterRadiusDelaunayPolytopeGeneral(GramMat, EXT);
  TheCVP:=CVPVallentinProgram(GramMat, CP.Center{[2..n+1]});
  NewListEXT:=List(TheCVP.ListVect, x->Concatenation([1], x));
  return Set(NewListEXT)=Set(EXT);
end;


GetAdjacentVertices:=function(EXT)
  local n, EXTrelevant, FuncInsert, i, eVert;
  n:=Length(EXT[1])-1;
  EXTrelevant:=[];
  Append(EXTrelevant, EXT);
  FuncInsert:=function(eVert)
    local pos;
    pos:=Position(EXTrelevant, eVert);
    if pos=fail then
      Add(EXTrelevant, eVert);
    fi;
  end;
  for i in [1..n]
  do
    eVert:=ListWithIdenticalEntries(n+1, 0);
    eVert[1]:=1;
    eVert[i+1]:=1;
    FuncInsert(ShallowCopy(eVert));
    eVert[i+1]:=-1;
    FuncInsert(ShallowCopy(eVert));
  od;
  return EXTrelevant;
end;

OurFunc_GetVerticesAdjacentDelaunay:=function(GramMat, EXT, GRP)
  local ListOrbFac, LVERT, eInc, EXTadj, LnonOrb, GRPmat, ListSoughtVertices, eVert, O, ListMatGens, eMatrGen, eGen, TmpDir;
  TmpDir:=Filename(POLYHEDRAL_tmpdir, "");
  ListOrbFac:=__DualDescriptionCDD_Reduction(EXT, GRP, TmpDir);
#  ListOrbFac:=DualDescriptionStandard(EXT, GRP);
  Print("We have the dual description\n");
  LVERT:=[];
  for eInc in ListOrbFac
  do
#    Print("|eInc|=", Length(eInc), "\n");
    EXTadj:=FindAdjacentDelaunayPolytope(GramMat, EXT, eInc);
    Append(LVERT, Difference(Set(EXTadj), Set(EXT)));
  od;
  Print("We have the adjacent Delaunay polytopes\n");
  LnonOrb:=Set(LVERT);
  ListMatGens:=[];
  for eGen in GeneratorsOfGroup(GRP)
  do
    eMatrGen:=__LemmaFindTransformation(EXT, EXT, eGen);
    if IsIntegralMat(eMatrGen)=false then
      Print("The matrix is not integral, sorry for that\n");
      Print(NullMat(5));
    fi;
    Add(ListMatGens, eMatrGen);
  od;
  GRPmat:=Group(ListMatGens);
  ListSoughtVertices:=[];
  while(true)
  do
    eVert:=LnonOrb[1];
    O:=Orbit(GRPmat, eVert);
    Append(ListSoughtVertices, O);
    LnonOrb:=Difference(LnonOrb, Set(O));
    if Length(LnonOrb)=0 then
      break;
    fi;
  od;
  Print("We have the complete list of relevant (so far) vertices\n");
  Append(ListSoughtVertices, EXT);
  return ListSoughtVertices;
end;



__FindAdjacentPerfectDelaunayPolytopes:=function(EXT, GRP)
  local n, GramMat, ListOrbFac, eInc, EXTrelevant, ListValues, eLine, ListOrbHinge, ListIncidentValues, ListSelect, iAdj, eAdj, TheLift, ListAdjacent, FCR, GramMatNew, TheCentNew, TheSquareRadius, CP, TheCVP, EXTtotalNew, RPLift, EXTnew, GetValueVector, NSP, NewAffineBasis, NewBasis, eVect, eSubset, TheDiff, iCol, EXTret, nRet, EXTfiniteRet, GramMatRet, GramMatFiniteRet, NewListEXT, TheDiffEXT;
  n:=Length(EXT[1])-1;
  if Set(List(EXT, x->x[1]))<>[1] then
    Print("The first coordinate should be 1\n");
    Print(NullMat(5));
  fi;
  GramMat:=FuncExtreme(EXT).GramMatrix;
  EXTrelevant:=OurFunc_GetVerticesAdjacentDelaunay(GramMat, EXT, GRP);
# NOTE: We cannot use function below, because it does not guarantees
# definiteness which is a must for us.
#  EXTrelevant:=GetAdjacentVertices(EXT);
  Print("Find ", Length(EXTrelevant), " possibly relevant vertices\n");
  GetValueVector:=function(eVect)
    return SymmetricMatrixToVector(TransposedMat([eVect])*[eVect]);
  end;
  Print("|Relevant vertices|=", Length(EXTrelevant), "\n");
  ListIncidentValues:=List(EXT, GetValueVector);
  ListValues:=List(EXTrelevant, GetValueVector);
  ListSelect:=List(ListIncidentValues, x->Position(ListValues, x));
  Print("Before call to MyDualDescriptionStandard\n");
  ListOrbHinge:=DualDescriptionStandard(ListIncidentValues, GRP);
  Print(" After call to MyDualDescriptionStandard\n");
  RPLift:=__ProjectionLiftingFramework(ListValues, ListSelect);
  ListAdjacent:=[];
  Print("|Cases|=", Length(ListOrbHinge), "\n");
  for iAdj in [1..Length(ListOrbHinge)]
  do
    Print("iAdj=", iAdj, "\n");
    eAdj:=ListOrbHinge[iAdj];
    TheLift:=RPLift.FuncLift(eAdj);
    EXTnew:=EXTrelevant{TheLift};
    GramMatNew:=FuncExtreme(EXTnew).GramMatrix;
    if RankMat(GramMatNew) = n then
      while(true)
      do
        CP:=CenterRadiusDelaunayPolytopeGeneral(GramMatNew, EXTnew);
        TheCVP:=CVPVallentinProgram(GramMatNew, CP.Center{[2..n+1]});
        NewListEXT:=List(TheCVP.ListVect, x->Concatenation([1], x));
        if TheCVP.TheNorm = CP.SquareRadius then
          if IsExtremeDelaunayPolytope(NewListEXT)=false then
            Print("We have reached a problem at this spot\n");
            Print(NullMat(5));
          fi;
          Add(ListAdjacent, NewListEXT);
          break;
        fi;
        TheDiffEXT:=Filtered(NewListEXT, x->Position(EXTrelevant, x)=fail);
        Append(EXTrelevant, TheDiffEXT);
        ListValues:=List(EXTrelevant, GetValueVector);
        RPLift:=__ProjectionLiftingFramework(ListValues, ListSelect);
        TheLift:=RPLift.FuncLift(eAdj);
        EXTnew:=EXTrelevant{TheLift};
        GramMatNew:=FuncExtreme(EXTnew).GramMatrix;
      od;
    fi;
  od;
  return ListAdjacent;
end;


SetFirstCoordTo1:=function(EXT)
  local n;
  n:=Length(EXT[1])-1;
  return List(EXT, x->Concatenation([1], x{[2..n+1]}));
end;


FuncComputeAdjacencies:=function(EXTinput)
  local EXT, DistMat, GroupEXT;
  EXT:=SetFirstCoordTo1(EXTinput);
  DistMat:=FuncDistMat(EXT);
  GroupEXT:=AutomorphismGroupEdgeColoredGraph(DistMat);
  return FindAdjacentPerfectDelaunayPolytopes(EXT, GroupEXT);
end;








FuncDoAdjacencyChecks:=function(TheInitialList)
  local ListExtremeDelaunay, FuncInsert, eRec, IsFinished, iDelaunay, DistMat, GRP, ListEXT, TheEXT, RelEXT, iVal, jDelaunay, ListEdges, ListInfinites, TheRec;
  ListExtremeDelaunay:=[];
  ListInfinites:=[];
  FuncInsert:=function(EXT)
    local DistMat, eStr, iDelaunay;
    DistMat:=FuncDistMat(EXT);
    eStr:=CanonicalStringEdgeColoredGraph(DistMat);
    for iDelaunay in [1..Length(ListExtremeDelaunay)]
    do
      if ListExtremeDelaunay[iDelaunay].eStr=eStr then
        return iDelaunay;
      fi;
    od;
    Add(ListExtremeDelaunay, rec(EXT:=EXT, eStr:=eStr, Status:="NO"));
    Print("|ListExtremeDelaunay|=", Length(ListExtremeDelaunay), "\n");
    return Length(ListExtremeDelaunay);
  end;
  for eRec in TheInitialList
  do
    iVal:=FuncInsert(SetFirstCoordTo1(eRec.EXT));
  od;
  ListEdges:=[];
  while(true)
  do
    IsFinished:=true;
    for iDelaunay in [1..Length(ListExtremeDelaunay)]
    do
      if ListExtremeDelaunay[iDelaunay].Status="NO" then
        IsFinished:=false;
        RelEXT:=ListExtremeDelaunay[iDelaunay].EXT;
        DistMat:=FuncDistMat(RelEXT);
        GRP:=AutomorphismGroupEdgeColoredGraph(DistMat);
        ListEXT:=FindAdjacentPerfectDelaunayPolytopes(RelEXT, GRP);
        for TheRec in ListEXT
        do
          if TheRec.AddDim>0 then
            Add(ListInfinites, TheRec);
          else
            TheEXT:=TheRec.EXTfinite;
            jDelaunay:=FuncInsert(TheEXT);
            Add(ListEdges, [iDelaunay, jDelaunay]);
          fi;
        od;
        ListExtremeDelaunay[iDelaunay].Status:="YES";
      fi;
    od;
    if IsFinished=true then
      break;
    fi;
  od;
  return rec(ListExtremeDelaunay:=ListExtremeDelaunay,
             ListEdges:=ListEdges, 
             ListInfinites:=ListInfinites);
end;



FuncComputationsSelection:=function(n)
  local ListExtremeDelaunay, FuncInsert, eRec, IsFinished, iDelaunay, DistMat, GRP, ListEXT, TheEXT, RelEXT, iVal, jDelaunay, ListEdges, ListInfinites, TheRec, TheSelection, MaxIncidence, EXT, nbEXT, FileName;
  ListExtremeDelaunay:=[];
  FuncInsert:=function(EXT)
    local DistMat, DistMatCall, eInv, iDelaunay, nbDelaunay, TheOrb, test;
    DistMat:=FuncDistMat(EXT);
    eInv:=GetInvariantOfExtremeDelaunayPolytope(EXT);
    for iDelaunay in [1..Length(ListExtremeDelaunay)]
    do
      if ListExtremeDelaunay[iDelaunay].eInv=eInv then
        DistMatCall:=FuncDistMat(ListExtremeDelaunay[iDelaunay].EXT);
        test:=IsIsomorphicEdgeColoredGraph(DistMat, DistMatCall);
        if test<>false then
          return;
        fi;
      fi;
    od;
    nbDelaunay:=Length(ListExtremeDelaunay)+1;
    TheOrb:=rec(EXT:=EXT, eInv:=eInv, Status:="NO");
    Add(ListExtremeDelaunay, TheOrb);
    FileName:=Concatenation("./SAVING/DelaunayEXT", String(nbDelaunay));
    SaveDataToFilePlusTouch(FileName, TheOrb);
    Print("|ListExtremeDelaunay|=", Length(ListExtremeDelaunay), "\n");
  end;
  nbEXT:=0;
  while(true)
  do
    FileName:=Concatenation("./SAVING/DelaunayEXT", String(nbEXT+1));
    if IsExistingFilePlusTouch(FileName)=false then
      break;
    fi;
    Add(ListExtremeDelaunay, ReadAsFunction(FileName)());
    nbEXT:=nbEXT+1;
  od;
  if nbEXT=0 then
    EXT:=ErdahlInfiniteSequence(n);
    FuncInsert(EXT);
  fi;
  while(true)
  do
    IsFinished:=true;
    TheSelection:=0;
    MaxIncidence:=2^n;
    for iDelaunay in [1..Length(ListExtremeDelaunay)]
    do
      if ListExtremeDelaunay[iDelaunay].Status="NO" then
        IsFinished:=false;
        if Length(ListExtremeDelaunay[iDelaunay].EXT) < MaxIncidence then
          MaxIncidence:=Length(ListExtremeDelaunay[iDelaunay].EXT);
          TheSelection:=iDelaunay;
        fi;
      fi;
    od;
    if IsFinished=false then
      RelEXT:=ListExtremeDelaunay[TheSelection].EXT;
      DistMat:=FuncDistMat(RelEXT);
      GRP:=AutomorphismGroupEdgeColoredGraph(DistMat);
      ListEXT:=__FindAdjacentPerfectDelaunayPolytopes(RelEXT, GRP);
      for TheEXT in ListEXT
      do
        FuncInsert(TheEXT);
      od;
      ListExtremeDelaunay[TheSelection].Status:="YES";
      FileName:=Concatenation("./SAVING/DelaunayEXT", String(TheSelection));
      SaveDataToFile(FileName, ListExtremeDelaunay[TheSelection]);
    fi;
    if IsFinished=true then
      break;
    fi;
  od;
  return ListExtremeDelaunay;
end;




FuncDoAdjacencyChecksDim8:=function()
  local LE;
  LE:=ReadAsFunction("Data/ExtremeDim8")();
  return FuncDoAdjacencyChecks(LE);
end;
FuncDoAdjacencyChecksDim7:=function()
  local LE, LE7;
  LE:=ReadAsFunction("Data/LE78")();
  LE7:=Filtered(LE, x->Length(x.EXT[1])=8);
  return FuncDoAdjacencyChecks(LE7);
end;



MultiplyBy2diagonal:=function(SymMat)
  local n, RetMat, i, j;
  n:=Length(SymMat);
  RetMat:=NullMat(n,n);
  for i in [1..n]
  do
    for j in [1..n]
    do
      if i<>j then
        RetMat[i][j]:=SymMat[i][j];
      else
        RetMat[i][j]:=2*SymMat[i][j];
      fi;
    od;
  od;
  return RetMat;
end;

LowerQuadrangle:=function(Zmat)
  local n;
  n:=Length(Zmat);
  return List(Zmat{[2..n]}, x->x{[2..n]});
end;

ZeroFunctionExtremeDelaunay:=function(EXTfinite)
  local n, ListVect, eVect, NSP, TheMat, TheMatRed;
  n:=Length(EXTfinite[1])-1;
  ListVect:=[];
  for eVect in EXTfinite
  do
    Add(ListVect, SymmetricMatrixToVector(TransposedMat([eVect])*[eVect]));
  od;
  NSP:=NullspaceMat(TransposedMat(ListVect));
  if Length(NSP)>1 then
    Print("We have a problem somewhere\n");
    Print(NullMat(5));
  fi;
  TheMat:=VectorToSymmetricMatrix(NSP[1], n+1);
  TheMatRed:=MultiplyBy2diagonal(List(TheMat{[2..n+1]}, x->x{[2..n+1]}));
  if IsPositiveDefiniteSymmetricMatrix(TheMatRed)=false then
    return -TheMat;
  else
    return TheMat;
  fi;
end;

AddDimensions:=function(SymMat, kAdd)
  local n, RetMat, i, j;
  n:=Length(SymMat);
  RetMat:=NullMat(n+kAdd, n+kAdd);
  for i in [1..n]
  do
    for j in [1..n]
    do
      RetMat[i][j]:=SymMat[i][j];
    od;
  od;
  return RetMat;
end;


SingleInitial:=function(EXTfinite)
  local n, TheInitial, eEXT, GramMat, Zmat;
  n:=Length(EXTfinite[1]);
  TheInitial:=[];
  for eEXT in EXTfinite
  do
    Add(TheInitial, Concatenation(eEXT, [0]));
    Add(TheInitial, Concatenation(eEXT, [1]));
  od;
  Zmat:=NullMat(n+1, n+1);
  Zmat[n+1][n+1]:=1;
  return rec(EXTfinite:=TheInitial, Zmat:=Zmat);
end;

SingleInfiniteFromExtremeDelaunay:=function(EXTfinite)
  local Zmat, ZmatAdd, DM, GRP;
  Zmat:=ZeroFunctionExtremeDelaunay(EXTfinite);
  ZmatAdd:=AddDimensions(Zmat, 1);
  DM:=DistanceMatrixExtremeDelaunay(EXTfinite);
  GRP:=AutomorphismGroupEdgeColoredGraph(DM);
  return rec(EXT:=rec(EXTfinite:=EXTfinite, AddDim:=1),
             ListZmat:=[ZmatAdd],
             DM:=DM,
             GRP:=GRP);
end;

IsomorphyTools:=function(EXTfinite, k)
  local n, EXTret, eEXT;
  n:=Length(EXTfinite[1]);
  EXTret:=[];
  for eEXT in EXTfinite
  do
    Add(EXTret, eEXT{[1..n-k]});
  od;
  return EXTret;
end;
IsobarycenterSlow:=function(TheList)
  local TheRet, iObj;
  TheRet:=TheList[1];
  for iObj in [2..Length(TheList)]
  do
    TheRet:=TheRet+TheList[iObj];
  od;
  return TheRet/Length(TheList);
end;



CharacteristicDistMatSingleInfinite:=function(EXTfinite, EXTsingleInfinity)
  local n, ScalarMat, DistMat1, Zmat, GramMat, nbVert, nbVertTot, DistMat2, iVert, jVert, eDiff, eDist, DistMatRet, eList;
  n:=Length(EXTfinite[1])-1;
  ScalarMat:=__VectorConfigurationFullDim_ScalarMat(EXTfinite);
  DistMat1:=MappedScalarMatrixDistanceMatrix(ScalarMat);
  Zmat:=IsobarycenterSlow(EXTsingleInfinity.ListZmat);
  GramMat:=MultiplyBy2diagonal(List(Zmat{[2..n+1]}, x->x{[2..n+1]}));
  nbVert:=Length(EXTfinite);
  nbVertTot:=Length(EXTfinite)+2;
  DistMat2:=NullMat(nbVertTot, nbVertTot);
  for iVert in [1..nbVert-1]
  do
    for jVert in [iVert+1..nbVert]
    do
      eDiff:=EXTfinite[iVert]{[2..n+1]}-EXTfinite[iVert]{[2..n+1]};
      eDist:=eDiff*GramMat*eDiff;
      DistMat2[iVert][jVert]:=eDist;
    od;
  od;
  DistMatRet:=[];
  for iVert in [1..nbVertTot]
  do
    eList:=[];
    for jVert in [1..nbVertTot]
    do
      Add(eList, [DistMat1[iVert][jVert], DistMat2[iVert][jVert]]);
    od;
    Add(DistMatRet, eList);
  od;
  return DistMatRet;
end;


AutomorphismFiniteStructure:=function(EXTfinite, EXTsingleInfinite)
  local DistMat, GRP, eGen, eMat, EXTsingleInfinity;
  DistMat:=CharacteristicDistMatSingleInfinite(EXTfinite, EXTsingleInfinity);
  GRP:=AutomorphismGroupEdgeColoredGraph(DistMat);
  for eGen in GeneratorsOfGroup(GRP)
  do
    eMat:=__LemmaFindTransformation(EXTfinite, EXTfinite, eGen);
    if IsIntegralMat(eMat)=false then
      Print("Let's be fair, we are in shit now 1\n");
      Print(NullMat(5));
    fi;
  od;
  return GRP;
end;

EquivalenceFiniteStructure:=function(EXT1finite, EXT2finite, EXTsingleInfinite)
  local DistMat1, DistMat2, test, ePerm, eMat, EXTfinite, EXTsingleInfinity;
  DistMat1:=CharacteristicDistMatSingleInfinite(EXTfinite, EXTsingleInfinity);
  test:=IsIsomorphicEdgeColoredGraph(DistMat1, DistMat2);
  if test=false then
    return test;
  fi;
  ePerm:=PermList(test);
  eMat:=__LemmaFindTransformation(EXT1finite, EXT2finite, ePerm);
  if IsIntegralMat(eMat)=false then
    Print("Let's be fair, we are in shit now 2\n");
    Print(NullMat(5));
  fi;
  return test;
end;

GetValueVector:=function(eVect)
  return SymmetricMatrixToVector(TransposedMat([eVect])*[eVect]);
end;


IntersectionFiniteInfinite:=function(EXTfinite, EXTinfinite)
  local EXTret, eVert, eVertRed, InfDim;
  InfDim:=Length(EXTinfinite.EXTfinite[1]);
  EXTret:=[];
  for eVert in EXTfinite
  do
    eVertRed:=eVert{[1..InfDim]};
    if Position(EXTinfinite.EXTfinite, eVertRed)<>fail then
      Add(EXTret, eVert);
    fi;
  od;
  return EXTret;
end;


AdjacencyScheme_V1:=function(EXTsingleInfinite)
  local EXTfinite, TheInitial, ListCases, FuncInsert, nbCase, IsFinished, iCase, Zmat, GRP, ZmatRel, LVERT1, LVERT2, ListIncidentValues, ListValues2, ListSelect2, ListOrbHinge, RPLift2, iAdj, eAdj2, TheLift2, EXTnew, GramMat, EXT, LOrbVert, StructEXTfinite;
  EXTfinite:=EXTsingleInfinite.EXT.EXTfinite;
  TheInitial:=SingleInitial(EXTfinite);
  ListCases:=[TheInitial];
  FuncInsert:=function(EXTfinite, Zmat)
    local eCase, test;
    for eCase in ListCases
    do
      test:=EquivalenceFiniteStructure(eCase.EXTfinite, StructEXTfinite.EXTfinite, EXTsingleInfinite);
      if test<>false then
        return;
      fi;
    od;
    Add(ListCases, rec(EXTfinite:=EXTfinite, Zmat:=Zmat, Status:="NO"));
  end;
  while(true)
  do
    nbCase:=Length(ListCases);
    IsFinished:=true;
    for iCase in [1..nbCase]
    do
      if ListCases[iCase].Status="NO" then
        ListCases[iCase].Status:="YES";
        EXTfinite:=ListCases[iCase].EXTfinite;
        Zmat:=ListCases[iCase].Zmat;
        IsFinished:=false;
        GRP:=AutomorphismFiniteStructure(EXTfinite, EXTsingleInfinite);
        ZmatRel:=(EXTsingleInfinite.ListZmat[1]+Zmat)/2;
        GramMat:=MultiplyBy2diagonal(LowerQuadrangle(ZmatRel));
        LVERT1:=GetVerticesAdjacentDelaunay(GramMat, EXT, GRP);
        Print("Find ", Length(LOrbVert), " possibly relevant vertices\n");
        LVERT2:=IntersectionFiniteInfinite(LVERT1, EXTsingleInfinite);
        ListIncidentValues:=List(EXT, GetValueVector);
        ListValues2:=List(LVERT2, GetValueVector);
        ListSelect2:=List(ListIncidentValues, x->Position(ListValues2, x));
        ListOrbHinge:=DualDescriptionStandard(ListIncidentValues, GRP);
        RPLift2:=__ProjectionLiftingFramework(ListValues2, ListSelect2);
#wrong        ListValues1:=List(LVERT1, GetValueVector);
#wrong        ListSelect1:=List(ListIncidentValues, x->Position(ListValues1, x));
#wrong        RPLift1:=__ProjectionLiftingFramework(ListValues1, ListSelect1);
        Print("|Cases|=", Length(ListOrbHinge), "\n");
        for iAdj in [1..Length(ListOrbHinge)]
        do
          Print("iAdj=", iAdj, "\n");
          eAdj2:=ListOrbHinge[iAdj];
          TheLift2:=RPLift2.FuncLift(eAdj2);
          EXTnew:=LVERT2{TheLift2};
          # there is a critical missing spot here with respect to
          # the forms. If S is a rank 2 subset of Sch x Z, how do we
          # find the other forms containing it?
          # We seem to need a RPLlift for infinite structures.
        od;
      fi;
    od;
    if IsFinished=true then
      break;
    fi;
  od;
  return ListCases;
end;
