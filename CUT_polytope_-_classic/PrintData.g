TreatFile_MET:=function(eFile)
  local fFile, fRec, eRec, GRA, FileSave, K5gra, IsMinor, RecSave, eRecB, NameRed, FAC, GRP, FuncStabilizer, FuncIsomorphy, FuncInvariant, IsRespawn, TestNeedMoreSymmetry, IsBankSave, TheFunc, Data, DataBank, BF, LORB, BankPath, ThePathWork, PathSave, ThePath, LimitNbVert;
  fFile:=Concatenation("CutPolytopes/", eFile);
  if IsExistingFile(fFile)=false then
    Print("Missing fFile=", fFile, "\n");
    Error("Please correct");
  fi;
  fRec:=ReadAsFunction(fFile)();
  eRec:=fRec.eRec;
  GRA:=eRec.GRA;
  FileSave:=Concatenation("MetPolytopes/", eFile);
  if IsExistingFile(FileSave)=false then
    K5gra:=CompleteGraph(Group(()), 5);
    IsMinor:=SAGE_TestIfHasMinor(GRA, K5gra);
    RecSave:=rec(IsMinor:=IsMinor);
    if IsMinor then
      eRecB:=CMC_GetCutPolytope(GRA);
      NameRed:=eFile{[5..Length(eFile)]};
      FAC:=eRecB.MET_info.EXT;
      GRP:=LinPolytope_Automorphism(FAC);
      #
      ThePathWork:=Concatenation("SCRATCH_MET/SAVE_", NameRed, "/");
      Exec("mkdir -p ", ThePathWork);
      ThePath:=Concatenation(ThePathWork, "tmp/");
      Exec("mkdir -p ", ThePath);
      PathSave:=Concatenation(ThePathWork, "PathSAVE/");
      Exec("mkdir -p ", PathSave);
      #
      FuncStabilizer:=LinPolytope_Automorphism;
      FuncIsomorphy:=LinPolytope_Isomorphism;
      FuncInvariant:=LinPolytope_Invariant;
      #
      IsRespawn:=function(OrdGRP, EXT, TheDepth)
        local nbInc;
        nbInc:=Length(EXT);
        if nbInc <= 70 then
          return false;
        fi;
        return true;
      end;
      #
      TestNeedMoreSymmetry:=function(EXT)
        if Length(EXT)<=50 then
          return false;
        else
          return true;
        fi;
      end;
      #
      IsBankSave:=function(EllapsedTime, OrdGRP, EXT, TheDepth)
        if TheDepth=0 then
          return false;
        fi;
        if EllapsedTime>=600 then
          return true;
        fi;
        if Length(EXT)<=50 then
          return false;
        fi;
        return true;
      end;
      #
#      TheFunc:=__DualDescriptionCDD_Reduction;
      TheFunc:=__DualDescriptionPPL_Reduction;
      #
      LimitNbVert:=100;
      Data:=rec(TheDepth:=0,
            ThePath:=ThePath,
            GetInitialRays:=GetInitialRays_LinProg, 
            IsBankSave:=IsBankSave,
            GroupFormalism:=OnSetsGroupFormalism(LimitNbVert), 
            TestNeedMoreSymmetry:=TestNeedMoreSymmetry,
            DualDescriptionFunction:=TheFunc, 
            IsRespawn:=IsRespawn,
            Saving:=true,
            ThePathSave:=PathSave);
      #
      BankPath:=Concatenation(ThePathWork, "TheBank/");
      Exec("mkdir -p ", BankPath);
      DataBank:=rec(BankPath:=BankPath, Saving:=true);
      BF:=BankRecording(DataBank, FuncStabilizer, FuncIsomorphy, FuncInvariant, OnSetsGroupFormalism(LimitNbVert));
      #
      LORB:=__ListFacetByAdjacencyDecompositionMethod(FAC, GRP, Data, BF);
      RecSave.LORB:=LORB;
      RecSave.FAC:=FAC;
      RecSave.GRP:=GRP;
    fi;
    SaveDataToFile(FileSave, RecSave);
  fi;
end;


ListClassicalData:=[];
FuncInsertGraphCase:=function(eFile)
  local fFile, fRec, eRec, GRA, RecGraph, NameRed;
  fFile:=Concatenation("CutPolytopes/", eFile);
  if IsExistingFile(fFile)=false then
    Print("Missing fFile=", fFile, "\n");
    Error("Please correct");
  fi;
  fRec:=ReadAsFunction(fFile)();
  eRec:=fRec.eRec;
  GRA:=eRec.GRA;
  NameRed:=eFile{[5..Length(eFile)]};
  RecGraph:=rec(name:=NameRed, GRA:=GRA);
  Add(ListClassicalData, RecGraph);
end;



TreatFile:=function(eFile)
  local fRec, eRec, eFileB, GRA, n, ListEdges, i, LAdj, j, eEdge, nbEdge, ListPermSwitch, eMat, iEdge, nbOrb, GRPgra, eStab, eSize, FinalFileTex, FinalFileDVI, FinalFileTexB, FinalFileDVIB, eList, ePerm, output, eEnt, ListEnt, eAdj, iOrb, FinalFilePSB, SizeGRPgraSwitch, NewListInc, GRPcanonic, ListPermGraph, eGen, nEdge, ePermEdge, eListB, eVect, eVectRed, fVectRed, fVect, ePermGraph, FunctionPrint, fOrb, fFile, HypermetricRecognition, GetStringHypermetric, EdgeIneqRecognition, ScycleRecognition, RecognizeInequality;
  fFile:=Concatenation("CutPolytopes/", eFile);
  if IsExistingFile(fFile)=false then
    Print("Missing fFile=", fFile, "\n");
    Error("Please correct");
  fi;
  fRec:=ReadAsFunction(fFile)();
  eRec:=fRec.eRec;
  eFileB:=Concatenation("TEX/", eFile, "_pre.tex");
  RemoveFileIfExist(eFileB);
  GRA:=eRec.GRA;
  n:=GRA.order;
  ListEdges:=[];
  for i in [1..n-1]
  do
    LAdj:=Adjacency(GRA, i);
    for j in LAdj
    do
      if j>i then
        eEdge:=[i, j];
        Add(ListEdges, eEdge);
      fi;
    od;
  od;
  nbEdge:=Length(ListEdges);
  #
  EdgeIneqRecognition:=function(eIneq)
    local NNZ;
    NNZ:=Filtered([1..nbEdge], x->eIneq[1+x]<>0);
    if Length(NNZ)=1 then
      return Concatenation("edge inequality e=", String(ListEdges[NNZ[1]]));
    fi;
    return false;
  end;
  ScycleRecognition:=function(eIneq)
    local ListPairP, ListPairM, TheList, iEdge, ePt, i, nbPairP, pos, posB, TheCycle, ListStatus, iter, eLast, iPairSel, iPair;
    ListPairP:=[];
    ListPairM:=[];
    TheList:=ListWithIdenticalEntries(n,0);
    for iEdge in [1..nbEdge]
    do
      if AbsInt(eIneq[1+iEdge]) > 1 then
        return false;
      fi;
      if eIneq[1+iEdge]=1 then
	Add(ListPairP, ListEdges[iEdge]);
      fi;
      if eIneq[1+iEdge]=-1 then
	Add(ListPairM, ListEdges[iEdge]);
      fi;
      if AbsInt(eIneq[1+iEdge])=1 then
        for ePt in ListEdges[iEdge]
        do
          TheList[ePt]:=TheList[ePt]+1;
        od;
      fi;
    od;
    for i in [1..n]
    do
      if Position([0,2], TheList[i])=fail then
        return false;
      fi;
    od;
    if Length(ListPairM)<>1 or Length(ListPairP) < 2 then
      return false;
    fi;
    nbPairP:=Length(ListPairP);
    ListStatus:=ListWithIdenticalEntries(nbPairP,0);
    TheCycle:=[ListPairM[1][1], ListPairM[1][2]];
    for iter in [1..nbPairP]
    do
      eLast:=TheCycle[iter+1];
      iPairSel:=-1;
      for iPair in [1..nbPairP]
      do
        if ListStatus[iPair]=0 then
          if Position(ListPairP[iPair], eLast)<>fail then
            if iPairSel<>-1 then
              return false;
            fi;
            iPairSel:=iPair;
          fi;
        fi;
      od;
      #
      pos:=Position(ListPairP[iPairSel], eLast);
      posB:=3-pos;
      ListStatus[iPairSel]:=1;
      #
      ePt:=ListPairP[iPairSel][posB];
      if iter<nbPairP then
        if Position(TheCycle, ePt)<>fail then
          return false;
        fi;
        Add(TheCycle, ePt);
      else
        if TheCycle[1]<>ePt then
          return false;
        fi;
      fi;
    od;
    if Position(ListStatus, 0)<>fail then
      return false;
    fi;
    return Concatenation(String(Length(TheCycle)), "-cycle inequality, C=", String(TheCycle), " F=", String(ListPairM[1]));
  end;
  HypermetricRecognition:=function(eIneq)
    local eIneqCorr, TheList, i, j, eFirst, iIdx, jIdx, ePair, pos, nbUndone, nbOper, eVal, TheSum, TheSign;
    eIneqCorr:=-eIneq{[2..nbEdge+1]};
    TheList:=[];
    for i in [1..n]
    do
      Add(TheList, fail);
    od;
    eFirst:=First([1..nbEdge], x->AbsInt(eIneqCorr[x])=1);
    if eFirst=fail then
      return false;
    fi;
    iIdx:=ListEdges[eFirst][1];
    jIdx:=ListEdges[eFirst][2];
    if eIneqCorr[eFirst]=1 then
      TheList[iIdx]:=1;
      TheList[jIdx]:=1;
    else
      TheList[iIdx]:=1;
      TheList[jIdx]:=-1;
    fi;
    while(true)
    do
      # Check Coherency
      for i in [1..n]
      do
        for j in [i+1..n]
        do
          ePair:=[i,j];
          pos:=Position(ListEdges, ePair);
          if pos<>fail then
            if TheList[i]<>fail and TheList[j]<>fail then
              eVal:=TheList[i]*TheList[j];
              if eVal<>eIneqCorr[pos] then
                return false;
              fi;
            fi;
          fi;
        od;
      od;
      # Counting number of undone. If zero exiting
      nbUndone:=0;
      for i in [1..n]
      do
        if TheList[i]=fail then
          nbUndone:=nbUndone+1;
        fi;
      od;
      if nbUndone=0 then
        TheSum:=Sum(TheList);
        if TheSum<0 then
          TheSign:=-1;
        else
          TheSign:=1;
        fi;
        if TheSum=0 then
          return Concatenation("Negative type inequality, b=", String(TheList * TheSign));
        else
          return Concatenation("Hypermetric, b=", String(TheList * TheSign));
        fi;
      fi;
      # Iterative completion
      nbOper:=0;
      for i in [1..n]
      do
        for j in [i+1..n]
        do
          ePair:=[i,j];
          pos:=Position(ListEdges, ePair);
          if pos<>fail then
            if TheList[i]<>fail and TheList[j]=fail then
              if TheList[i]<>0 then
                eVal:=eIneqCorr[pos] / TheList[i];
                TheList[j]:=eVal;
                nbOper:=nbOper+1;
              fi;
            fi;
            if TheList[i]=fail and TheList[j]<>fail then
              if TheList[j]<>0 then
                eVal:=eIneqCorr[pos] / TheList[j];
                TheList[i]:=eVal;
                nbOper:=nbOper+1;
              fi;
            fi;
          fi;
        od;
      od;
      if nbOper=0 then
        Print("We have not been able to make progress\n");
#        Error("Please debug");
        return false;
      fi;
    od;
  end;
  RecognizeInequality:=function(eIneq)
    local eRec1, eRec2, eRec3;
    eRec1:=EdgeIneqRecognition(eIneq);
    if eRec1<>false then
      return eRec1;
    fi;
    eRec2:=ScycleRecognition(eIneq);
    if eRec2<>false then
      return eRec2;
    fi;
    eRec3:=HypermetricRecognition(eIneq);
    if eRec3<>false then
      return eRec3;
    fi;
    return "unknown";
  end;
  #
  ListPermSwitch:=[];
  for i in [1..n]
  do
    eMat:=NullMat(nbEdge+1,nbEdge+1);
    eMat[1][1]:=1;
    for iEdge in [1..nbEdge]
    do
      eEdge:=ListEdges[iEdge];
      if Position(eEdge, i)<>fail then
        eMat[1][iEdge+1]:=1;
        eMat[iEdge+1][iEdge+1]:=-1;
      else
        eMat[iEdge+1][iEdge+1]:=1;
      fi;
    od;
    eList:=List(eRec.ListVect, x->Position(eRec.ListVect, x*eMat));
    ePerm:=PermList(eList);
    if ePerm=fail then
      Error("Consistency error in ePerm");
    fi;
    if ePerm=() then
      Error("Should not be identity");
    fi;
    Add(ListPermSwitch, ePerm);
  od;
  ListPermGraph:=[];
  for eGen in GeneratorsOfGroup(AutGroupGraph(GRA))
  do
    eList:=[];
    for eEdge in ListEdges
    do
      nEdge:=OnSets(eEdge, eGen);
      Add(eList, Position(ListEdges, nEdge));
    od;
    ePermEdge:=PermList(eList);
    eListB:=[];
    for eVect in eRec.ListVect
    do
      eVectRed:=eVect{[2..nbEdge+1]};
      fVectRed:=Permuted(eVectRed, ePermEdge);
      fVect:=Concatenation([1], fVectRed);
      Add(eListB, Position(eRec.ListVect, fVect));
    od;
    ePermGraph:=PermList(eListB);
    Add(ListPermGraph, ePermGraph);
  od;
  GRPcanonic:=Group(Concatenation(ListPermGraph, ListPermSwitch));
  #
  output:=OutputTextFile(eFileB, true);
  AppendTo(output, "Number of vertices \$n=", n, "\$.\n\n");
  AppendTo(output, "Adjacencies of Graph\n");
  AppendTo(output, "\\begin{enumerate}\n");
  for i in [1..n]
  do
    AppendTo(output, "\\item vertex ", i, " adjacent to");
    for eAdj in Adjacency(GRA, i)
    do
      AppendTo(output, " ", eAdj);
    od;
    AppendTo(output, "\n");
  od;
  AppendTo(output, "\\end{enumerate}\n");
  nbOrb:=Length(fRec.LOrb);
  GRPgra:=AutGroupGraph(GRA);
  SizeGRPgraSwitch:=Order(GRPgra)*2^(n-1);
  #
  #
  AppendTo(output, "Size of automorphism group of the graph=", Order(GRPgra), "\n\n");
  AppendTo(output, "Full group: \$\\vert Aut(polytope) \\vert  =", Order(fRec.GRP), "\$\n\n");
  AppendTo(output, "Restricted group: \$\\vert Aut(G) \\times switch \\vert = ", SizeGRPgraSwitch, "\$\n\n");
  AppendTo(output, "Number of orbits for the full group : ", nbOrb, "\n\n");
  AppendTo(output, "List of orbits of facets for the full group:\n");
  FunctionPrint:=function(PreListOrb, GRPrel)
    local ListLen, ePerm, ListOrb, nbOrb, eOrb, eStab, eSize, iSet, zPerm, eVal, nOrb, eIneq, iEdge, eEdge, eEnt, nbCol, NbSymPerItem, TotalNumberFacet;
    ListLen:=List(PreListOrb, x->1/(Length(x)));
    ePerm:=SortingPerm(ListLen);
    ListOrb:=Permuted(PreListOrb, ePerm);
    nbOrb:=Length(ListOrb);
    AppendTo(output, "Total number of orbits = ", nbOrb, "\n");
    TotalNumberFacet:=0;
    for iOrb in [1..nbOrb]
    do
      eOrb:=ListOrb[iOrb];
      eStab:=Stabilizer(GRPrel, eOrb, OnSets);
      eSize:=Order(GRPrel)/Order(eStab);
      TotalNumberFacet:=TotalNumberFacet + eSize;
    od;
    AppendTo(output, "Total number of facets = ", TotalNumberFacet, "\n");
    AppendTo(output, "\\begin{enumerate}\n");
    for iOrb in [1..nbOrb]
    do
      eOrb:=ListOrb[iOrb];
      eStab:=Stabilizer(GRPrel, eOrb, OnSets);
      eSize:=Order(GRPrel)/Order(eStab);
      iSet:=eOrb[1];
      zPerm:=();
      for eVal in eRec.ListSet[iSet]
      do
        zPerm:=zPerm*ListPermSwitch[eVal];
      od;
      nOrb:=OnSets(eOrb, zPerm);
      eIneq:=__FindFacetInequality(fRec.EXT, nOrb);
      if eIneq[1]<>0 then
        Error("BUG clearly\n");
      fi;
      AppendTo(output, "\\item Inequality ", iOrb, " with incidence ", Length(eOrb), " and stabilizer of size ", Order(eStab), ". Orbit size is ", eSize, " nature: ", RecognizeInequality(eIneq), "\n");
      ListEnt:=[];
      for iEdge in [1..nbEdge]
      do
        eEdge:=ListEdges[iEdge];
        eVal:=eIneq[1+iEdge];
        eEnt:=Concatenation("(", String(eEdge[1]), ",", String(eEdge[2]), ") : ", String(eVal));
        Add(ListEnt, eEnt);
      od;
      nbCol:=6;
      NbSymPerItem:=1;
      AppendTo(output, "\\begin{center}\n");
      LatexPrintInTab(output, ListEnt, nbCol, NbSymPerItem);
      AppendTo(output, "\\end{center}\n");
    od;
    AppendTo(output, "\\end{enumerate}\n");
  end;
  FunctionPrint(fRec.LOrb, fRec.GRP);
  if Order(fRec.GRP)<>Order(GRPcanonic) then
    NewListInc:=[];
    for fOrb in fRec.LOrb
    do
      Append(NewListInc, __IndividualLifting(fOrb, GRPcanonic, fRec.GRP));
    od;
    AppendTo(output, "Number of orbits for the restricted group : ", Length(NewListInc), "\n\n");
    AppendTo(output, "List of orbits of facets for the restricted group:\n");
    FunctionPrint(NewListInc, GRPcanonic);
  fi;
  CloseStream(output);
  #
  FinalFileTex:=Concatenation("TEX/", eFile, ".tex");
  FinalFileDVI:=Concatenation("TEX/", eFile, ".dvi");
  FinalFileTexB:=Concatenation(eFile, ".tex");
  FinalFileDVIB:=Concatenation(eFile, ".dvi");
  FinalFilePSB:=Concatenation(eFile, ".ps");
  Exec("cat TEX/Header.tex ", eFileB, " TEX/Footer.tex > ", FinalFileTex);
  Exec("(cd TEX && latex ", FinalFileTexB, ")");
  Exec("(cd TEX && dvips ", FinalFileDVIB, " -o )");
  Exec("(cd TEX && ps2pdf ", FinalFilePSB, ")");
end;

ListFile1:=["DataK22", "DataK23", "DataK24", "DataK25", "DataK26", "DataK27", "DataK28"];
ListFile2:=["DataK33", "DataK34", "DataK35", "DataK36"];
ListFile3:=["DataK44", "DataK45", "DataK46", "DataK47"];
ListFile4:=["DataK122", "DataK123", "DataK124", "DataK125", "DataK126", "DataK133", "DataK134", "DataK135", "DataK136", "DataK144", "DataK333"];

ListFile5:=["DataCube", "DataOctahedron", "DataHeawood", "DataPetersen", "DataIcosahedron"];
ListFile6:=["DataPrism3", "DataPrism4", "DataPrism5", "DataPrism6", "DataPrism7"];
ListFile7:=["DataAntiPrism4", "DataAntiPrism5", "DataAntiPrism6"];
ListFile8:=["DataTruncatedTetrahedron", "DataCuboctahedron"];
ListFile9:=["DataMobiusLadder4", "DataMobiusLadder5", "DataMobiusLadder6", "DataMobiusLadder7"];
ListFile10:=["DataK113", "DataK114", "DataK1112", "DataK1113"];
ListFile11:=["DataK1122", "DataK1123", "DataK1124", "DataK1125", "DataK1133", "DataK7-K2", "DataK6-K2", "DataK7-K3", "DataK8-K3", "DataK8-K4", "DataK9-K5", "DataK7-C4", "DataK7-C7", "DataK7-K4", "DataK8-K5", "DataK9-K6"];
ListFile12:=["DataProjGC11cube"];
ListFile13:=["DataPyrPrism3", "DataPyrPrism4", "DataPyrOctahedron", "DataPyrPyrCycl4", "DataPyrPyrCycl5", "DataPyrPyrCycl6", "DataPyrPyrCycl7", "DataPyrAntiPrism4", "DataPyrOctahedron", "DataPyrPrism5"];
ListFile14:=["DataK112", "DataK115", "DataK13", "DataK14", "DataK15"];
ListFile:=Concatenation(ListFile1, ListFile2, ListFile3, ListFile4, ListFile5,
      ListFile6, ListFile7, ListFile8, ListFile9, ListFile10, ListFile11, ListFile11,
      ListFile12, ListFile13, ListFile14);
#  ListFile:=ListFile11;
#  ListFile:=ListFile4;
#  ListFile:=ListFile5;
#ListFile:=["DataK7-K2"];
#ListFile:=["DataAntiPrism4"];
#ListFile:=["DataK55"];
ListFile:=["DataK11122"];

for eFile in ListFile
do
  TreatFile(eFile);
#  TreatFile_MET(eFile);
#  FuncInsertGraphCase(eFile);
od;

