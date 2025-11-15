GetDualDescInfo:=function(ListSize)
  local GRA, eRec, EXT, GRP, Lorb, nbFacet, eOrb, eStab, eSize, eVal, eFile, fRec;
  GRA:=GRAPH_GetMultiComplement(ListSize);
  eRec:=CMC_GetCutPolytope(GRA);
  EXT:=eRec.ListVect;
  GRP:=LinPolytope_Automorphism(EXT);
  Lorb:=DualDescriptionStandard(EXT, GRP);
  nbFacet:=0;
  for eOrb in Lorb
  do
    eStab:=Stabilizer(GRP, eOrb, OnSets);
    eSize:=Order(GRP)/Order(eStab);
    nbFacet:=nbFacet + eSize;
  od;
  fRec:=rec(eRec:=eRec, 
            EXT:=eRec.ListVect,
            LOrb:=Lorb, 
            GRP:=GRP);
  eFile:="DataK";
  for eVal in ListSize
  do
    eFile:=Concatenation(eFile, String(eVal));
  od;
  SaveDataToFile(eFile, fRec);
  return fRec;
end;
