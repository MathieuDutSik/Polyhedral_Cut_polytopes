FincTheN:=function(nb)
  local TheN, nb1;
  TheN:=1;
  while(true)
  do
    nb1:=1+TheN*(TheN-1)/2;
    if nb1=nb then
      return TheN;
    fi;
    TheN:=TheN+1;
    if TheN>100 then
      Print("There is shit here\n");
      Print(NullMat(5));
    fi;
  od;
end;






ListOfVertices:=function(TheDistVect)
  local n, TheDistMat, TheGramMat, EXT, CP, TheCVP, ListVectReturn, eVect, TheSum;
  n:=FincTheN(Length(TheDistVect));
  TheDistMat:=DistanceVectorToDistanceMatrix(TheDistVect, n);
  TheGramMat:=RemoveFractionMatrix(DistanceMatrixToGramMatrix(TheDistMat));
  EXT:=[ListWithIdenticalEntries(n, 0)];
  Append(EXT, IdentityMat(n));
  CP:=CenterRadiusDelaunayPolytopeGeneral(TheGramMat, EXT);
  TheCVP:=CVPVallentinProgram(TheGramMat, CP.Center{[2..n]});
  if TheCVP.TheNorm<>CP.SquareRadius then
    Print("We find a shorter vector than it should have been\n");
    Print(NullMat(5));
  fi;
  ListVectReturn:=[];
  for eVect in TheCVP.ListVect
  do
    TheSum:=1-Sum(eVect);
    Add(ListVectReturn, Concatenation([TheSum], eVect));
  od;
  return ListVectReturn;
end;



