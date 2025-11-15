FuncDistMat:=function(EXT)
  local AffBas, Hcoord, Mat, HypDim, DimHyp, V, VZ, VZ2;
  AffBas:=CreateAffineBasis(EXT);
  Hcoord:=HypermetricCoordinates(AffBas, EXT);
  Mat:=MatrixHypermetricInequalities(Hcoord);
  HypDim:=RankMat(EXT);
  DimHyp:=1+HypDim*(HypDim-1)/2;
  V:=ListWithIdenticalEntries(DimHyp, 0);
  V[1]:=1;
  Add(Mat, V);
  VZ:=NullspaceMat(TransposedMat(Mat))[1];
  VZ2:=RemoveFraction(VZ);
  return RemoveFractionMatrix(DistanceConstructionDelaunay(VZ2, Hcoord));
end;

GetGramMatrix:=function(EXT)
  local AffBas, Hcoord, HcoordIneq, V, DistVect, DistMat;
  AffBas:=CreateAffineBasis(EXT);
  Hcoord:=HypermetricCoordinates(AffBas, EXT);
  HcoordIneq:=List(Hcoord, FromHypermetricVectorToHypermetricFace);
  V:=ListWithIdenticalEntries(Length(HcoordIneq[1]), 0);
  V[1]:=1;
  Add(HcoordIneq, V);
  DistVect:=NullspaceMat(TransposedMat(HcoordIneq))[1];
  DistMat:=DistanceVectorToDistanceMatrix(DistVect, Length(EXT[1]));
  return RemoveFractionMatrix(DistanceMatrixToGramMatrix(DistMat));
end;
