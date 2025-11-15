EXT:=ReadAsFunction("EXT64.gap")();

ScalMat:=__VectorConfigurationFullDim_ScalarMat(EXT);

eLineTot:=[];
for eLine in ScalMat
do
  Append(eLineTot, eLine);
od;
eColl:=Collected(eLineTot);
