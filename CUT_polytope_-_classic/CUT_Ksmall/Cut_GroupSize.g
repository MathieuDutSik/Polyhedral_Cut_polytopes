GetDualDescInfo:=function(n)
  local GRA, eRec, EXT;
  GRA:=GRAPH_GetMultiComplement([2,n]);
  eRec:=CMC_GetCutPolytope(GRA);
  EXT:=eRec.ListVect;
  GRP:=LinPolytope_Automorphism(EXT);
  eRestrictSize:=2^(n+1)*Factorial(n);
  eConjSize:=eRestrictSize*2^(n)*Factorial(n);
  eQuot:=Order(GRP)/eConjSize;
  Print("n=", n, " |GRP|=", Order(GRP), " eQuot=", eQuot, "\n");

end;

for n in [3..10]
do
  GetDualDescInfo(n);
od;
