Read("FunctionSet.g");
n:=11;

ListFacet:=ReadAsFunction("SomeFacetCUT11")();

ExhaustiveSearchFromLowIncidence(n, ListFacet);

