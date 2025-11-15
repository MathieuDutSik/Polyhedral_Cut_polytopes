Read("FunctionSet.g");
n:=12;

ListFacet:=ReadAsFunction("SomeFacetCUT12")();

ExhaustiveSearchFromLowIncidence(n, ListFacet);

