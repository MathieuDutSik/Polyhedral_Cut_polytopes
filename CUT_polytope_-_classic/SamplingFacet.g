Read("FunctionSet.g");
#
# there is a conjecture that every facet of CUTn is adjacent
# to a triangle inequality facet.
# we are searching a counterexample to this conjecture by doing
# random movement in the graph of facets

for n in [11..20]
do
  eRec:=GetRecordForOperations(n);
  eFileSave:=Concatenation("DATA/RecordDesc", String(n));
  SaveDataToFile(eFileSave, eRec);
  #
  ListCand:=SamplingStandard(eRec.EXT);
  eFileSave:=Concatenation("DATA/RecordFAC", String(n));
  SaveDataToFile(eFileSave, ListCand);
od;

