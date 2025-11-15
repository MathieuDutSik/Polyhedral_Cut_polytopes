#
#
# We want to check whether the symmetry of the faces are entirely
# determined by the stabilizers, i.e. whether it make sense to
# have additional stabilizer computations.
BankPath:="./TheWork/TheBank/";

nbRec:=Bank_GetSize(BankPath);
ListGRPsizeTotal:=[];
ListGRPsizeStab:=[];

FileSaveGRP:="TheData_GRP";
GRP:=ReadAsFunction(FileSaveGRP)();

PL:=ArchimedeanPolyhedra("Icosahedron");
GRA:=PlanGraphToGRAPE(PL);

eRec:=CMC_GetCutPolytope(GRA);
EXTtotal:=eRec.ListVect;


ListExceptional:=[];
for iRec in [1..nbRec]
do
  Print("iRec=", iRec, "/", nbRec, "\n");
  eFileEXT:=Concatenation(BankPath, "AccountEXT_", String(iRec));
  eInfoEXT:=ReadAsFunction(eFileEXT)();
  eSizeTot:=Order(eInfoEXT.Group);
  Add(ListGRPsizeTotal, eSizeTot);
  #
  eSet:=Set(List(eInfoEXT.EXT, x->Position(EXTtotal, x)));
  eStab:=Stabilizer(GRP, eSet, OnSets);
  eStabRed:=SecondReduceGroupAction(eStab, eSet);
  eSizeStab:=Order(eStabRed);
  Add(ListGRPsizeStab, eSizeStab);
  eQuot:=eSizeTot/eSizeStab;
  if eSizeStab<>eSizeTot then
    eExcept:=rec(iRec:=iRec, eSizeTot:=eSizeTot, eSizeStab:=eSizeStab, eQuot:=eQuot);
    Add(ListExceptional, eExcept);
  fi;
od;
