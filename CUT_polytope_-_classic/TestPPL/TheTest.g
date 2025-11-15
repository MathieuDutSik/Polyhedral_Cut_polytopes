n:=6;
eRec:=CutPolytope(n);

EXT:=eRec.EXT;

GroupExt:=LinPolytope_Automorphism(EXT);

ThePath:="./TheWork/";
TheProg:="PPL";
#TheProg:="CDD";

TheList:=__DualDescriptionDoubleDescMethod_Reduction(EXT, GroupExt, ThePath, TheProg);
