#
#
# finding symmetry group
R:=ReadAsFunction("DATA/NewExt7")();
EXT:=R.EXT;
FAC:=R.FAC;
GraphExt:=R.GraphExt;
GraphFac:=R.GraphFac;
distvector:=R.distvector;


SymGrp:=AutomorphismExtremeDelaunay(EXT, FAC, GraphExt, GraphFac, distvector);
SaveDataToFile("SymGrp", SymGrp);