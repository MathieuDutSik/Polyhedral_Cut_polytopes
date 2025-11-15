n:=6;

eRec:=Cone_QMET(n);

EXT:=ColumnReduction(eRec.FAChom).EXT;
DDA:=DualDescriptionAdjacencies(EXT);
Print("We have DDA\n");

TheDiamSkel:=Diameter(DDA.SkeletonGraph);
TheDiamRidge:=Diameter(DDA.RidgeGraph);

