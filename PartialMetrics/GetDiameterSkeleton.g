n:=5;
eRec:=ConePoly_PartialMetrics(n);


FAC:=eRec.ListIneqPoly;
DDA:=DualDescriptionAdjacencies(FAC);
nbVert:=DDA.RidgeGraph.order;
eDiam:=Diameter(DDA.RidgeGraph);
Print("nbVert=", nbVert, " eDiam=", eDiam, "\n");
