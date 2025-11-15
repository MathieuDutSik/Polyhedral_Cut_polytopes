GRA1:=GRAPH_Cycle(4);
GRA2:=GRAPH_Pyramid(GRA1);
GRAfinal1:=GRAPH_Pyramid(GRA2);




ListSize:=[1,1,1,2,2];
GRAfinal2:=GRAPH_GetMultiComplement(ListSize);


test:=IsIsomorphicGraph(GRAfinal1, GRAfinal2);
