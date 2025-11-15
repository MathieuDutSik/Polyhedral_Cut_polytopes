GRA1:=GRAPH_Cycle(4);
GRAfinal1:=ComplementGraph(GRAPH_ZeroExtension(GRA1, 3));


ListSize:=[1,1,1,2,2];
GRAfinal2:=GRAPH_GetMultiComplement(ListSize);


test:=IsIsomorphicGraph(GRAfinal1, GRAfinal2);
