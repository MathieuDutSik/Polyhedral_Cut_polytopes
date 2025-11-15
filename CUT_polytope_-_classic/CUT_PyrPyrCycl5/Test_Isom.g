GRA1:=GRAPH_Cycle(5);
GRAfinal1:=ComplementGraph(GRAPH_ZeroExtension(GRA1, 2));


GRA1:=GRAPH_Cycle(5);
GRA2:=GRAPH_Pyramid(GRA1);
GRAfinal2:=GRAPH_Pyramid(GRA2);


test:=IsIsomorphicGraph(GRAfinal1, GRAfinal2);
