#PL:=ArchimedeanPolyhedra("Dodecahedron");
PL:=ArchimedeanPolyhedra("Cube");
#PL:=ArchimedeanPolyhedra("Octahedron");
#PL:=ArchimedeanPolyhedra("Icosahedron");
TheGRA:=PlanGraphToGRAPE(PL);

#TheGRA:=CompleteGraph(Group(()), 5);
#TheGRA:=CompleteGraph(Group(()), 4);
#TheGRA:=CompleteGraph(Group(()), 3);


DoTestCycle:=false;
if DoTestCycle then
  GRPvert:=AutGroupGraph(TheGRA);
  RecOpt:=rec(TheLen:=-1,
              RequireIrreducible:=true,
              RequireIsometricCycleFixedLength:=false,
              RequireIsometric:=false);
  TheRes:=GRAPH_EnumerateCycles(TheGRA, GRPvert, RecOpt);
fi;

DoCreatePolytopes:=true;
if DoCreatePolytopes then
  RecTotal:=CMC_GetCutPolytopeOriented(TheGRA);
          

fi;
