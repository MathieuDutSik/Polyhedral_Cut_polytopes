PL:=ArchimedeanPolyhedra("Dodecahedron");
#PL:=ArchimedeanPolyhedra("Icosahedron");

GRA:=PlanGraphToGRAPE(PL);

RecCMC:=CMC_GetCutPolytope(GRA);

ListInfo:=RecCMC.GetMET_info_all();

RecListOrbitCycle:=ListInfo.RecListOrbitCycle;

ListSizes:=List(RecListOrbitCycle.FinalListOrb, Length);

