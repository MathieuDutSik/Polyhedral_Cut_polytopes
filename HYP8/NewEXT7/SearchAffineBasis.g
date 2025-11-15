R:=ReadAsFunction("DATA/NewExt7")();
EXT:=R.EXT;
SymGrp:=ReadAsFunction("SymGrp")();

OrbAffBas:=FindAllOrbitAffineBasis(EXT, SymGrp);
SaveDataToFile("DATA/OrbitAffineBasis", OrbAffBas);
