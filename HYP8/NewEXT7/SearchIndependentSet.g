R:=ReadAsFunction("DATA/NewExt7")();
EXT:=R.EXT;
GRP:=LinPolytope_Automorphism(EXT);

OrbIndep:=FindAllOrbitIndependent(EXT, GRP);
SaveDataToFile("DATA/OrbitIndependentSet", OrbIndep);
