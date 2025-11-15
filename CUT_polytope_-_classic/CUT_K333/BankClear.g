ThePrefix:="./TheWork/TheBank/";

FuncKeep:=function(nbVert, i)
  if nbVert<50 then
    return false;
  else
    return true;
  fi;
end;

Bank_Restructing(ThePrefix, FuncKeep);
