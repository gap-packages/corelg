corelg.BuildManual:=function()
   MakeGAPDocDoc( ".", "manual", [], "corelg" );;
end;

corelg.makeManual:=function(path)
   MakeGAPDocDoc( path, "manual", [], "corelg" );;
end;
