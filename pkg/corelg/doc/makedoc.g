corelg.BuildManual:=function()
   MakeGAPDocDoc( ".", "manual", [], "corelg", ".", "MathJax" );;
end;

corelg.makeManual:=function(path,pathgap)
   MakeGAPDocDoc( path, "manual", [], "corelg", pathgap, "MathJax" );;
end;
