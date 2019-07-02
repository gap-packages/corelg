#####################################################################################
#
#  init.g
#
#
# The package CoReLG is free software; you can redistribute it and/or modify it under the 
# terms of the GNU General Public License as published by the Free Software Foundation; 
# either version 2 of the License, or (at your option) any later version. 


################################
# INFO CLASSES
DeclareInfoClass("InfoSqrtField");
SetInfoLevel(InfoSqrtField,0);
DeclareInfoClass("InfoCorelg");
SetInfoLevel(InfoCorelg,0);

################################
# record for functions:
corelg := rec();


################################
# READ gd FILES
ReadPackage( "corelg", "gap/sqrt.gd" );
ReadPackage( "corelg", "gap/realforms.gd" );
ReadPackage( "corelg", "gap/carrierrtsys.gd" );
ReadPackage( "corelg", "gap/nilpotentOrbits.gd" );
ReadPackage( "corelg", "gap/cartandecomp.gd" );
ReadPackage( "corelg", "gap/realtheta.gd" );
ReadPackage( "corelg", "gap/realweyl.gd" );
