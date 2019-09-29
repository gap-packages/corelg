# CoReLG, chapter 4
#
# DO NOT EDIT THIS FILE - EDIT EXAMPLES IN THE SOURCE INSTEAD!
#
# This file has been generated by AutoDoc. It contains examples extracted from
# the package documentation. Each example is preceded by a comment which gives
# the name of a GAPDoc XML file and a line range from which the example were
# taken. Note that the XML file in turn may have been generated by AutoDoc
# from some other input.
#
gap> START_TEST( "corelg03.tst");

# doc/manual.xml:875-881
gap> L:= RealFormById( "F", 4, 3 );;
gap> no:= NilpotentOrbitsOfRealForm( L );;
#I CoReLG: read database of real triples ... done
gap> no[1];
<nilpotent orbit in Lie algebra>

# doc/manual.xml:891-905
gap> L:= RealFormById( "F", 4, 2 );;
gap> no:= NilpotentOrbitsOfRealForm( L );;
gap> o:= no[10];
<nilpotent orbit in Lie algebra>
gap> t:=RealCayleyTriple(o);;
gap> theta:= CartanDecomposition(L).CartanInv;
function( v ) ... end
gap> theta(t[1]) = -t[3];
true
gap> theta(t[2]) = -t[2];
true
gap> t[3]*t[1] = t[2];
true

#
gap> STOP_TEST("corelg03.tst", 1 );