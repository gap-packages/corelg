DeclareAttribute( "NilpotentOrbitsOfRealForm", IsLieAlgebra);

##############################################################################
##
#A   RealCayleyTriple( <O> )
##
##   returns the real Cayley triple [f,h,e] attached to the real nilpotent 
##   orbit <O>, where ef=h, he=2e, hf=-2f
##
DeclareAttribute( "RealCayleyTriple", IsNilpotentOrbit);


##############################################################################
##
#A   Invariants( <O> )
##
##   returns certain invariants of the real nilpotent orbits <O>
DeclareAttribute( "Invariants", IsNilpotentOrbit);             

##############################################################################
##
#F  RealNilpotentOrbitsFromDatabase( <L> )
##  reads and returns all real nilpotent orbits in the real form L;
##  at the moment, the database contains A2-A8, B2-B10, C2-C10, D4-D8, F4, G2, E6-E8
##  
DeclareGlobalFunction( "RealNilpotentOrbitsFromDatabase" );

##############################################################################
##
#F  CarrierAlgebraOfNilpotentOrbit( <L>, <O> )
##  
## returns the carrier algebra of a nilpotent orbit (wrt h-grading)
##  
DeclareGlobalFunction( "CarrierAlgebraOfNilpotentOrbit" );