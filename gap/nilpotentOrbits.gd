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
#F  CarrierAlgebraOfNilpotentOrbit( <L>, <O> )
##  
## returns the carrier algebra of a nilpotent orbit (wrt h-grading)
##  
DeclareGlobalFunction( "CarrierAlgebraOfNilpotentOrbit" );