##############################################################################
##
#F   RootsystemOfCartanSubalgebra( <L> )
#F   RootsystemOfCartanSubalgebra( <L>, <H> )
##
##   <L> is a semisimple lie algebra over Gaussian rationals or SqrtField;
##   this function returns a rootsystem of <L> with respect to <H>, and
##   with repect to CartanSubalgebra(<L>) if <H> is not provided
##
DeclareGlobalFunction( "RootsystemOfCartanSubalgebra" );

##############################################################################
##
#F   RootSystemOfZGradedLieAlgebra( <L>, <gr> );
#F   RootSystemOfZGradedLieAlgebra( <L>, <gr>, <H> );
##
##   <L> is a semisimple lie algebra over Gaussian rationals or SqrtField
##   with Z-grading <gr>, which is a record with entries g0, gp, gn;
##   this function returns a rootsystem of <L> with respect to <H>, and
##   with repect to CartanSubalgebra(<L>) if <H> is not provided, such that
##   the simple roots lie in <gr>.gp[1] and <gr>.g0
##
DeclareGlobalFunction( "RootSystemOfZGradedLieAlgebra" );


##############################################################################
##
#F LieAlgebraIsomorphismByCanonicalGenerators( <L1>, <R1>, <L2>, <R2> )
##
## <L1> and <L2> both are semisimple lie algebras over Gaussian rationals or SqrtField
## and either <R1> and <R2> both are canonical generators of <L1> and <L2> defining 
## the same Cartan Matrix, of <R1> and <R2> are rootsystems or Cartan subalgebras of 
## <L1> and <L2>, respectively; this functions constructs an isomorphism from
## <L1> to <L2> by mapping canonical generators onto canonical generators.
## Attention: This function does not check whether the map actually is a
##            Lie isomorphisms! 
##
DeclareGlobalFunction( "LieAlgebraIsomorphismByCanonicalGenerators" );

##############################################################################
##
#F RegularCarrierAlgebraOfSL2Triple( <L>, <sl2> )
##
## <L> is a semisimple lie algebra over Gaussian rationals or SqrtField
## and <sl2> is an SL2-triple in <L> of the form [f,h,e] with ef=h, he=2e, hf=-2f;
## this function returns the Z-graded carrier algebra of <sl2> normalised by
## CartanSubalgebra(<L>).
##
DeclareGlobalFunction( "RegularCarrierAlgebraOfSL2Triple" );


##############################################################################
##
#A ChevalleyBasis( <R> )
#
DeclareAttribute( "ChevalleyBasis", IsRootSystem );

##############################################################################
##
#P IsIsomorphismOfLieAlgebras( <iso> )
#
DeclareProperty( "IsIsomorphismOfLieAlgebras", IsAlgebraHomomorphism );