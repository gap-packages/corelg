
##############################################################################
##
#A   CartanDecomposition( <L> )
##
##   returns the record with entries <K>, <P>, and <CartanInv>
##
DeclareAttribute( "CartanDecomposition", IsLieAlgebra );

##############################################################################
##
#A   RealStructuren( <L> )
##
##   returns the complex conjugation wrt the real Lie algebra <L>
##
DeclareAttribute( "RealStructure", IsLieAlgebra );

##############################################################################
##
#A MaximallyNonCompactCartanSubalgebra ( <L> )
##
##
DeclareAttribute( "MaximallyNonCompactCartanSubalgebra", IsLieAlgebra );

##############################################################################
##
#A MaximallyCompactCartanSubalgebra ( <L> )
##
##
DeclareAttribute( "MaximallyCompactCartanSubalgebra", IsLieAlgebra );

##############################################################################
##
#A corelgCompactDimOfCSA( <L> )
##
##
DeclareAttribute( "corelgCompactDimOfCSA", IsLieAlgebra );

##############################################################################
##
#F CompactDimensionOfCartanSubalgebra( <L> )
#F CompactDimensionOfCartanSubalgebra( <L>, <H> )
##
##
DeclareGlobalFunction( "CompactDimensionOfCartanSubalgebra" );