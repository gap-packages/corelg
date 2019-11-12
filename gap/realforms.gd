
##############################################################################
##
## Vogan diagram
##
DeclareAttribute( "VoganDiagram", IsLieAlgebra );
DeclareAttribute( "CartanName", IsLieAlgebra );
DeclareCategory( "IsVoganDiagramOfRealForm", IsObject );
DeclareCategoryCollections( "IsVoganDiagramOfRealForm" );
DeclareCategoryFamily( "IsVoganDiagramOfRealForm" );
DeclareAttribute("CanonicalGenerators",IsVoganDiagramOfRealForm);
DeclareAttribute("BasisOfSimpleRoots", IsVoganDiagramOfRealForm);
DeclareAttribute("MovedPoints",IsVoganDiagramOfRealForm);
DeclareAttribute("PermInvolution",IsVoganDiagramOfRealForm);
DeclareAttribute("Signs",IsVoganDiagramOfRealForm);
DeclareAttribute("CartanMatrix",IsVoganDiagramOfRealForm);
DeclareAttribute("CoefficientsOfSigmaAndTheta",IsVoganDiagramOfRealForm);
#DeclareGlobalFunction("ConstructVoganDiagramOfRealForm");

################################################################################
##
## Satake diagram
##
DeclareAttribute( "SatakeDiagram", IsLieAlgebra );
#DeclareAttribute( "STKDTA", IsLieAlgebra );
DeclareCategory( "IsSatakeDiagramOfRealForm", IsObject );
DeclareAttribute("BasisOfSimpleRoots", IsSatakeDiagramOfRealForm);
DeclareAttribute("ThetaInvolution",IsSatakeDiagramOfRealForm);
DeclareAttribute("CompactSimpleRoots",IsSatakeDiagramOfRealForm);
DeclareAttribute("CartanMatrix",IsSatakeDiagramOfRealForm);


##############################################################################
##
#F   IdRealForm( <L> )
##
##   returns the id of the real form
##
DeclareGlobalFunction( "IdRealForm" );

##############################################################################
##
#F   RealFormsInformation( <t>, <r> )
##
##   lists information about real forms of type t and rank r
##
DeclareGlobalFunction( "RealFormsInformation" );

##############################################################################
##
#F   NumberRealForms( <t>, <r> )
##
##   lists information about real forms of type t and rank r
##
DeclareGlobalFunction( "NumberRealForms" );

##############################################################################
##
#F   RealFormById( <t>, <r>, <id> [ <F> ] )
#F   RealFormById( <[t,r,id]> )
##
##   returns real form with given ID, defined over <F>
##
DeclareGlobalFunction( "RealFormById");

##############################################################################
##
#F   AllRealForms( <t>, <r> )
##
##   returns all real forms of type t and rank r
##
DeclareGlobalFunction( "AllRealForms");


##############################################################################
##
#A   RealFormParameters( <L> )
##
##   returns the list [type, n, signs, perm] describing the real form
##
DeclareAttribute( "RealFormParameters", IsLieAlgebra );


##############################################################################
##
#P   IsCompactForm( <L> )
##
##   returns if the real form is compact
##
DeclareProperty( "IsCompactForm", IsLieAlgebra );

##############################################################################
##
#P   IsRealification( <L> )
##
##   returns if the real form is a realification of a complex simple LA
##
DeclareProperty( "IsRealification", IsLieAlgebra );

##############################################################################
##
#P   IsRealFormOfInnerType( <L> )
##
##   returns if the real form is defined by an inner involutive automorphism
##
DeclareProperty( "IsRealFormOfInnerType", IsLieAlgebra );


##############################################################################
##
#F  IsomorphismOfRealSemisimpleLieAlgebras( <L>, <K> );
##
DeclareGlobalFunction( "IsomorphismOfRealSemisimpleLieAlgebras" );

##############################################################################
##
#A   CartanSubalgebrasOfRealForm( <L> )
##
DeclareAttribute( "CartanSubalgebrasOfRealForm", IsLieAlgebra );

##############################################################################
##
#A   NameRealForm( <L> )
##
DeclareAttribute( "NameRealForm", IsLieAlgebra );

#############################################################################
##
#O   MaximalReductiveSubalgebras( <type>, <rank>, <no> )
##
DeclareOperation( "MaximalReductiveSubalgebras", [ IsString, IsInt, IsInt ] );
