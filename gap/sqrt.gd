
################################
# SQRTFIELD ELEMENTS
DeclareCategory( "IsSqrtFieldElement", 
                 IsMultiplicativeElementWithInverse
                 and IsAdditiveElementWithInverse and IsZDFRE);

DeclareCategoryFamily( "IsSqrtFieldElement" );
DeclareCategoryCollections( "IsSqrtFieldElement" );
DeclareCategoryFamily( "IsSqrtFieldElementCollection" );
DeclareCategoryCollections( "IsSqrtFieldElementCollection" );
DeclareRepresentation( "IsSqrtFieldElementRep", 
                        IsPositionalObjectRep,
                        [1] );

fam_SqrtFieldElt := NewFamily("SqrtFieldElt",IsSqrtFieldElement);
SqrtFieldType    := NewType( fam_SqrtFieldElt,
                    IsSqrtFieldElement and IsSqrtFieldElementRep );


################################
# SQRTFIELD
DeclareCategory( "IsSqrtField", IsField );

BindGlobal( "SqrtField",  Objectify( NewType( 
                 CollectionsFamily( fam_SqrtFieldElt ),
                 IsField and
                 IsAttributeStoringRep and 
                 IsSqrtField), rec() ) );

SetName( SqrtField, "SqrtField" );
SetIsLeftActedOnByDivisionRing( SqrtField, true );
SetSize( SqrtField, infinity );
SetIsFiniteDimensional( SqrtField, false);
SetLeftActingDomain( SqrtField, GaussianRationals );
SetCharacteristic( SqrtField, 0);
SetPrimeField( SqrtField, Rationals );
SqrtFieldFam:= ElementsFamily( FamilyObj( SqrtField ) );


##################################
# GLOBAL FUNCTIONS AND ATTRIBUTES
DeclareGlobalFunction("Sqroot");
DeclareGlobalFunction("SqrtFieldIsGaussRat");
DeclareGlobalFunction("SqrtFieldMakeRational");
DeclareGlobalFunction("SqrtFieldMinimalPolynomial");
DeclareGlobalFunction("SqrtFieldEltByRationalSqrt");
DeclareGlobalFunction("SqrtFieldEltRealAndComplexPart");
DeclareGlobalFunction("IsPosSqrtFieldElt");
DeclareGlobalFunction("CoefficientsOfSqrtFieldElt");
DeclareGlobalFunction("SqrtFieldEltByCoefficients");
DeclareGlobalFunction("SqrtFieldEltToCyclotomic");
DeclareGlobalFunction("SqrtFieldEltByCyclotomic");
DeclareGlobalFunction("SqrtFieldPolynomialToRationalPolynomial");
DeclareGlobalFunction("SqrtFieldRationalPolynomialToSqrtFieldPolynomial");
DeclareAttribute( "AbsoluteValue" ,  IsSqrtFieldElement  );
