# CoReLG, chapter 2
#
# DO NOT EDIT THIS FILE - EDIT EXAMPLES IN THE SOURCE INSTEAD!
#
# This file has been generated by AutoDoc. It contains examples extracted from
# the package documentation. Each example is preceded by a comment which gives
# the name of a GAPDoc XML file and a line range from which the example were
# taken. Note that the XML file in turn may have been generated by AutoDoc
# from some other input.
#
gap> START_TEST( "corelg01.tst");

# doc/manual.xml:226-245
gap> F := SqrtField;
SqrtField
gap> IsField( F ); LeftActingDomain( F ); Size( F ); Characteristic( F );
true
GaussianRationals
infinity
0
gap> one := One( F );
1
gap> 2 in F; 2*one in F; 2*E(4)*one in F;
false
true
true
gap> a := 2/3*E(4)*one;; 
gap> a in SqrtField; a in GaussianRationals; SqrtFieldIsGaussRat( a );
true
false
true

# doc/manual.xml:274-285
gap> Sqroot(-(2*3*4)/(11*13)); Sqroot(245/15); Sqroot(16/9);
2/143*E(4)*Sqroot(858)
7/3*Sqroot(3)
4/3
gap> a := 2+Sqroot(7)+Sqroot(99);
2 + Sqroot(7) + 3*Sqroot(11)
gap> CoefficientsOfSqrtFieldElt(a);
[ [ 2, 1 ], [ 1, 7 ], [ 3, 11 ] ]
gap> SqrtFieldEltByCoefficients([[2,9],[1,7],[E(4),13]]);
6 + Sqroot(7) + E(4)*Sqroot(13)

# doc/manual.xml:301-312
gap> SqrtFieldEltToCyclotomic( Sqroot(2) );
E(8)-E(8)^3
gap> SqrtFieldEltToCyclotomic( Sqroot(2)+E(4)*Sqroot(7) );
E(56)^5+E(56)^8+E(56)^13-E(56)^15+E(56)^16-E(56)^23-E(56)^24+E(56)^29-E(56)^31
 +E(56)^32+E(56)^37-E(56)^39-E(56)^40+E(56)^45-E(56)^47-E(56)^48+E(56)^53
 -E(56)^55
gap> SqrtFieldEltByCyclotomic( E(8)-E(8)^3 );
Sqroot(2)
gap> SqrtFieldEltByCyclotomic( 3*E(4)*Sqrt(11)-2/4*Sqrt(-13/7) );
3*E(4)*Sqroot(11) + (-1/14*E(4))*Sqroot(91)

# doc/manual.xml:320-341
gap> a := Sqroot( 2 ) + 3 * Sqroot( 3/7 ); b := Sqroot( 21 ) - Sqroot( 2 );
Sqroot(2) + 3/7*Sqroot(21)
(-1)*Sqroot(2) + Sqroot(21)
gap> a + b; a * b; a - b;
10/7*Sqroot(21)
7 + 4/7*Sqroot(42)
2*Sqroot(2) + (-4/7)*Sqroot(21)
gap> c := ( a - b )^-2;
91/8 + 7/4*Sqroot(42)
gap> a := Sum( List( [2,3,5,7], x -> Sqroot( x ) ) );
Sqroot(2) + Sqroot(3) + Sqroot(5) + Sqroot(7)
gap> b := a^-1; a*b;                                  
37/43*Sqroot(2) + (-29/43)*Sqroot(3) + (-133/215)*Sqroot(5) + 27/43*Sqroot(
7) + 62/215*Sqroot(30) + (-10/43)*Sqroot(42) + (-34/215)*Sqroot(70) + 22/
215*Sqroot(105)
1
gap> ComplexConjugate(Sqroot(17)+Sqroot(-7));
(-E(4))*Sqroot(7) + Sqroot(17)
gap> Random( SqrtField );
E(4) + (-7/6+1/4*E(4))*Sqroot(2) + (-3/2*E(4))*Sqroot(3)

# doc/manual.xml:344-358
gap> m:=[[Sqroot(2),Sqroot(3)],[Sqroot(2),Sqroot(5)],[1,0]]*One(SqrtField);
[ [ Sqroot(2), Sqroot(3) ], [ Sqroot(2), Sqroot(5) ], [ 1, 0 ] ]
gap> NullspaceMat(m);
[ [ (-5/4)*Sqroot(2) + (-1/4)*Sqroot(30), 3/4*Sqroot(2) + 1/4*Sqroot(30), 1 ] ]
gap> RankMat(m);
2
gap> m := [[Sqroot(2),Sqroot(3)],[Sqroot(2),Sqroot(5)]];  
[ [ Sqroot(2), Sqroot(3) ], [ Sqroot(2), Sqroot(5) ] ]
gap> Determinant( m );  DefaultFieldOfMatrix( m );
(-1)*Sqroot(6) + Sqroot(10)
SqrtField
gap> x := Indeterminate( SqrtField, "x" );; f := x^2+x+1;
x^2+x+1

# doc/manual.xml:388-400
gap> F := SqrtField;; one := One( SqrtField );;                 
gap> x := Indeterminate( F, "x" );; f := x^5 + 4*x^3 + E(4)*one*x;
x^5+4*x^3+E(4)*x
gap> SqrtFieldPolynomialToRationalPolynomial(f);
x_1^5+4*x_1^3+E(4)*x_1
gap> SqrtFieldRationalPolynomialToSqrtFieldPolynomial(last);
x^5+4*x^3+E(4)*x
gap> f := x^2-1;; Factors(f);
[ x-1, x+1 ]
gap> f := x^2+1;; Factors(f);
[ x^2+1 ]
gap> STOP_TEST("corelg01.tst", 1 );
