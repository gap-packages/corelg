# Main functions:
#
#
#F   CarrierAlgsForNilpOrbsInZGrading( type, rank, d )
#F   CarrierAlgsForNilpOrbsInZmGrading( type, rank, m0, str, num )
#
#   Gives a record containing the carrier algebras of the real theta group
#   specified by the input. Explanation of the input:
#
#   type: type of the Lie algebra where everything happens,
#   rank: its rank,
#   d : (for Z grading) the degrees of the simple roots,
#   m0 : the order of the automorphism defining the grading,
#   str: "inner" or "outer", the first when an inner automorphism
#        defines the grading, the second otherwise,
#   num : the num-th automorphism in the list 
#                 FiniteOrderInnerAutomorpisms( type, rank, m0 ),
#         or 
#                 FiniteOrderOuterAutomorphisms( type, rank, m0, 2 ) 
#   is used to define the grading.
#
#   The output is a record with the following fields:
# 
#
#   L : the Lie algebra,
#   grad : the grading that was used (different format for Z-grading and 
#          Z/mZ-grading),
#   Hs : the Cartan subalgebras of g_0 that are used,
#   L0 : subalgebra g_0,
#   cars : the carrier algebras, this is a list of lists; for each Cartan 
#          subalgebra of g_0 there is one list: the first corresponds to the 
#          split Cartan subalgebra, and so has just the complex carrier 
#          algebras (which are also real), the other lists contain lists as 
#          well, for each complex carrier algebra (i.e., for each element of 
#          the first list) there is a list containing the real carrier 
#          algebras which are strongly H_i-regular, and over the complexes
#          conjugate to the given complex carrier algebra.
#          Furthermore, a carrier algebra is given by a record, containing 
#          the fields g0, gp (positive degree), and gn (negative degree). 
#
#
#
#
#F    corelg.torusparam( L, H )
#
#     Here L is the "ambient" Lie algebra, H is a torus in it
#     (example, with r the record output by the previous function,
#     L=r.L, H = Intersection( r.Hs[2], CartanDecomposition(L).K )).
#     The output is a toral subgroup of Aut(L), parametrised.
#     Not to be looked at, really... (but to be used in the next function).
#
#F    corelg.resmat( L, T, c )
#
#     Here L is as previous, T is the output of the previous function, c
#     is a carrier algebra. The output of this function is a record, with 
#     fields:  bas:= vv, mats:= mats, densep
#     
#     mats : the matrices of the torus T, restricted to the space c_1,
#     bas : the basis of the space c_1 used for this,
#     densep: an element in c_1 is in general position iff its coefficients
#             wrt the basis bas, when substituted in this polynomial 
#             give non-zero.
#
#F    corelg.IsSupport( L, L0, c, e )
#
#     Here L = r.L, L0 = r.L0, c a carrier algebra, e an element in c_1,
#     in general position. This function checks whether c can arise as
#     a carrier algebra of e.  
#
#F    corelg.expmat( L, u, c1 )
# 
#     Here L is as before, c1 a basis of the 1-component of a carrier algebra,
#     u a nilpotent element of L stabilising c1 (for example coming from
#     c0). This function returns the exp of the matrix of ad u restricted
#     to c1, with parameter s (so exp( s ad u )).
#
#
#
#   EXAMPLE:
#
#gap> r:= CarrierAlgsForNilpOrbsInZmGrading( "F", 4, 2, "inner", 2 );;
#gap> c:=r.cars[1][25];
#rec( g0 := [ v.3, v.27, v.49, v.50+v.52, v.51, v.52 ], 
#  gn := [ [ v.6, v.9, v.13, v.35, v.38, v.40 ], 
#      [ v.25, v.28, v.29, v.31, v.34 ], [ v.2, v.41, v.42, v.43 ], 
#      [ v.32, v.36, v.48 ], [ v.44, v.45, v.46 ], [ v.39 ], [ v.47 ] ], 
#  gp := [ [ v.11, v.14, v.16, v.30, v.33, v.37 ], 
#      [ v.1, v.4, v.5, v.7, v.10 ], [ v.17, v.18, v.19, v.26 ], 
#      [ v.8, v.12, v.24 ], [ v.20, v.21, v.22 ], [ v.15 ], [ v.23 ] ] )
#
#  We want to classify the nilpotent orbits in c_1; first of all,
#  this is a carrier algebra relative the first (split) CSA, so we take that
#  one, and parametrise the corresponding group. In this case it is not
#  necessary to split H into compact/noncompact parts, as there only is
#  the noncompact part.
#
#gap> H:= r.Hs[1];
#<Lie algebra of dimension 4 over SqrtField>
#gap> T:= corelg.torusparam(r.L,H);;

# Now we look at the action of the torus on c_1:
#gap> rs:= corelg.resmat(r.L,T,c);
#rec( bas := [ v.11, v.14, v.16, v.30, v.33, v.37 ], 
#  densep := -2*x1^2*x3*x4^2*x6+2*x1^2*x3*x4*x5^2-4*x1*x2*x3*x4*x5*x6+
#4*x1*x2*x3*x5^3-2*x2^2*x3*x4*x6^2+2*x2^2*x3*x5^2*x6, 
#  mats := 
#    [ [ [ 1, 0, 0, 0, 0, 0 ], [ 0, a1^-1, 0, 0, 0, 0 ], [ 0, 0, 1, 0, 0, 0 ], 
#          [ 0, 0, 0, 1, 0, 0 ], [ 0, 0, 0, 0, a1, 0 ], 
#          [ 0, 0, 0, 0, 0, a1^2 ] ], 
#      [ [ 1, 0, 0, 0, 0, 0 ], [ 0, a2, 0, 0, 0, 0 ], [ 0, 0, a2, 0, 0, 0 ], 
#          [ 0, 0, 0, 1, 0, 0 ], [ 0, 0, 0, 0, 1/a2, 0 ], 
#          [ 0, 0, 0, 0, 0, 1/a2^2 ] ], 
#      [ [ 1/a3, 0, 0, 0, 0, 0 ], [ 0, 1, 0, 0, 0, 0 ], [ 0, 0, 1, 0, 0, 0 ], 
#          [ 0, 0, 0, a3, 0, 0 ], [ 0, 0, 0, 0, 1, 0 ], 
#          [ 0, 0, 0, 0, 0, 1/a3 ] ], 
#      [ [ 1/a4, 0, 0, 0, 0, 0 ], [ 0, a4^3, 0, 0, 0, 0 ], 
#          [ 0, 0, 1, 0, 0, 0 ], [ 0, 0, 0, a4^2, 0, 0 ], 
#          [ 0, 0, 0, 0, 1/a4^2, 0 ], [ 0, 0, 0, 0, 0, 1/a4^6 ] ] ] )
#
#   First we factor the densep in Magma:
#
#> P<x1,x2,x3,x4,x5,x6>:= PolynomialRing( Rationals(), 6 );
#> densep := -2*x1^2*x3*x4^2*x6+2*x1^2*x3*x4*x5^2-4*x1*x2*x3*x4*x5*x6+
#> 4*x1*x2*x3*x5^3-2*x2^2*x3*x4*x6^2+2*x2^2*x3*x5^2*x6;
#> Factorization(densep);
#[
#    <x4*x6 - x5^2, 1>,
#    <x3, 1>,
#    <x1^2*x4 + 2*x1*x2*x5 + x2^2*x6, 1>
#]
#
#  If xi also denotes the coefficients of an elt in general position,
#  then we see that x3\neq 0, for example. Also we see that we cannot have
#  x1=x2=0.
#
# Now we look at some exp-s of nilpotent elements (L.3, L.27 are in c_0):
#gap> corelg.expmat( r.L, r.L.3, rs.bas );;
#gap> Display(last);               
#[ [    1,    0,    0,    0,    0,    0 ],
#  [   -s,    1,    0,    0,    0,    0 ],
#  [    0,    0,    1,    0,    0,    0 ],
#  [    0,    0,    0,    1,  2*s,  s^2 ],
#  [    0,    0,    0,    0,    1,    s ],
#  [    0,    0,    0,    0,    0,    1 ] ]
#
#  We call this ex1, and
#gap> corelg.expmat( r.L, r.L.27, rs.bas );;
#gap> Display(last);
#[ [    1,   -s,    0,    0,    0,    0 ],
#  [    0,    1,    0,    0,    0,    0 ],
#  [    0,    0,    1,    0,    0,    0 ],
#  [    0,    0,    0,    1,    0,    0 ],
#  [    0,    0,    0,    s,    1,    0 ],
#  [    0,    0,    0,  s^2,  2*s,    1 ] ]
#
#  which we call ex2. Now if x1=0, then in ex2 we choose a nonzero s,
#  and get the new x1 nonzero, so we may assume x1\neq 0. Then by 
#  choosing an appropriate s in ex1, we get x2=0.
#  Then from the factorisation of densep we get that x4\neq 0.
#  Then again we can choose an s in ex2 such that x5 is mapped to 0
#  (and x2 remains 0). So we may assume x2=x5=0, which also
#  implies that the remaining coordinates are nonzero.
#
#  Now we act by the torus:
#gap> Product(rs.mats);
#[ [ 1/(a3*a4), 0, 0, 0, 0, 0 ], [ 0, a2*a4^3/a1, 0, 0, 0, 0 ], 
#  [ 0, 0, a2, 0, 0, 0 ], [ 0, 0, 0, a3*a4^2, 0, 0 ], 
#  [ 0, 0, 0, 0, a1/(a2*a4^2), 0 ], [ 0, 0, 0, 0, 0, a1^2/(a2^2*a3*a4^6) ] ]
# 
#  The elements on positions (2,2) and (5,5) are of no interest to us.
#  In Magma we perform the following calculation:
#
#> F<s,t,u,v>:= RationalFunctionField(Rationals(),4);
#> R<a1,a2,a3,a4>:= PolynomialRing( F, 4 );
#> r:=[a3*a4-s,a2-t,a3*a4^2-u,a1^2-v*a2^2*a3*a4^6];
#> GroebnerBasis(r);
#[
#    a1^2 - t^2*u^5*v/s^4,
#    a2 - t,
#    a3 - s^2/u,
#    a4 - u/s
#]
#
#  We see that if we want a general matrix of T to restrict to 
#  diag(s^-1,*,t,u,*,v) then we can take the choices given by the above
#  calculation for the ai. However, we see also that u*v needs to be a
#  square, otherwise there is no solution for the ai. So we get at most
#  two representatives, with coefficients [1,0,1,1,0,1], and [1,0,1,1,0,-1].
#
#  Now we have to show that they are not G_0-conjugate:
#gap> e1:= [1,0,1,1,0,1]*rs.bas;
#v.11+v.16+v.30+v.37
#gap> e2:= [1,0,1,1,0,-1]*rs.bas;
#v.11+v.16+v.30+(-1)*v.37
#gap> corelg.IsSupport( r.L, r.L0, c, e1 );
#rank of Z_0(c): 0
#rank of Z_0(sl2): 0
#true
#gap> corelg.IsSupport( r.L, r.L0, c, e2 );
#rank of Z_0(c): 0
#rank of Z_0(sl2): 0
#true
#
#  So c can be the carrier algebra of both. We compute two sl2-triples:
#
#gap> L0:= r.L0; L:= r.L;
#<Lie algebra of dimension 24 over SqrtField>
#<Lie algebra of dimension 52 over SqrtField>
#gap> t1:= SL2Triple(L,e1);
#[ (18)*v.6+(8)*v.13+(10)*v.35+(14)*v.40, 
#  (10)*v.49+(8)*v.50+(16)*v.51+(22)*v.52, v.11+v.16+v.30+v.37 ]
#gap> t2:= SL2Triple(L,e2);
#[ (18)*v.6+(-8)*v.13+(10)*v.35+(14)*v.40, 
#  (10)*v.49+(8)*v.50+(16)*v.51+(22)*v.52, v.11+v.16+v.30+(-1)*v.37 ]
#gap> Intersection( L0, LieCentralizer(L,Subalgebra(L,t1)));  
#<Lie algebra of dimension 0 over SqrtField>
#
#   No luck there....
#
#gap> Intersection( L0, LieCentralizer(L,Subalgebra(L,[t1[3]-t1[1]])));
#<Lie algebra of dimension 1 over SqrtField>
#gap> BasisVectors( Basis(last));
#[ v.1+(8)*v.3+v.7+(-120)*v.25+(-12)*v.27+(-168)*v.31 ]
#gap> h1:= last[1];
#v.1+(8)*v.3+v.7+(-120)*v.25+(-12)*v.27+(-168)*v.31
#gap> Intersection( L0, LieCentralizer(L,Subalgebra(L,[t2[3]-t2[1]])));
#<Lie algebra of dimension 1 over SqrtField>
#gap> BasisVectors( Basis(last));
#[ v.1+(-8)*v.3+v.7+(120)*v.25+(-12)*v.27+(168)*v.31 ]
#gap> h2:= last[1];
#v.1+(-8)*v.3+v.7+(120)*v.25+(-12)*v.27+(168)*v.31
#gap> ad:= AdjointMatrix( Basis(L0), h1 );;
#gap> MinimalPolynomial(ad);
#s^5+1920*s^3+589824*s
#gap> Factors(last);
#[ s, s^2+384, s^2+1536 ]
#gap> ad:= AdjointMatrix( Basis(L0), h2 );;
#gap> MinimalPolynomial(ad);
#s^5+(-1920)*s^3+589824*s
#gap> Factors(last);
#[ s, s^2+(-1536), s^2+(-384) ]
#
#  So in the second case the eigenvalues are real, whereas in the 
#  first case they are not, therefore the elements are not conjugate.

#########################################################################
#########################################################################
##############################GAP code###################################

corelg.normalizer:= function( L, U )


     local R, B, T, n, s, v, A, i, j, l, k, pos, A0, b, bas, cij;

     if not LeftActingDomain(L) = SqrtField then
        return LieNormalizer( L, U );
     fi;

    if Dimension(L) = Dimension(U) then 
       return L;
    fi;

    # We need not work if `U' knows to be an ideal in its parent `L'.
    if HasParent( U ) and IsIdenticalObj( L, Parent( U ) )
       and HasIsLeftIdealInParent( U ) and IsLeftIdealInParent( U ) then
      return L;
    fi;

    R:= LeftActingDomain( L );
    B:= Basis( L );
    T:= StructureConstantsTable( B );
    n:= Dimension( L );
    s:= Dimension( U );

    if s = 0 or n = 0 then
      return L;
    fi;

    v:= List( BasisVectors( Basis( U ) ),
              x -> Coefficients( B, x ) );

    # The equations.
    # First the normalizer part, \ldots

    A:= NullMat( n + s*s, n*s, R );
    for i in [ 1..n ] do
      for j in [ 1..n ] do
        cij:= T[i][j];
        for l in [ 1..s ] do
          for k in [ 1..Length( cij[1] ) ] do
            pos:= (l-1)*n+cij[1][k];
            A[i][pos]:= A[i][pos]+v[l][j]*cij[2][k];
          od;
        od;
      od;
    od;

    # \ldots and then the "superfluous" part.

    for k in [1..n] do
      for l in [1..s] do
        for i in [1..s] do
          A[ n+(l-1)*s+i ][ (l-1)*n+k ]:= -v[i][k];
        od;
      od;
    od;

    # Solve the equation system.
    A0:= SqrtFieldMakeRational(A);
    if A0 = false then
       b:= NullspaceMat(A);
    else
       b:= NullspaceMat(A0)*One(SqrtField);
    fi;

    # Extract the `normalizer part` of the solution.
    l:= Length(b);
    bas:= NullMat( l, n, R );
    for i in [ 1..l ] do
      for j in [ 1..n ] do
        bas[i][j]:= b[i][j];
      od;
    od;

    # Construct the generators from the coefficients list.
    bas:= List( bas, x -> LinearCombination( B, x ) );

    # Return the subalgebra.
    return SubalgebraNC( L, bas, "basis" );

end;





#################################################Real Weyl Group########
corelg.coroots:= function(L,H)

   # the coroots of the root sys of L wrt H

   local R, ch, h;

   R:= RootsystemOfCartanSubalgebra( L, H );
   ch:= ChevalleyBasis(R);
   h:= List( [1..Length(ch[1])], i -> ch[1][i]*ch[2][i] );
   Append( h, -h ); # for the negative roots...
   return h;

end;

######################

corelg.refl:= function( R, i ) # i : index of a pos root, return the perm of the corr
                        # reflection

   local rts, r, B, pm, j, s;

   rts:= Concatenation( PositiveRootsNF(R), -PositiveRootsNF(R) );
   r:= rts[i];
   B:= BilinearFormMatNF(R);
   pm:= [ ];
   for j in [1..Length(rts)] do
       s:= rts[j] - (2*rts[j]*B*r/(r*B*r))*r;
       pm[j]:= Position( rts, s );
   od;
   return PermList( pm );

end;


######################

corelg.isreal:= function( L, HK, R, i ) # test whether the i-th pos root of R is real, L a real Lie algebra
                                # R root sys wrt CSA H

    local cd, ch, x;

    cd:= CartanDecomposition(L);
    #HK:= BasisVectors( Basis( Intersection( cd.K, H ) ) );
    ch:= ChevalleyBasis(R);
    x:= ch[1][i];
    return ForAll( HK, h -> IsZero(h*x) );

end;

######################


corelg.isimag:= function( L, HP, R, i ) # test whether the i-th pos root of R is real, L a real Lie algebra
                                # R root sys wrt CSA H

    local cd, ch, x;

    cd:= CartanDecomposition(L);
    #HP:= BasisVectors( Basis( Intersection( cd.P, H ) ) );
    ch:= ChevalleyBasis(R);
    x:= ch[1][i];
    return ForAll( HP, h -> IsZero(h*x) );

end;

######################

corelg.base:= function(R,inds) # inds : indices of pos rts of a subsys, find basis

    local pr, sums, r, s, posR;

    posR:= PositiveRootsNF(R);
    pr:= posR{inds};
    sums:= [ ];
    for r in pr do 
        for s in pr do
            Add( sums, r+s );
        od;
    od;
    return List( Filtered( pr, x -> not x in sums ), y -> Position(posR,y) );

end;

######################


corelg.centrzer:= function( L, H ) 

   local R, re, br, gr, im, bi, gi, ch, hr, hi, ex, bex, gex, h, theta,
         p, gex0, cpt, gcpt, bcpt, HP, HK;

   HK:= BasisVectors( Basis( Intersection( H, CartanDecomposition(L).K ) ) );
   HP:= BasisVectors( Basis( Intersection( H, CartanDecomposition(L).P ) ) );

   R:= RootsystemOfCartanSubalgebra(L,H);

   re:= Filtered( [1..Length(PositiveRoots(R))], i ->  corelg.isreal(L,HK,R,i) );     
   br:= corelg.base( R, re );
   gr:= List( br, x -> corelg.refl( R, x ) );

   im:= Filtered( [1..Length(PositiveRoots(R))], i ->  corelg.isimag(L,HP,R,i) );     
   bi:= corelg.base( R, im );
   gi:= List( bi, x -> corelg.refl( R, x ) );

   ch:= ChevalleyBasis(R);

   cpt:= Filtered( im, i -> ch[1][i] in CartanDecomposition(L).K );

   bcpt:= corelg.base( R, cpt );
   gcpt:= List( bcpt, x -> corelg.refl( R, x ) );

   hr:= Sum( re, i -> ch[1][i]*ch[2][i] );
   hi:= Sum( im, i -> ch[1][i]*ch[2][i] );
   ex:= Filtered( [1..Length(PositiveRoots(R))], i -> 
              IsZero(hr*ch[1][i]) and IsZero(hi*ch[1][i]) );

   if Length(ex) = 0 then 
      return [ gr, gi, ex, gcpt ];
   fi;

   bex:= corelg.base( R, ex );
   gex:= List( bex, x -> corelg.refl( R, x ) );

   h:= corelg.coroots(L,H);
   theta:= CartanDecomposition(L).CartanInv;
   p:= PermList( List(h,x -> Position(h,theta(x)) ) );

   gex0:= GeneratorsOfGroup( Centralizer( Group(gex), p ) );

   return [ gr, gi, gex0, gcpt ];

   # we have that the centralizer is generated by the first three lists,
   # the second may not be contained in the real weyl group; the fourth
   # list generates a subgroup of the second, that is contained in the
   # real Weyl group.

end;



############################################################################
##
#F   RealWeylGroup( <L> )
#F   RealWeylGroup( <L>, <H> )
##
##   <L> is a semisimple lie algebra over Gaussian rationals or SqrtField;
##   <H> is a Cartan subalgebra of <L>; returns the real Weyl group wrt <H>
##   if <H> is not given as input, then CartanSubalgebra(<L>) is used
##
InstallGlobalFunction(  RealWeylGroup, function( arg )

    local q, W0, G, cd, R, rank, ch, n, rvcs, rts, xx, theta, N, nu, lammat,
          i, j, k, l, u, cf, cfs, tr, Hm, P, rhs, pp, pos, elms, rhsj, ex, 
          good, g, tmp, L, H;


    L := arg[1];
    if Length(arg) = 2 then H := arg[2]; else H := CartanSubalgebra(L); fi;

    if HascorelgRealWG(H) then return corelgRealWG(H); fi;

    q:= corelg.centrzer(L,H);
    if Length(q[2]) = 0 then
       tmp := Group(Flat(q));
       SetcorelgRealWG(H,tmp);
       return tmp;
    fi; 
    W0:= Group(q[2]);
    if Length(q[4])=0 then
       G:= Group( [()] );
    else 
       G:= Group(q[4]); 
    fi;
    cd:= CosetDecomposition(W0,G);

    #Print("LENGTH CD:  ",Length(cd),"\n");

    R:= RootsystemOfCartanSubalgebra(L,H);
    rank:= Length( CartanMatrix(R) );
    ch:= ChevalleyBasis(R);
    n:= Length(ch[1]);
    rvcs:= Concatenation( ch[1], ch[2] );
    xx:= List( rvcs, u -> Basis( Subspace( L, [u] ), [u] ) );
    theta:= CartanDecomposition(L).CartanInv;
    rts:= Concatenation( PositiveRootsNF(R), -PositiveRootsNF(R) );

    N:= function( i , j )  # N_{alpha_i,alpha_j}

        local pos;

        pos:= Position( rts, rts[i]+rts[j] );
        if pos = fail then
           return Zero( LeftActingDomain(L) );
        fi;
        return Coefficients( xx[pos], rvcs[i]*rvcs[j])[1];
    end;

    nu:= function( w, i )

         local j, k, alpha, beta, gamma, pos, na, nb, sign, nn;

         if i <= n then
            j:= i; sign:= 1;
         else
            j:= i-n; sign:= -1;
         fi;
         if j <= rank then return One( LeftActingDomain(L) ); fi;

         gamma:= rts[j];
         for k in [1..rank] do
             pos:= Position( rts, gamma-rts[k] );
             if pos <> fail then
                alpha:= k; beta:= pos;
                break;
             fi;
         od;

         na:= nu( w, alpha );
         nb:= nu( w, beta );
         nn:= na*nb*N(alpha^w,beta^w)/N(alpha,beta);
         return nn^sign;

    end;          

    lammat:= [];
    for i in [1..rank] do
        u:= theta(ch[1][i]);
        for j in [1..Length(xx)] do
            cf:= Coefficients( xx[j], u );
            if cf <> fail then 
               cfs:= ShallowCopy( rts[j] );
               cfs[i]:= cfs[i]-1;
               Add( lammat, cfs );
               break;
            fi;
        od;
    od;

    tr:= HermiteNormalFormIntegerMatTransform(lammat);
    Hm:= tr.normal;
    P:= tr.rowtrans;

    rhs:= List( [1..rank], i -> One( LeftActingDomain(L) ) );
    pp:= [ ]; # contains positions of alpha_i^theta...
    for i in [1..rank] do
        u:= theta( ch[1][i] );
        for l in [1..Length(xx)] do
            cf:= Coefficients( xx[l], u );
            if cf <> fail then 
               rhs[i]:= rhs[i]*cf[1]^-1;
               pos:= l;
               break;
            fi;
        od;
        pp[i]:= pos;
    od;

    elms:= ShallowCopy( q[4] );
    for j in [2..Length(cd)] do
        rhsj:= ShallowCopy( rhs );
        for i in [1..rank] do
            k:= i^cd[j][1];
            u:= theta(rvcs[k]);
            for l in [1..Length(xx)] do
                cf:= Coefficients( xx[l], u );
                if cf <> fail then 
                   rhsj[i]:= rhsj[i]*cf[1];
                   break;
                fi;
            od;
            rhsj[i]:= rhsj[i]*nu( cd[j][1], pp[i] )^-1;
        od;

        good:= true;
        for i in [1..rank] do
            if IsZero(Hm[i]) then
               ex:= P[i];
               cf:= One( LeftActingDomain(L) );
               for l in [1..rank] do
                   cf:= cf*rhsj[l]^ex[l];
               od;
               if not IsOne(cf) then
                  good:= false;
                  break;
               fi;
            fi;
        od;
        if good then
           Add( elms, cd[j][1] );
        fi;
    od;
     
    Append( elms, q[1] ); Append( elms, q[3] );
    tmp := Group(elms);
    SetcorelgRealWG(H,tmp);
    return tmp;

end);


###################################Helper functions######################


corelg.setcd:= function(L,M) # M stable under theta of L...

   local cd;

   cd:= CartanDecomposition(L);
   SetCartanDecomposition( M, rec( CartanInv:= cd.CartanInv,
                                   K:= Intersection(cd.K,M), P:= Intersection(cd.P,M) ) );

end;

######################

corelg.cartdecsplit:= function(L)
local c1, c2, b1, b2, f, gr;
    # assume L is split
    if HasCartanDecomposition(L) then Error("a Cartan Decomposition is already set"); fi;

    c1:= CanonicalGenerators( RootSystem(L) );
    c2:= [ List(c1[2],x->-x),List(c1[1],x->-x),List(c1[3],x->-x) ];
    b1:= SLAfcts.canbas( L, c1 );
    b2:= SLAfcts.canbas( L, c2 );

    f:= AlgebraHomomorphismByImagesNC( L, L, Flat(b1), Flat(b2) );
    gr:= Grading(f);
    SetCartanDecomposition(L, rec(CartanInv:= function(u) 
                                     return Image(f,u); end,
                K:= Subalgebra(L,gr[1]), P:= Subspace(L,gr[2]) ) );

end;

######################

corelg.betterbasis:= function(L,K)

  local sigma, v, alpha, alphabar, a, b, k, l, eqns, eq, sol, bas, M;

  sigma:= function(u)

      local cf;

      cf:= Coefficients( Basis(L), u );
      cf:= List( cf, ComplexConjugate );
      return cf*Basis(L);
   end;

   v:= BasisVectors(Basis(K));
   alpha:= List( v, x -> Coefficients( Basis(K), sigma(x) ) );
   alphabar:= ComplexConjugate(alpha);

   a:= (alpha+alphabar)/2;
   b:= (alpha-alphabar)/(2*E(4)*One(SqrtField));

   eqns:= [ ];
   for k in [1..Dimension(K)] do
       eq:= [ ];
       for l in [1..Dimension(K)] do
           eq[l]:= a[l][k];
           eq[l+Dimension(K)]:= b[l][k];
       od;
       eq[k]:= eq[k]-1;
       Add( eqns, eq );
       eq:= [ ];
       for l in [1..Dimension(K)] do
           eq[l]:= b[l][k];
           eq[l+Dimension(K)]:= -a[l][k];
       od;
       eq[Dimension(K)+k]:= eq[Dimension(K)+k]-1;
       Add( eqns, eq );
   od;
   sol:= NullspaceMat( TransposedMat(eqns) );
   bas:= [ ];
   for v in sol do
       Add( bas, Sum( [1..Dimension(K)], k -> 
                                   (v[k]+E(4)*v[Dimension(K)+k])*Basis(K)[k]) );
   od;
   M:= Subalgebra(L,bas,"basis");
   corelg.setcd(L,M);

   return M;

end;

######################


corelg.IsSupport:= function( L, L0, c, e )

     # c: carrier alg, L: Lie alg, L0: zero comp..., e nilp elt in c_1

     local Cs, Ds, ds, t, Ca, Da, da, cr, c0, found, co, k, j, th, CCa;

     cr:= Subalgebra(L, Concatenation(c.g0,Flat(c.gp),Flat(c.gn)) );
     Cs:= Intersection( L0, LieCentralizer(L,cr) );
     Ds:= LieDerivedSubalgebra(Cs);
     corelg.setcd(L,Ds);
     ds:= Dimension( CartanSubspace(Ds) ) + Dimension( Intersection( LieCentre(Cs), CartanDecomposition(L).P ) );

     t:= SL2Triple( L, e ); 

     if not t[2] in L0 then
        Print("ERROR!!! h not in L0.\n");
     fi;

     # later make program that finds h in L0 - should be easy...

     Ca:= Intersection( L0, LieCentralizer(L,Subalgebra(L,t)) );
     Da:= LieDerivedSubalgebra(Ca);

     th:= CartanDecomposition(L).CartanInv;
     CCa:= LieCentre(Ca);
     if not ForAll( Basis(CCa), x -> th(x) in CCa ) then
        Print("ERROR!!! centre not theta-stable.\n");
     fi;

     da:= Dimension( CartanSubspace(Da) ) + Dimension( Intersection( CCa, CartanDecomposition(L).P ) );

Print("rank of Z_0(c): ",ds,"\nrank of Z_0(sl2): ",da,"\n");

     return ds=da;

end;

##
##
##  First part: functions for carrier algebras...
##
##  A carrier algebra is represented by a record with entries
##        g0: basis vectors of zero component
##        gp: list of lists, containing the basis vectors of the k-th component
##        gn: same, buth then -k-th component.
##
##  Main functions: 
##  
##         corelg.ZgradOrbs( L, grading )
##
##             Z-grading of a split semisimple Lie algebra, grading is a
##             list of integers, returns a record with components sl2s
##             and carr, the latter containing the carrier algebras.
##
##         corelg.CarrAlg( L, gr, t )
##
##             for a Z/mZ grading contained in gr, which is a list of lists,
##             on the k-th position a basis of the (k-1)-th component;
##             t is a homogeneous sl_2-triple. Returns a carrier algebra of t.
##
##
corelg.ZgradOrbs:= function( L, grading )

   # L: Lie algebra, gr: grading (0,1,-1 components).
   # 


   local R, B, ch, posR, N, rts, rr, pi, r1, zero, stack, res, r, 
         start, rrr, ips, i, vv, u, h, C, CT, pi_0, pi_1, t, s, pos,
         ct, eqns, rhs, eqn, j, sol, h0, psi0, psi1, good, x, y, es, fs, 
         valmat, val, chars, u0, v, done, gr1, gr2, g2, h_mats1, h_mats2, 
         mat, sl2s, id1, id2, Omega, V, e, ff, found, co, k, sp, extended,
         zz, bas, sim, Bw, W0, types, weights, wrts, tp, a, c, comb, hZ, hs,
         info, posRv, negRv, g0, g1, gm, CM, rr0, l0, l1, gr, deg, hs0, pis, pis0,
         cars, inds, K, gp, gn;


   ch:= ChevalleyBasis(L);

   R:= RootSystem(L);

   posR:= PositiveRootsNF(R);
   posRv:= PositiveRootVectors(R);
   negRv:= NegativeRootVectors(R);
   N:= Length( posR );
   rts:= ShallowCopy(posR);
   Append( rts, -posR );

   B:= BilinearFormMatNF(R);

   rr:= [ rec( pr0:= [ ], pv0:= [ ], nv0:= [] ), rec( r1:= [ ], rv1:= [ ] ), rec( rvm:= [ ] ) ];  
   for i in [1..Length(posR)] do
         v:= posR[i]*grading;
         if v = 0 then
            Add( rr[1].pr0, posR[i] );
            Add( rr[1].pv0, posRv[i] );
            Add( rr[1].nv0, negRv[i] );
         elif v = 1 then
            Add( rr[2].r1, posR[i] );
            Add( rr[2].rv1, posRv[i] );
            Add( rr[3].rvm, negRv[i] );
         fi;
   od;

   zz:= SLAfcts.zero_systems_Z( B, rr[1].pr0 );
   pi:= zz.subs;

   # now see how we can extend each element in pi with roots of
   # weight 1... and compute the maximal ones first!

   bas:= zz.bas;
   sim:= [ ];
   for a in bas do
       pos:= Position( posR, a );
       Add( sim, PositiveRootsAsWeights( R )[pos] );
   od;

   Bw:= SLAfcts.bilin_weights( R );
   W0:= rec( roots:= sim, wgts:= List( sim, x -> List( sim, y ->
                   2*x*(Bw*y)/( y*(Bw*y) ) ) ) );


   r1:= rr[2].r1;

   zero:= 0*r1[1];

   res:= [ ];
   for k in [1..Length(pi)] do

       types:= [ ];
       weights:= [ ];

       stack:= [ rec( rts0:= pi[k], rts1:= [ ], start:= 0,
                      sp:= VectorSpace( Rationals, pi[k], zero ) ) ];
       while Length(stack) > 0 do
           r:= stack[Length(stack)];
           RemoveElmList( stack, Length(stack) );
           start:= r.start+1;
           rrr:= Concatenation( r.rts0, r.rts1 );
           extended:= false;
           for i in [start..Length(r1)] do
               ips:= List( rrr, x -> x - r1[i] ); 
               if ForAll( ips, x -> not ( x in rts ) ) and
                           not r1[i] in r.sp then
                  vv:= ShallowCopy( BasisVectors( Basis(r.sp) ) );
                  Add( vv, r1[i] );
                  u:= ShallowCopy( r.rts1 );
                  Add( u, r1[i] );
                  Add( stack, rec( rts0:= r.rts0, rts1:= u, start:= i,
                          sp:= VectorSpace( Rationals, vv ) ) );
                  extended:= true;
               fi;
           od;
           if not extended then # see whether we can extend by
                                # adding something "smaller"
              for i in [1..start-1] do
                  if not r1[i] in rrr then
                     ips:= List( rrr, x -> x - r1[i] ); 
                     if ForAll( ips, x -> not ( x in rts ) ) and
                                    not r1[i] in r.sp then
                        extended:= true; break;
                     fi;
                  fi;
              od;
           fi;

           if not extended then 
              C:= List( rrr, x -> List( rrr, y -> 2*x*(B*y)/(y*(B*y)) ) );
              tp:= CartanType( C );
              SortParallel( tp.types, tp.enumeration );
              wrts:= [ ];
              for i in [1..Length(tp.enumeration)] do
                  for j in tp.enumeration[i] do
                      pos:= Position( rts, rrr[j] );
                      if pos <= N then
                         Add( wrts, PositiveRootsAsWeights(R)[pos] );
                      else
                         Add( wrts, -PositiveRootsAsWeights(R)[pos-N] );
                      fi;
                  od;
              od;
              found:= false;
              if tp.types in types then
                 for i in [1..Length(types)] do
                     if tp.types = types[i] then
                        if SLAfcts.my_are_conjugate( W0, R, Bw, wrts, weights[i] ) then
                           found:= true;
                           break;
                        fi;
                     fi;
                 od;
              fi;
              if not found then
                 Add( types, tp.types );
                 Add( weights, wrts );
                 Add( res, r );
              fi; 
           fi;
       od;

   od;

   stack:= [ ];
   for r in res do

       comb:= Combinations( [1..Length(r.rts1)] );
       comb:= Filtered( comb, x -> x <> [ ] );
       for c in comb do
           Add( stack, rec( rts0:= r.rts0, rts1:= r.rts1{c} ) );
       od;

   od;

   res:= stack;

info:= "Constructed ";
Append( info, String(Length(res)) );
Append( info, " root bases of possible flat subalgebras, now checking them...");
Info( InfoSLA, 2, info );

   h:= BasisVectors( Basis( CartanSubalgebra(L) ) );

   C:= CartanMatrix(R);
   CT:= TransposedMat( C );   

   good:= [ ];
   pis:= [ ];
   for r in res do

       pi_0:= r.rts0;
       pi_1:= r.rts1;
       pi:= Concatenation( pi_0, pi_1 );

       CM:= List( pi, x -> List( pi, y -> 2*x*(B*y)/( y*(B*y) ) ) );
       rr0:= SLAfcts.CartanMatrixToPositiveRoots( CM );
       l0:= 0; l1:= 0;
       gr:= Concatenation( List( pi_0, x -> 0 ), List( pi_1, x -> 1 ) );
       for s in rr0 do 
           deg:= s*gr;
           if deg=0 then
              l0:= l0+1;
           elif deg=1 then
              l1:= l1+1;
           fi;
       od;

       if 2*l0+Length(pi) = l1 then

          t:= [ ];
          for s in pi do
              pos:= Position( rts, s );
              if pos <= N then
                 Add( t, ch[1][pos]*ch[2][pos] );
              else
                 Add( t, ch[2][pos-N]*ch[1][pos-N] );
              fi;
          od; 

          t:= BasisVectors( Basis( Subspace( L, t ) ) );

          ct:= List( t, x -> Coefficients( Basis(CartanSubalgebra(L)), x ) );

          # i.e. t is a Cartan subalgebra of s

          # find h0 in t such that a(h0)=1 for all a in pi_1, a(h0)=0
          # for all a in pi_0

          eqns:=[ ];
          rhs:= [ ];
          for j in [1..Length(pi_0)] do
              eqn:= [ ];
              for i in [1..Length(t)] do
                  eqn[i]:= pi_0[j]*( C*ct[i] );
              od;
              Add( eqns, eqn ); Add( rhs, 0 );
          od;
          for j in [1..Length(pi_1)] do
              eqn:= [ ];
              for i in [1..Length(t)] do
                  eqn[i]:= pi_1[j]*( C*ct[i] );
              od;
              Add( eqns, eqn ); Add( rhs, 1 );
          od;

          sol:= SolutionMat( TransposedMat(eqns), rhs );
          h0:= sol*t;

          # Find a basis of the subspace of h consisting of u with 
          # a(u) = 0, for a in pi = pi_0 \cup pi_1.

          eqns:= [ ];
          for i in [1..Length(h)] do
              eqns[i]:= [ ];
              for j in [1..Length(pi_0)] do
                  Add( eqns[i], pi_0[j]*CT[i] );
              od;
              for j in [1..Length(pi_1)] do
                  Add( eqns[i], pi_1[j]*CT[i] );
              od;
          od;
          sol:= NullspaceMat( eqns );
          hZ:= List( sol, u -> u*h );

          # Now we compute |Psi_0| and |Psi_1|...

          psi0:= [ ];
          for a in rr[1].pv0 do 
              if h0*a = 0*a and ForAll( hZ, u -> u*a = 0*a ) then
                 Add( psi0, a );
              fi;
          od;

          psi1:= [ ];
          for a in rr[2].rv1 do
              if h0*a = a and ForAll( hZ, u -> u*a = 0*a ) then
                 Add( psi1, a );
              fi;
          od;

          if Length(pi_0)+Length(pi_1) + 2*Length(psi0) = Length(psi1) then

             if not 2*h0 in good then
                Add( good, 2*h0 );
                Add( pis, Concatenation(pi_0,pi_1) );
             fi;

          fi;
       fi;
   od;

info:= "Obtained ";
Append( info, String( Length(good) ) );
Append( info, " Cartan elements, weeding out equivalent copies...");
Info(InfoSLA,2,info);

# NEXT can be obtained from Kac diagram!!

   x:= ChevalleyBasis(L)[1];
   y:= ChevalleyBasis(L)[2];
   es:= [ ];
   fs:= [ ];
   g0:= Subspace( L, Concatenation( Basis(CartanSubalgebra(L)), rr[1].pv0, rr[1].nv0 ) );

   for i in [1..Length(CartanMatrix(R))] do
       if x[i] in g0 then
          Add( es, x[i] );
          Add( fs, y[i] );
       fi;
   od;
   hs:= List( [1..Length(es)], i -> es[i]*fs[i] );

   valmat:= [ ];
   for i in [1..Length(hs)] do
       val:= [ ];
       for j in [1..Length(hs)] do
           Add( val, Coefficients( Basis( Subspace(L,[es[j]]), [es[j]] ), 
                       hs[i]*es[j] )[1] );
       od;
       Add( valmat, val );
   od;


   chars:= [ ];
   hs0:= [ ];
   pis0:= [ ];
   for i in [1..Length(good)] do

       u0:= good[i];
       v:= List( es, z -> Coefficients( Basis(Subspace(L,[z]),[z]), u0*z )[1] );
       done:= ForAll( v, z -> z >= 0 );

       while not done do
           pos:= PositionProperty( v, z -> z < 0 );
           u0:= u0 - v[pos]*hs[pos];
           v:= v - v[pos]*valmat[pos];
           done:= ForAll( v, z -> z >= 0 );
       od;

       if not u0 in chars then
          Add( chars, u0 );
          Add( hs0, good[i] );
          Add( pis0, pis[i] );
       fi;
   od;

     sl2s:= [ ];
     cars:= [ ];
     Omega:= [-1,0,1,1];
     for i in [1..Length(hs0)] do

         # first we make the carrier...

         inds:= List( pis0[i], x -> Position( posR, x ) );
         K:= Subalgebra( L, Concatenation( posRv{inds}, negRv{inds} ) );
         mat:= List( Basis(K), x -> Coefficients( Basis(K), (hs0[i]/2)*x ) );
         g0:= List( NullspaceMat(mat), x -> x*Basis(K) );
         gp:= [ ]; gn:= [ ];
         k:= 1;
         while Length(g0) + Sum( gp, Length ) + Sum( gn, Length ) < Dimension(K) do
             gp[k]:= List( NullspaceMat( mat-k*mat^0 ), x -> x*Basis(K) ); 
             gn[k]:= List( NullspaceMat( mat+k*mat^0 ), x -> x*Basis(K) );
             k:= k+1; 
         od;
         Add( cars, rec( g0:= g0, gp:= gp, gn:= gn ) );

         # now get sl2 triple...
         found:= false;
         while not found do

             co:= List( gp[1], x -> Random(Omega) );
             x:= co*gp[1];
             sp:= Subspace( L, List( gn[1], y -> x*y) );

             if Dimension(sp) = Length(gp[1]) and hs0[i] in sp then

                # look for a nice one...
                for j in [1..Length(co)] do
                    k:= 0;
                    found:= false;
                    while not found do
                        co[j]:= k;
                        x:= co*gp[1];
                        sp:= Subspace( L, List( gn[1], y -> x*y) );

                        if Dimension(sp) = Length(gn[1]) and hs0[i] in sp then
                           found:= true;
                        else
                           k:= k+1;
                        fi;
                    od;
                od;

                mat:= List( gn[1], u -> Coefficients( Basis(sp), x*u ) );
                sol:= SolutionMat( mat, Coefficients( Basis(sp), hs0[i] ) );

                Add( sl2s, [sol*gn[1],hs0[i],x] );

                found:= true;

                       
             fi;
      
         od;
         
     od;
   
   return rec( sl2:= sl2s, carr:= cars );

end;



############################################



corelg.gradedSubalgByChar:= function( L, gr, h )

    # taken from corelg, modified a bit, for more general gradings...
   
    # here L is a Z/m-graded Lie algebra, grading in gr, m element list...
    # h nuetral elt of sl2 triple. We get the Z-graded subalgebra such that
    # g_k = { x\in L \mid x in gr[k mod m], [h,x] = 2*k*x}

    local adh, id, g0, g1, grad, gp, gn, k, done, cf, sp, m;

    m:= Length( gr );

    adh:= TransposedMat( AdjointMatrix( Basis(L), h ) );
    id:= adh^0;
    grad:= List( gr, u -> SubspaceNC( L, u, "basis" ) );
    gp:= [ ];
    k:= 1;
    done:= false;
    while not done do
       cf:= NullspaceMat( adh-2*k*id );
       if cf <> [] then
          sp:= Intersection( grad[(k mod m)+1], SubspaceNC( L, List( cf, c -> c*Basis(L) ) ) );
          Add( gp, BasisVectors( Basis(sp) ) );
          k:= k+1;
       else
          done:= true;
       fi;
    od;

    gn:= [ ];
    k:= 1;
    done:= false;
    while not done do
       cf:= NullspaceMat( adh+2*k*id );
       if cf <> [] then
          sp:= Intersection( grad[(-k mod m)+1], SubspaceNC( L, List( cf, c -> c*Basis(L) ) ) );
          Add( gn, BasisVectors( Basis(sp) ) );
          k:= k+1;
       else
          done:= true;
       fi;
    od;

    cf:= NullspaceMat( adh );
    sp:= Intersection( grad[1], SubspaceNC( L, List( cf, c -> c*Basis(L) ) ) );
    return rec( g0:= BasisVectors( Basis(sp) ), gp:= gp, gn:= gn );

end;


##################################################################################


corelg.carrierAlgBySL2:= function( L, H0, grad, sl2 )

   local R, B, ch, posR, N, rts, rr, pi, r1, zero, stack, res, r, 
         start, rrr, ips, i, vv, u, h, C, CT, pi_0, pi_1, t, s, pos,
         ct, eqns, rhs, eqn, j, sol, h0, psi0, psi1, good, x, y, es, fs, 
         valmat, val, chars, u0, v, done, gr1, gr2, g2, h_mats1, h_mats2, 
         mat, sl2s, id1, id2, Omega, V, e, ff, found, co, k, sp, extended,
         zz, bas, sim, Bw, W0, types, weights, wrts, tp, a, c, comb, hZ, hs,
         info, posRv, negRv, g0, g1, gm, CM, rr0, l0, l1, gr, deg, R0, gs, grading,
         cardat, U, gsp, grr, r0, gp, gn, L0, rvs, F, fct, rsp,H;

   # H0 is a Cartan subalgebra of the zero component, the carrier algebra
   # should be normalised by that one.

   gs:= corelg.gradedSubalgByChar( L, grad, sl2[2] );

   F:= LeftActingDomain(L);

   L0:= SubalgebraNC( L, Concatenation( gs.g0, Flat( gs.gp ), Flat( gs.gn ) ) );
   L0:= LieDerivedSubalgebra( L0 );
   gs.g0:= BasisVectors( Basis( Intersection( L0, SubspaceNC( L, gs.g0,"basis" ) ) ) );

   H:= Intersection(L0,H0);
   R0:= RootsystemOfCartanSubalgebra(L0,H);
   rvs:= Concatenation( PositiveRootVectors(R0), NegativeRootVectors(R0) );
   
   R0:= corelg.rtsys_withgrad( L0, rvs, H, gs );
   
   grading:= [ ];
   for v in CanonicalGenerators(R0)[1] do
       sp:= Basis( SubspaceNC( L, [v],"basis" ), [v] );
       Add( grading, Coefficients( sp, sl2[2]*v )[1]/2 );
   od;

   posR:= PositiveRootsNF(R0);
   posRv:= PositiveRootVectors(R0);
   negRv:= NegativeRootVectors(R0);
   N:= Length( posR );
   rts:= ShallowCopy(posR);
   Append( rts, -posR );

   B:= BilinearFormMatNF(R0);

   rr:= [ rec( pr0:= [ ], pv0:= [ ], nv0:= [] ), rec( r1:= [ ], rv1:= [ ] ), rec( rvm:= [ ] ) ];  
   for i in [1..Length(posR)] do
         v:= posR[i]*grading;
         if IsZero(v) then
            Add( rr[1].pr0, posR[i] );
            Add( rr[1].pv0, posRv[i] );
            Add( rr[1].nv0, negRv[i] );
         elif IsOne(v) then
            Add( rr[2].r1, posR[i] );
            Add( rr[2].rv1, posRv[i] );
            Add( rr[3].rvm, negRv[i] );
         fi;
   od;

   zz:= SLAfcts.zero_systems_Z( B, rr[1].pr0 );
   pi:= zz.subs;

   # now see how we can extend each element in pi with roots of
   # weight 1... and compute the maximal ones first!

   bas:= zz.bas;
   sim:= [ ];
   for a in bas do
       pos:= Position( posR, a );
       Add( sim, PositiveRootsAsWeights( R0 )[pos] );
   od;

   Bw:= SLAfcts.bilin_weights( R0 );
   W0:= rec( roots:= sim, wgts:= List( sim, x -> List( sim, y ->
                   2*x*(Bw*y)/( y*(Bw*y) ) ) ) );


   r1:= rr[2].r1;

   zero:= 0*r1[1];

   res:= [ ];
   for k in [1..Length(pi)] do

       types:= [ ];
       weights:= [ ];

       stack:= [ rec( rts0:= pi[k], rts1:= [ ], start:= 0,
                      sp:= VectorSpace( Rationals, pi[k], zero ) ) ];
       while Length(stack) > 0 do
           r   := stack[Length(stack)];
           rsp := BasisVectors(Basis(r.sp));
           if rsp = [] then 
              rsp := r.sp; 
           else 
              rsp := VectorSpace(Rationals,IdentityMat(Length(rsp[1]))); 
           fi;
           RemoveElmList( stack, Length(stack) );
           start:= r.start+1;
           rrr:= Concatenation( r.rts0, r.rts1 );
           extended:= false;
           for i in [start..Length(r1)] do
               ips:= List( rrr, x -> x - r1[i] ); 
               if ForAll( ips, x -> not ( x in rts ) ) and
                           not r1[i] in r.sp then
                  vv:= ShallowCopy( BasisVectors( Basis(r.sp) ) );
                  Add( vv, r1[i] );
                  u:= ShallowCopy( r.rts1 );
                  Add( u, r1[i] );
                  Add( stack, rec( rts0:= r.rts0, rts1:= u, start:= i,
                          sp:= VectorSpace( Rationals, vv ) ) );
                  extended:= true;
               fi;
           od;
           if not extended then # see whether we can extend by
                                # adding something "smaller"
              for i in [1..start-1] do
                  if not r1[i] in rrr then
                     ips:= List( rrr, x -> x - r1[i] ); 
                     if ForAll( ips, x -> not ( x in rts ) ) and not r1[i] in r.sp then
                          
                        extended:= true; break;
                     fi;
                  fi;
              od;
           fi;

           if not extended then 
              C:= List( rrr, x -> List( rrr, y -> 2*x*(B*y)/(y*(B*y)) ) );
              tp:= CartanType( C );
              SortParallel( tp.types, tp.enumeration );
              wrts:= [ ];
              for i in [1..Length(tp.enumeration)] do
                  for j in tp.enumeration[i] do
                      pos:= Position( rts, rrr[j] );
                      if pos <= N then
                         Add( wrts, PositiveRootsAsWeights(R0)[pos] );
                      else
                         Add( wrts, -PositiveRootsAsWeights(R0)[pos-N] );
                      fi;
                  od;
              od;
              found:= false;
              if tp.types in types then
                 for i in [1..Length(types)] do
                     if tp.types = types[i] then
                        if SLAfcts.my_are_conjugate( W0, R0, Bw, wrts, weights[i] ) then
                           found:= true;
                           break;
                        fi;
                     fi;
                 od;
              fi;
              if not found then
                 Add( types, tp.types );
                 Add( weights, wrts );
                 Add( res, r );
              fi; 
           fi;
       od;

   od;

   stack:= [ ];
   for r in res do

       comb:= Combinations( [1..Length(r.rts1)] );
       comb:= Filtered( comb, x -> x <> [ ] );
       for c in comb do
           Add( stack, rec( rts0:= r.rts0, rts1:= r.rts1{c} ) );
       od;

   od;

   res:= stack;

   C:= CartanMatrix(R0);
   CT:= TransposedMat( C );   

   sp:= Basis( SubspaceNC( L, CanonicalGenerators(R0)[3],"basis" ), CanonicalGenerators(R0)[3] );
   h:= BasisVectors( sp );

   good:= [ ];
   cardat:= [ ];
   for r in res do

       pi_0:= r.rts0;
       pi_1:= r.rts1;
       pi:= Concatenation( pi_0, pi_1 );

       CM:= List( pi, x -> List( pi, y -> 2*x*(B*y)/( y*(B*y) ) ) );
       rr0:= SLAfcts.CartanMatrixToPositiveRoots( CM );
       l0:= 0; l1:= 0;
       gr:= Concatenation( List( pi_0, x -> 0 ), List( pi_1, x -> 1 ) );
       for s in rr0 do 
           deg:= s*gr;
           if deg=0 then
              l0:= l0+1;
           elif deg=1 then
              l1:= l1+1;
           fi;
       od;

       if 2*l0+Length(pi) = l1 then

          t:= [ ];
          for s in pi do
              pos:= Position( rts, s );
              if pos <= N then
                 Add( t, posRv[pos]*negRv[pos] );
              else
                 Add( t, negRv[pos-N]*posRv[pos-N] );
              fi;
          od; 

          t:= BasisVectors( Basis( Subspace( L, t ) ) );

          ct:= List( t, x -> Coefficients( sp, x ) );

          # i.e. t is a Cartan subalgebra of s

          # find h0 in t such that a(h0)=1 for all a in pi_1, a(h0)=0
          # for all a in pi_0

          eqns:=[ ];
          rhs:= [ ];
          for j in [1..Length(pi_0)] do
              eqn:= [ ];
              for i in [1..Length(t)] do
                  eqn[i]:= pi_0[j]*( C*ct[i] );
              od;
              Add( eqns, eqn ); Add( rhs, Zero(F) );
          od;
          for j in [1..Length(pi_1)] do
              eqn:= [ ];
              for i in [1..Length(t)] do
                  eqn[i]:= pi_1[j]*( C*ct[i] );
              od;
              Add( eqns, eqn ); Add( rhs, One(F) );
          od;

          sol:= SolutionMat( TransposedMat(eqns), rhs );
          h0:= sol*t;

          # Find a basis of the subspace of h consisting of u with 
          # a(u) = 0, for a in pi = pi_0 \cup pi_1.

          eqns:= [ ];
          for i in [1..Length(h)] do
              eqns[i]:= [ ];
              for j in [1..Length(pi_0)] do
                  Add( eqns[i], pi_0[j]*CT[i] );
              od;
              for j in [1..Length(pi_1)] do
                  Add( eqns[i], pi_1[j]*CT[i] );
              od;
          od;
          sol:= NullspaceMat( eqns );
          hZ:= List( sol, u -> (u*One(F))*h );

          # Now we compute |Psi_0| and |Psi_1|...

          psi0:= [ ];
          for a in rr[1].pv0 do 
              if h0*a = 0*a and ForAll( hZ, u -> u*a = 0*a ) then
                 Add( psi0, a );
              fi;
          od;

          psi1:= [ ];
          for a in rr[2].rv1 do
              if h0*a = a and ForAll( hZ, u -> u*a = 0*a ) then
                 Add( psi1, a );
              fi;
          od;

          if Length(pi_0)+Length(pi_1) + 2*Length(psi0) = Length(psi1) then

             if not 2*h0 in good then
                Add( good, 2*h0 );
                Add( cardat, [ hZ, h0 ] );
             fi;

          fi;
       fi;
   od;

# NEXT can be obtained from Kac diagram!!

   x:= CanonicalGenerators(R0)[1];
   y:= CanonicalGenerators(R0)[2];
   es:= [ ];
   fs:= [ ];
   g0:= SubspaceNC( L, Concatenation( Basis(CartanSubalgebra(L)), rr[1].pv0, rr[1].nv0 ) );
 ##g0:= Subspace( L, Concatenation( Basis(MaximallyCompactCartanSubalgebra(L)), rr[1].pv0, rr[1].nv0 ) );

   for i in [1..Length(CartanMatrix(R0))] do
       if x[i] in g0 then
       #if corelg.eltInSubspace(L,BasisVectors(Basis(g0)),x[i]) then
          Add( es, x[i] );
          Add( fs, y[i] );
       fi;
   od;
   hs:= List( [1..Length(es)], i -> es[i]*fs[i] );

   valmat:= [ ];
   for i in [1..Length(hs)] do
       val:= [ ];
       for j in [1..Length(hs)] do
           Add( val, Coefficients( Basis( SubspaceNC(L,[es[j]]), [es[j]] ), 
                       hs[i]*es[j] )[1] );
       od;
       Add( valmat, val );
   od;


   chars:= [ ];
   fct:= function(x) if IsGaussRat(x) then return x; else return x![1][1][1]; fi; end;
   for i in [1..Length(good)] do

       u0:= good[i];
       v:= List( es, z -> Coefficients( Basis(SubspaceNC(L,[z]),[z]), u0*z )[1] );
       v:= List( v, fct );
       done:= ForAll( v, z -> z >= 0 );

       while not done do
           pos:= PositionProperty( v, z -> z < 0 );
           u0:= u0 - v[pos]*hs[pos];
           v:= v - v[pos]*valmat[pos];
           v:= List( v, fct );
           done:= ForAll( v, z -> z >= 0 );
       od;

       if not u0 in chars then
          Add( chars, u0 );
          if u0 = sl2[2] then
             U:= LieCentralizer( L, SubalgebraNC( L, cardat[i][1] ) );
             gsp:= List( grad, u -> SubspaceNC( L, u, "basis" ) );
             grr:= SL2Grading( L, cardat[i][2] );
             g0:= Intersection( U, gsp[1], SubspaceNC( L, grr[3] ) );
             g0:= SubalgebraNC( L, BasisVectors(Basis(g0)), "basis" );
             r0:= rec( g0:= BasisVectors( Basis( g0 ) ) );
             gp:= [ ];
             for j in [1..Length(grr[1])] do
                 g1:= Intersection( U, gsp[ (j mod Length(grad)) +1 ],SubspaceNC( L, grr[1][j]));
                 Add( gp, BasisVectors( Basis( g1 ) ) );
             od;
             gn:= [ ];
             for j in [1..Length(grr[2])] do
                 g1:= Intersection( U, gsp[(-j mod Length(grad)) +1 ],SubspaceNC( L, grr[2][j]));
                 Add( gn, BasisVectors( Basis( g1 ) ) );
             od;
             # remove trailing []-s...
             k:= Length(gp);
             while Length(gp[k]) = 0 do k:= k-1; od;
             gp:= gp{[1..k]};
             k:= Length(gn);
             while Length(gn[k]) = 0 do k:= k-1; od;
             gn:= gn{[1..k]};
             r0.gp:= gp; r0.gn:= gn;
             U:= SubalgebraNC( L, Concatenation( r0.g0, Flat(r0.gp), Flat(r0.gn) ), "basis" );
             U:= LieDerivedSubalgebra(U);
             r0.g0:= BasisVectors( Basis( Intersection( U, SubspaceNC( L, r0.g0, "basis" ) ) ) );

             r0.defelt := sl2[2]*(1/2);

             return r0;
          fi;
       fi;
   od;

   return "not found!!";

end;


############################################


corelg.CarrAlg0:= function( L, gr, sl2 )

   local h, lams, sp, i, gp, gn, eigensp, g0, g1, gm, m, K, k, dim,t0,e;

   e:= sl2[3];
   sp:= SubalgebraNC( L, gr[1] );
   sp:= Intersection( sp, LieCentralizer(L,SubalgebraNC(L,sl2)));
   if Dimension(sp) > 0 then
      h:= BasisVectors(CanonicalBasis(CartanSubalgebra(sp))); 
   else
      h:=[ ];
   fi;
   h:= Concatenation( [sl2[2]], h );
   lams:= [ ];
   sp:= BasisNC( SubspaceNC( L, [e],"basis" ), [e] );
   for i in [1..Length(h)] do
       Add( lams, Coefficients( sp, h[i]*e )[1] );
   od;

   gp:= [ ]; gn:= [ ];

    eigensp:= function( uu, t )

         local m, s, sp, eqns, i, j, k, c, sol;

         m:= Length(h);
         s:= Length(uu);
         sp:= Basis( SubspaceNC( L, uu ), uu );
         eqns:= NullMat( s, s*m );
         for j in [1..m] do
             for i in [1..s] do
                 c:= Coefficients( sp, h[j]*uu[i] );
                 for k in [1..s] do
                     eqns[i][(k-1)*m+j]:= c[k];
                 od;
             od;
         od;
         for k in [1..s] do
             for j in [1..m] do
                 eqns[k][(k-1)*m+j]:= eqns[k][(k-1)*m+j]-t*lams[j];
             od;
         od;

         sol:= NullspaceMat( eqns );
         return List( sol, x -> x*uu );
      end;

   m:= Length(gr);
   g0:= eigensp( gr[1], 0 );
   g1:= eigensp( gr[2], 1 );
   gm:= eigensp( gr[ m ], -1 );

   K:= LieDerivedSubalgebra( SubalgebraNC( L, Concatenation( gm, g0, g1 ) ) );

   g0:= BasisVectors( Basis( Intersection( SubspaceNC( L, g0,"basis" ), K ) ) );

   dim:= Length(g0);
   k:= 1;
   while dim < Dimension(K) do
      g1:= BasisVectors( Basis( Intersection( SubspaceNC( L, 
              eigensp( gr[ (k mod m) +1 ], k ) ), K ) ) );
      Add( gp, g1 );
      dim:= dim+Length(g1);
      gm:= BasisVectors( Basis( Intersection( SubspaceNC( L, 
              eigensp( gr[ (-k mod m) +1 ], -k ) ), K ) ) );
      Add( gn, gm );
      dim:= dim+Length(gm);
      k:= k+1;
   od;
 
   return rec( g0:= g0, gp:= gp, gn:= gn );
   
end;

############################################

corelg.CarrAlg:= function( L, H0, gr, sl2 )


   local cr, good, V, u, v;

   cr:= corelg.CarrAlg0( L, gr, sl2 );
   cr.defelt := sl2[2]*(1/2);
   good:= true;
   V:= Subspace( L, cr.g0, "basis" );
   for u in Basis(H0) do 
       for v in Basis(V) do
           if not u*v in V then
              good:= false; break;
           fi;
       od;
       if not good then break; fi;
   od;
   if good then
      V:= Subspace( L, cr.gp[1], "basis" );
      for u in Basis(H0) do 
          for v in Basis(V) do
              if not u*v in V then
                 good:= false; break;
              fi;
          od;
          if not good then break; fi;
      od;
   fi;

   if good then
      return cr;
   else
      return corelg.carrierAlgBySL2( L, H0, gr, sl2 );
   fi;

end;


##############################################################################
##
##
##
##     Second part: isomorphisms. 
##
##     Currently only for Z-gradings.
##     For Z/mZ gradings: map the characteristics, get sl2 triple, and
##     compute the carriers from there...
##
corelg.ZgradIsom:= function( L, H1, H2, grad )

    # MUST have f that respects grading!!!
    # For Z-gradings, take in both cases pos roots that have pos degree...

    # grad a grading in carrier form, ie a record with components g0, gp, gn.

    # can be much more efficient: can assume that H1 is the standard 
    # Cartan, so no work there, and running over the full symmetry group
    # must be avoided (although one can hope that it is not often necessary).

    local b1, b2, c1, c2, R1, R2, t, tp, en, i, d1, d2, g0, spc, C2, en0, rk, 
          sym, p, p0;

    R1:= RootsystemOfCartanSubalgebra(L,H1);
    R1:= corelg.rtsys_withgrad( L, Concatenation(PositiveRootVectors(R1),NegativeRootVectors(R1)),
                        H1, grad );
    R2:= RootsystemOfCartanSubalgebra(L,H2);
    R2:= corelg.rtsys_withgrad( L, Concatenation(PositiveRootVectors(R2),NegativeRootVectors(R2)),
                        H2, grad );

    g0:= Subspace( L, grad.g0, "basis" );
    spc:= Concatenation( [g0], List( grad.gp, u ->Subspace(L,u,"basis")) );

    t:= CartanType( CartanMatrix(R1) );
    tp:= ShallowCopy( t.types );
    en:= ShallowCopy( t.enumeration );
    SortParallel( tp, en );
    c1:= [ [], [], [] ];
    for i in [1..Length(en)] do
        Append( c1[1], CanonicalGenerators(R1)[1]{en[i]} );
        Append( c1[2], CanonicalGenerators(R1)[2]{en[i]} );
        Append( c1[3], CanonicalGenerators(R1)[3]{en[i]} );
    od;
    d1:= List( c1[1], x -> Filtered([0..Length(spc)-1], i -> 
            x in spc[i+1] ) );

    t:= CartanType( CartanMatrix(R2) );
    tp:= ShallowCopy( t.types );
    en:= ShallowCopy( t.enumeration );
    SortParallel( tp, en );
    c2:= [ [], [], [] ];
    for i in [1..Length(en)] do
        Append( c2[1], CanonicalGenerators(R2)[1]{en[i]} );
        Append( c2[2], CanonicalGenerators(R2)[2]{en[i]} );
        Append( c2[3], CanonicalGenerators(R2)[3]{en[i]} );
    od;

    d2:= List( c2[1], x -> Filtered([0..Length(spc)-1], i -> 
            x in spc[i+1] ) );

    if d1 <> d2 then # find permutation... (QD-way...)
Print("TRY TO FIND PERM -- better check!!!\n");
       C2:= CartanMatrix(R2);
       en0:= Flat(en);
       rk:= Length(en0);
       C2:= C2{en0}{en0};
       sym:= Elements( SymmetricGroup( rk ) );
       p0:= fail;
       for p in sym do
           if ForAll( [1..rk], i -> ForAll( [1..rk], j -> C2[i][j] = 
                   C2[i^p][j^p] ) ) then
              if ForAll( [1..rk], i -> d2[i^p] = d1[i] ) then
                 p0:= p;
                 break;
              fi;
           fi;
       od;
       if p0=fail then
          Print("NO perm found ERROR ERROR!!\n");
       fi;
       c2[1]:= List( [1..rk], i -> c2[1][i^p0] );
       c2[2]:= List( [1..rk], i -> c2[2][i^p0] );
       c2[3]:= List( [1..rk], i -> c2[3][i^p0] );
    fi;

    b1:= SLAfcts.canbas( L, c1 );
    b2:= SLAfcts.canbas( L, c2 );

    return AlgebraHomomorphismByImagesNC( L, L, Flat(b1), Flat(b2) );

end;

# THINGS THAT can go worng (and sometimes do):
# * in mapping the characteristic, there still is a diagram automorphism
#   possible, that does not leave the set of characteristics invariant
#   (of course, if the zero-comp does not have a diagram aut, then no prob),
#
# * in finding the carrier algebra, the characteristic found in the 
#   dominant Weyl chamber is not the one given (when using a different 
#   dominant Weyl chamber - have not seen this yet, but it may happen).
#
# So maybe compute the sl2-s anew for each CSA. (Maybe only in case there
# are diagram auts.)

corelg.mapsl2:= function( L, grad, L0, Z0, H1, H2, sl2 )

      # here L0 is the derived subalgebra of the zero-component,
      # Z0 is its centre, H1, H2 are CSA-s of the
      # zero component, so including Z0, and sl2s is a list of sl2 triples
      # wrt H1, get their images wrt H2.
      # grad is a list of bases of subspaces of L, giving the grading.

      local U1, U2, R1, R2, en, h1, h2, b, B, b2, gr1, m, grm, B1, Bm, 
            sl2s, t, h, adh, e, f, found, co, x, sp, i, k, mat, sol, tp,
            C2, rk, sym, p, good, goodperms, perms;
     
      U1:= Intersection( L0, H1 );
      U2:= Intersection( L0, H2 );
      R1:= RootsystemOfCartanSubalgebra( L0, U1 );
      R2:= RootsystemOfCartanSubalgebra( L0, U2 );
  
      tp:= CartanType( CartanMatrix(R1) );
      SortParallel( tp.types, tp.enumeration );

      en:= Flat( tp.enumeration );
      h1:= CanonicalGenerators( R1 )[3]{en};

      tp:= CartanType( CartanMatrix(R2) );
      SortParallel(tp.types, tp.enumeration );
      en:= Flat( tp.enumeration );

      h2:= CanonicalGenerators( R2 )[3]{en};

      b:= Concatenation( h1, Basis(Z0) );
      B:= Basis( Subspace( L, b ), b );

      b2:= Concatenation( h2, Basis(Z0) );
      gr1:= grad[2];
      m:= Length(grad);
      grm:= grad[m];
      B1:= Basis( Subspace( L, gr1 ), gr1 );
      Bm:= Basis( Subspace( L, grm ), grm );

       perms:= [ ];
       C2:= CartanMatrix(R2);
       rk:= Length(en);
       C2:= C2{en}{en};
       sym:= Elements( SymmetricGroup( rk ) );
       for p in sym do
           if ForAll( [1..rk], i -> ForAll( [1..rk], j -> C2[i][j] = 
                   C2[i^p][j^p] ) ) then
              Add( perms, p );
           fi;
       od;

       # now check for which permutation all characteristics work
       goodperms:= [ ];
       for p in perms do
           good:= true;
           h2:= CanonicalGenerators( R2 )[3]{en};
           h2:= List( [1..rk], i -> h2[i^p] );
           b2:= Concatenation( h2, Basis(Z0) );
           for t in sl2 do
               h:= Coefficients( B, t[2] )*b2;
               adh:= List( gr1, x -> Coefficients( B1, h*x ) );
               e:= List( NullspaceMat( adh-2*adh^0 ), u -> u*gr1 );
               adh:= List( grm, x -> Coefficients( Bm, h*x ) );
               f:= List( NullspaceMat( adh+2*adh^0 ), u -> u*grm );
               sp:= Subspace( L, Concatenation( List( e, x -> List( f, y -> x*y ) ) ) );
               if not h in sp then
                  good:= false;
                  break;
               fi;
           od;
           if good then Add( goodperms, p ); fi;
       od;

       if Length(goodperms) > 1 then 
          Print("MORE THAN ONE PERM POSSIBILE, taking the first one...\n");
       fi;
       p:= goodperms[1];

       if p <> () then Print("permutation used ",p,"\n"); fi;
       h2:= CanonicalGenerators( R2 )[3]{en};
       h2:= List( [1..rk], i -> h2[i^p] );
       b2:= Concatenation( h2, Basis(Z0) );
              
      sl2s:= [ ];
      for t in sl2 do
          h:= Coefficients( B, t[2] )*b2;
          adh:= List( gr1, x -> Coefficients( B1, h*x ) );
          e:= List( NullspaceMat( adh-2*adh^0 ), u -> u*gr1 );
          adh:= List( grm, x -> Coefficients( Bm, h*x ) );
          f:= List( NullspaceMat( adh+2*adh^0 ), u -> u*grm );

          found:= false;

          while not found do

                co:= List( e, x -> Random([-20..20]) );
                x:= co*e;
                sp:= Subspace( L, List( f, y -> x*y) );

                if h in sp then

                   # look for a nice one...
                   for i in [1..Length(co)] do
                       k:= 0;
                       found:= false;
                       while not found do
                           co[i]:= k;
                           x:= co*e;
                           sp:= Subspace( L, List( f, y -> x*y) );

                           if h in sp then
                              found:= true;
                           else
                              k:= k+1;
                           fi;
                       od;
                   od;

                   mat:= List( f, u -> Coefficients( Basis(sp), x*u ) );
                   sol:= SolutionMat( mat, Coefficients( Basis(sp), h ) );

                   Add( sl2s, [sol*f,h,x] );

                   found:= true;
                fi;
          od;
      od;
      return sl2s;

end;


##############################################################################
##
##
##    real Weyl group, and weight decomposition of the graded Lie algebra:
##    the i-th component is the sum of weight spaces (i,lambda) with
##    respect to the Cartan subalgebra H0 of g_0; we represent the
##    real Weyl group as a perm group on all weights.
##
##

corelg.weightvecdec:= function( F, h, vecs )

    # taken from RoootsysOFCSA....

    local B, i, j, newB, V, Mold, M, f, facs, facs0, num, fam, l, cf, b, c, r,
          one;

    one:= One(F);

    B:= [ vecs ];
    for i in h do
      newB := [ ];
      for j in B do

        if Length(j) = 1 then
           Add( newB, j ); 
        else
           V    := Basis( VectorSpace( F, j, "basis" ), j );
           Mold := List( j, x -> Coefficients( V, i*x ) );

           if fail in Flat(Mold) then
              Print("Extension of base field!, have to return fail\n");
              return fail;
           fi;

           if IsSqrtField(F) then
              M    := SqrtFieldMakeRational(Mold);
              if M = false then 
                 Print(" matrix cannot be made rations\n");
                 M    := Mold;
                 f    := CharacteristicPolynomial( M );
                 facs := Set(Factors( f ));
              else
                 f    := CharacteristicPolynomial( M );
                 facs := Set(Factors( f ));
                 f    := SqrtFieldRationalPolynomialToSqrtFieldPolynomial(f);
                 facs := Set(List(facs,SqrtFieldRationalPolynomialToSqrtFieldPolynomial));
              fi;
           else
              M    := Mold;
              f    := CharacteristicPolynomial( M );
              facs := Set(Factors( f ));
           fi;

           num  := IndeterminateNumberOfUnivariateLaurentPolynomial(f);
           fam  := FamilyObj( f );

           facs0:= [ ];

           for l in facs do
               if Degree(l) = 1 then
                  Add( facs0, l );
               elif Degree(l) = 2 then # we just take square roots...
                  cf := CoefficientsOfUnivariatePolynomial(l);
                  b  := cf[2];
                  c  := cf[1];
                  r  := (-b+Sqrt(b^2-4*c))/2;  
                  if not r in F then 
                     Print("cannot do this over ",F,"\n"); 
                     return fail;
                  fi;
                  Add( facs0, PolynomialByExtRep( fam, [ [], -r, [num,1], one] ) );
                  r  := (-b-Sqrt(b^2-4*c))/2;
                  if not r in F then 
                     Print("cannot do this over ",F,"\n"); 
                     return fail;
                  fi;
                  Add( facs0, PolynomialByExtRep( fam, [ [], -r, [num,1], one] ) );

               else
                  Print("not split\n");
                  return fail;
               fi;
           od;

           for l in facs0 do
              V := NullspaceMat( Value( l, Mold ) );
              Add( newB, List( V, x -> LinearCombination( j, x ) ) );
           od;
        fi;

      od;
      B:= newB;
   od;

   return B;


end;


####################################################

corelg.RealWeylGroupWts:= function( L, grad, m0, L0, Z0, H0 )

       # if m = infinity then grad is a record of carrier type (g0, gp, gn),
       # otherwise it is a list of lists...
       # L0 is the derived subalgebra of the zero-component, Z0 is its centre,
       # H0 is a Cartan subalgebra of the zero-component.

       local U0, W, gW, R, h, rk, ms, e, sp, vv, ch, N, g, m, i, j, z, b, v, 
             wts, wvecs, perms, p, w0, w, alpha, cfwts, gWC, msC, pC;

       U0:= Intersection( L0, H0 );
       W:= RealWeylGroup( L0, U0 );
Print("size real Weyl group:  ",Size(W),"\n");
       gW:= GeneratorsOfGroup(W);
       R:= RootsystemOfCartanSubalgebra( L0, U0 );

       gWC:= SLAfcts.perms(R); # also the big Weyl group....
       
       h:= Concatenation( CanonicalGenerators(R)[3], Basis(Z0) );
       # we represent a weight as a vector of values rel to h
       rk:= Length( CartanMatrix(R) );
       # so the first rk components of such a vector belong to the can gens,
       # note that the real Weyl group does not act on the remaining
       # coordinates.
       # first we get matrix reps of the generators on the first rk
       # coordinates...

       ms:= [ ];
       e:= CanonicalGenerators( R )[1];
       sp:= List( e, x -> Basis( Subspace(L,[x]),[x]) );
       alpha:= List( [1..rk], i -> List( h{[1..rk]}, u -> Coefficients( sp[i],
                   u*e[i] )[1] ) );
       sp:= Basis( VectorSpace( LeftActingDomain(L), alpha ), alpha );
       ch:= ChevalleyBasis(R);
       N:= Length(ch[1]);
       for g in gW do
           m:= [ ];
           for i in [1..rk] do
               j:= i^g;
               if j > N then
                  z:= ch[2][j-N];
               else
                  z:= ch[1][j];
               fi;
               b:= Basis( Subspace( L, [z] ), [z] );
               v:= List( h{[1..rk]}, u -> Coefficients( b, u*z )[1] );
               Add( m, Coefficients( sp, v ) );
           od;
           Add( ms, m );
       od;

       msC:= [ ];
       for g in gWC do
           m:= [ ];
           for i in [1..rk] do
               j:= i^g;
               if j > N then
                  z:= ch[2][j-N];
               else
                  z:= ch[1][j];
               fi;
               b:= Basis( Subspace( L, [z] ), [z] );
               v:= List( h{[1..rk]}, u -> Coefficients( b, u*z )[1] );
               Add( m, Coefficients( sp, v ) );
           od;
           Add( msC, m );
       od;
       
       # now compute all weights, the weights (0,lambda) are the 
       # roots in L0
       # we leave out the weights of the form (i,0), as they may have
       # multiplicity >1, and they will never occur in a carrier algebra.

       wts:= [ ];
       wvecs:= [ ];
       for z in ch[1] do 
           b:= Basis( Subspace( L, [z] ), [z] );
           v:= List( h, u -> Coefficients( b, u*z )[1] );      
           Add( wts, [0,v] );
           Add( wvecs, z );
       od;
       for z in ch[2] do 
           b:= Basis( Subspace( L, [z] ), [z] );
           v:= List( h, u -> Coefficients( b, u*z )[1] );      
           Add( wts, [0,v] );
           Add( wvecs, z );
       od;

       if m0 = infinity then
  
          for i in [1..Length(grad.gp)] do

              vv:= corelg.weightvecdec( LeftActingDomain(L), h, grad.gp[i] );
              for z in vv do
                  z:= z[1];
                  b:= Basis( Subspace( L, [z] ), [z] );
                  v:= List( h, u -> Coefficients( b, u*z )[1] );      
                  Add( wts, [i,v] );
                  Add( wvecs, z );
              od;

          od;

          for i in [1..Length(grad.gn)] do

              vv:= corelg.weightvecdec( LeftActingDomain(L), h, grad.gn[i] );
              for z in vv do
                  z:= z[1];
                  b:= Basis( Subspace( L, [z] ), [z] );
                  v:= List( h, u -> Coefficients( b, u*z )[1] );      
                  Add( wts, [i,v] );
                  Add( wvecs, z );
              od;

          od;

       else

          for i in [2..Length(grad)] do

              vv:= corelg.weightvecdec( LeftActingDomain(L), h, grad[i] );
              for z in vv do
                  z:= z[1];
                  b:= Basis( Subspace( L, [z] ), [z] );
                  v:= List( h, u -> Coefficients( b, u*z )[1] );      
                  if not IsZero(v) then
                     Add( wts, [i-1,v] );
                     Add( wvecs, z );
                  fi;
              od;

          od;

       fi;

       perms:= [ ];
       cfwts:= List( wts, w -> Coefficients( sp, w[2]{[1..rk]}) );
       for m in ms do
           p:= [ ];
           for i in [1..Length(wts)] do
               w0:= (cfwts[i]*m)*alpha;
               Append( w0, wts[i][2]{[rk+1..Length(wts[i][2])]} );
               w:= [ wts[i][1], w0 ];
               Add( p, Position( wts, w ) );
           od;
           Add( perms, PermList( p ) );
       od;

       pC:= [ ];
       for m in msC do
           p:= [ ];
           for i in [1..Length(wts)] do
               w0:= (cfwts[i]*m)*alpha;
               Append( w0, wts[i][2]{[rk+1..Length(wts[i][2])]} );
               w:= [ wts[i][1], w0 ];
               Add( p, Position( wts, w ) );
           od;
           Add( pC, PermList( p ) );
       od;

       return rec( WR:= Group(perms), WC:= Group(pC), wts:= wts, wvecs:= wvecs );

end; 


#############################################################################
##
##  Main functions....
##


corelg.realcarriers:= function( L, grad, m0, g0, H1, H2, data )

        # m0 is an integer or infinity,
        # grad is as in previous function
        # g0 is the zero-component of L, as a subalgebra,
        # H1, H2 are Cartan subalgebras of g0
        # H1 is the "standard" or "split" one with respect to which 
        # the nilpotent orbits are already computed,
        # if m < infinity then data is a list of sl2-s,
        # if m = infinity then data is a list of carrier algebras,
        # both wrt H1.

        local sig, L0, Z0, cars, f, c, r, sl2, V, i, j, cos, elms, WR, res,
              wv, P, ncptH2, realf, g, gg0, ggp, ggn, U, N, DN, csp, found,
              p, u0, u1, s, weyldata, HN, M0, Mm, M1, S;

    sig:= function(u)
       return List( Coefficients( Basis(L), u ), ComplexConjugate )*Basis(L); 
    end;

        L0:= LieDerivedSubalgebra(g0);
        Z0:= LieCentre(g0);

        # first we get all complex carrier algebras wrt H2...
        cars:= [ ];
        if m0 = infinity then
           
           f:= corelg.ZgradIsom( L, H1, H2, grad );
           for c in data do
               r:= rec( g0:= List( c.g0, x -> Image(f,x) ) );
               r.gp:= List( c.gp, u -> List( u, x -> Image(f,x) ) );
               r.gn:= List( c.gn, u -> List( u, x -> Image(f,x) ) );
               Add( cars, r );
           od;

        else

           sl2:= corelg.mapsl2( L, grad, L0, Z0, H1, H2, data );

           for c in sl2 do 
               Add( cars, corelg.CarrAlg( L, H2, grad, c ) );
           od;

        fi;

        # get Weyl groups, and weight dec
        weyldata:= corelg.RealWeylGroupWts( L, grad, m0, L0, Z0, H2 );

        # in each carrier algebra we add components g0w, gpw, gnw,
        # consisting of lists of integers, indicating the index of 
        # the weight vecs occurring in the corresponding components -
        # makes computing the action easy...

        for c in cars do

            c.g0w:= [ ];
            V:= Subspace( L, c.g0 );
            for i in [1..Length(weyldata.wvecs)] do
                if weyldata.wvecs[i] in V then
                   Add( c.g0w, i );
                fi;
            od;

            c.gpw:= List( c.gp, x -> [] );
            for j in  [1..Length(c.gp)] do
                V:= Subspace( L, c.gp[j] );
                for i in [1..Length(weyldata.wvecs)] do
                    if weyldata.wvecs[i] in V then
                       Add( c.gpw[j], i );
                    fi;
                 od;
                 if not Length( c.gpw[j] ) = Length( c.gp[j] ) then
                    Print("grr1\n");
                 fi;
            od;

            c.gnw:= List( c.gn, x -> [] );
            for j in  [1..Length(c.gn)] do
                V:= Subspace( L, c.gn[j] );
                for i in [1..Length(weyldata.wvecs)] do
                    if weyldata.wvecs[i] in V then
                       Add( c.gnw[j], i );
                    fi;
                 od;
                 if not Length( c.gnw[j] ) = Length( c.gn[j] ) then
                    Print("grr2\n");
                 fi;
            od;
        od;

        cos:= CosetDecomposition( weyldata.WC, weyldata.WR );     
        elms:= List( cos, x -> x[1]^-1 );    
        WR:= Elements( weyldata.WR );

        res:= [ ];
        wv:= weyldata.wvecs;
        P:= CartanDecomposition(L).P;
        ncptH2:= Dimension( Intersection( H2, P ) );

        for c in cars do
            # Find all real carrier algebras that are WC-conjugate to c... 

            realf:= [ ];
            for g in elms do
                gg0:= List( c.g0w, i -> i^g );
                ggp:= List( c.gpw, u -> List( u, i -> i^g ) );
                ggn:= List( c.gnw, u -> List( u, i -> i^g ) );

                M1:= Subspace(L,wv{ggp[1]});
                Mm:= Subspace(L,wv{ggn[1]});
                M0:= Subspace(L,wv{gg0} );
                U:= Subalgebra( L, Concatenation( 
                         wv{gg0}, wv{Flat(ggp)}, wv{Flat(ggn)} ) );
                if ForAll( Basis(M1), x -> sig(x) in M1 ) and
                   ForAll( Basis(Mm), x -> sig(x) in Mm ) and
                   ForAll( Basis(M0), x -> sig(x) in M0 ) then

                   # so it is a real subalgebra - just need to check the
                   # generators!

                   # check whether it is strongly H2-regular...

                   U:= corelg.betterbasis(L,U);

                   N:= Intersection( g0, corelg.normalizer( L, U ) );
                   DN:= LieDerivedSubalgebra(N);
                   corelg.setcd(L,DN);
                   csp:= CartanSubspace(DN);

#Print("csp... ",Dimension(csp),"\n");
#HN:= MaximallyNonCompactCartanSubalgebra(DN);
#HN:= Subalgebra( L, Concatenation( Basis(HN), Basis(LieCentre(N)) ) );
#Print(Dimension(Intersection(HN,P)),"  ",Dimension(Intersection(HN,
#CartanDecomposition(L).K)),"  ",
#Dimension(csp)+Dimension(Intersection(LieCentre(N),P)),"  ",
#ncptH2,"\n");
#Print(gg0,"\n",ggp,"\n",ggn,"\n");



                   if Dimension(csp)+Dimension(Intersection(LieCentre(N),P))=
                      ncptH2 then  #yes!

                      # see whether it is conjugate to something seen before

                      found:= false;
                      for p in WR do
                          u0:= Set( gg0, i -> i^p );
                          u1:= Set( ggp[1], i -> i^p );
                          for s in realf do
                              if Set(s.g0w) = u0 and Set(s.gpw[1]) = u1 then
                                 found:= true; break;
                              fi;
                          od;
                          if found then break; fi;
                      od;
                      if not found then

                          r:= rec(g0:= BasisVectors(Basis(M0)));
                          r.gp:= List( ggp, u -> BasisVectors( 
                       CanonicalBasis(Intersection(U,Subspace(L,wv{u})))));
                          r.gn:= List( ggn, u -> BasisVectors( 
                       CanonicalBasis(Intersection(U,Subspace(L,wv{u})))));
                          r.g0w:= gg0; 
                          r.gpw:= ggp;
                          r.gnw:= ggn;
                          
                          S:= Subalgebra( L, Concatenation( Concatenation( 
                                List( r.gp[1], x -> List( r.gn[1], y -> x*y))),
                                Concatenation( List( r.g0, x -> List( r.g0,
                                 y-> x*y)))));
                          r.g0:= BasisVectors(CanonicalBasis(
                                                        Intersection(U,S)));
                          Add( realf, r );
                      fi;
                   fi;
                fi;
            od;
            Add( res, realf );
        od;           
        return res;

end; 





############################################################################
##
#F   CarrierAlgsForNilpOrbsInZmGrading( <type>, <rank>, <m0>, <str>, <num> )
##
#   Gives a record containing the carrier algebras of the real theta group
#   specified by the input. Explanation of the input:
#
#   type: type of the Lie algebra where everything happens,
#   rank: its rank,
#   m0 : the order of the automorphism defining the grading,
#   str: "inner" or "outer", the first when an inner automorphism
#        defines the grading, the second otherwise,
#   num : the num-th automorphism in the list 
#                 FiniteOrderInnerAutomorpisms( type, rank, m0 ),
#         or 
#                 FiniteOrderOuterAutomorphisms( type, rank, m0, 2 ) 
#   is used to define the grading.
#
#   The output is a record with the following fields:
# 
#   L : the Lie algebra,
#   grad : the grading that was used
#   Hs : the Cartan subalgebras of g_0 that are used,
#   L0 : subalgebra g_0,
#   cars : the carrier algebras, this is a list of lists; for each Cartan 
#          subalgebra of g_0 there is one list: the first corresponds to the 
#          split Cartan subalgebra, and so has just the complex carrier 
#          algebras (which are also real), the other lists contain lists as 
#          well, for each complex carrier algebra (i.e., for each element of 
#          the first list) there is a list containing the real carrier 
#          algebras which are strongly H_i-regular, and over the complexes
#          conjugate to the given complex carrier algebra.
#          Furthermore, a carrier algebra is given by a record, containing 
#          the fields g0, gp (positive degree), and gn (negative degree). 

InstallGlobalFunction(CarrierAlgsForNilpOrbsInZmGrading, function( type, rank, m0, isin, ind )

       # m0 is the order of the automorphism,
       # isin = "inner" or isin = "outer",
       # ind: the ind-th such automorphism in the SLA listing.
   
       local L, gr, K, one, cfs, grad, t, sl2, L0, D0, Z0, Hs, cars, i,
             f, tmp;

       if isin = "inner" then
          f:= FiniteOrderInnerAutomorphisms( type, rank, m0 )[ind];
       else
          f:= FiniteOrderOuterAutomorphisms( type, rank, m0, 2 )[ind];
       fi;

       L:= SimpleLieAlgebra(type,rank,SqrtField);
       corelg.cartdecsplit(L);

       gr:= Grading(f);
       K:= Source(f);
       one:= One( SqrtField );
       cfs:= List( gr, u -> List( u, x -> Coefficients( Basis(K), x )*one ) );
       grad:= List( cfs, u -> List( u, x -> x*Basis(L) ) );

       t:= NilpotentOrbitsOfThetaRepresentation(f);
       cfs:= List( t, u -> List( u, x -> Coefficients( Basis(K), x )*one ) );
       sl2:= List( cfs, u -> List( u, x -> x*Basis(L) ) );


       L0:= Subalgebra( L, grad[1] );
       D0:= LieDerivedSubalgebra( L0 );
       corelg.setcd(L,D0);
       Z0:= LieCentre( L0 );
       Hs:= CartanSubalgebrasOfRealForm(D0 );
       Hs:= List( Hs, U -> Subalgebra( L0, Concatenation(Basis(U),Basis(Z0))));

       cars:= [ ];
       # WE ASSUME that the first CSA is the split one!
       # so no work there..

       Add( cars, List( sl2, u -> corelg.CarrAlg( L, Hs[1], grad, u ) ) );
       for i in [2..Length(Hs)] do
           Add( cars, corelg.realcarriers( L, grad, m0, L0, Hs[1], Hs[i], sl2 ) );
       od;

       return rec(L:=L,cars:=cars,grad:= grad,Hs:=Hs,L0:=L0);

end);



##################################################################
############################################################################
##
#F   CarrierAlgsForNilpOrbsInZGrading( <type>, <rank>, <d> )
##
#   Gives a record containing the carrier algebras of the real theta group
#   specified by the input. Explanation of the input:
#
#   type: type of the Lie algebra where everything happens,
#   rank: its rank,
#   d : (for Z grading) the degrees of the simple roots,
#  
#
#   The output is a record with the following fields:
# 
#
#   L : the Lie algebra,
#   grad : the grading that was used 
#   Hs : the Cartan subalgebras of g_0 that are used,
#   L0 : subalgebra g_0,
#   cars : the carrier algebras, this is a list of lists; for each Cartan 
#          subalgebra of g_0 there is one list: the first corresponds to the 
#          split Cartan subalgebra, and so has just the complex carrier 
#          algebras (which are also real), the other lists contain lists as 
#          well, for each complex carrier algebra (i.e., for each element of 
#          the first list) there is a list containing the real carrier 
#          algebras which are strongly H_i-regular, and over the complexes
#          conjugate to the given complex carrier algebra.
#          Furthermore, a carrier algebra is given by a record, containing 
#          the fields g0, gp (positive degree), and gn (negative degree). 

InstallGlobalFunction( CarrierAlgsForNilpOrbsInZGrading,  function( type, rank, d )

   
       local L, gr, R, C, h, grad, LL0, r, toSqrtField, car, i, c, L0, D0, 
             Z0, Hs, cars;

       L:= SimpleLieAlgebra(type,rank,SqrtField);
       corelg.cartdecsplit(L);

       R:= RootSystem(L);
       C:= CartanMatrix(R);
       h:= (C^-1*d*One(LeftActingDomain(L)))*ChevalleyBasis(L)[3];
       gr:= SL2Grading( L, h );  # INEFFICIENT!
       grad:= rec( g0:= gr[3], gp:= gr[1], gn:= gr[2] );

       LL0:= SimpleLieAlgebra( type, rank, Rationals );
       r:=  corelg.ZgradOrbs( LL0, d );

       toSqrtField:= function( lst )
          return List( lst, x -> 
              (Coefficients(Basis(LL0),x)*One(SqrtField))*Basis(L) );
       end;

       car:= [ ];
       for i in [1..Length(r.carr)] do       
           c:= rec( g0:= toSqrtField(r.carr[i].g0),
                    gp:= List( r.carr[i].gp, toSqrtField ),
                    gn:= List( r.carr[i].gn, toSqrtField ) );
           Add( car, c );
       od;         

       L0:= Subalgebra( L, grad.g0 );
       D0:= LieDerivedSubalgebra( L0 );
       corelg.setcd(L,D0);
       Z0:= LieCentre( L0 );
       Hs:= CartanSubalgebrasOfRealForm(D0 );
       Hs:= List( Hs, U -> Subalgebra( L0, Concatenation(Basis(U),Basis(Z0))));

       cars:= [ ];
       # WE ASSUME that the first CSA is the split one!
       # so no work there..

       Add( cars, car );
       for i in [2..Length(Hs)] do
           Add( cars, corelg.realcarriers( L, grad, infinity, L0, Hs[1], Hs[i], car));
       od;

       return rec(L:=L,cars:=cars,grad:=grad,Hs:=Hs,L0:=L0);

end);



##########################################


corelg.finde:= function( L, L0, c )


      # c: carrier alg, L the Lie alg; find an sl2-triple
      # maybe c.g0 is incomplete!
      # L0: zero comp of L

      local C, C0, C1, n, Omega, found, cf, e, k, j;

      #C:= Subalgebra( L, Concatenation( Flat(c.gp), Flat(c.gn) ) );
      #C0:= Intersection( L0, C );
      C0:= Subalgebra( L, Concatenation( List( c.gp[1], x -> x*c.gn[1] ) ) ); 
      C1:= BasisVectors( CanonicalBasis( Subspace( L, c.gp[1] ) ) );
      n:= Length(C1);
      Omega:= [-Dimension(L)..Dimension(L)];
#Omega:= [-2..2];
      found:= false;
      while not found do 
            cf:= List( C1, x -> Random(Omega) );
            e:= cf*C1;
            if Dimension( Subspace( L, Basis(C0)*e ) ) = n then
               found:= true;
            fi;
      od;
#return e;
      for k in [1..Length(C1)] do
          j:= 0;
          found:= false;
          while not found do
               cf[k]:= j;
               e:= cf*C1;
               if Dimension( Subspace( L, Basis(C0)*e ) ) = n then
                  found:= true;
               else
                  j:= j+1;
               fi;
          od;
      od;
    
      return e;

end;



#############################################################################

corelg.issupp:= function( L, L0, c, e )


     local t, C;

 
     t:= SL2Triple(L,e);
     if not t[2] in L0 then
        Print("ERROR!!!! h not in L0.\n");
     fi;
     
     C:= Intersection( L0, LieCentralizer( L, Subalgebra(L,t) ) );
     

end;

#############################################################################

corelg.rankC:= function( L, L0, c )

    local C, A, D;

    C:= Subalgebra(L,Concatenation(Flat(c.gp),Flat(c.gn)) );
    A:= Intersection(L0,LieCentralizer(L,C));
    D:= LieDerivedSubalgebra(A);
    corelg.setcd(L,D);
    return Dimension(CartanSubspace(D))+Dimension(
       Intersection( CartanDecomposition(L).P, LieCentre(A)) ); 

end;

#############################################################################


corelg.rankT:= function( L, L0, c )

    local e, B, DB, d, t, theta;

    e:= corelg.finde(L,L0,c); t:= SL2Triple(L,e); 
    Print(t[1] in Subspace(L,c.gn[1]),"\n");
    B:= Intersection(L0,LieCentralizer(L,Subalgebra(L,t)));
    if Dimension(B) = 0 then return 0; fi;
    B:= corelg.betterbasis(L,B);
    DB:= LieDerivedSubalgebra(B);
    if Dimension(DB)>0 then
       d:= Dimension(CartanSubspace(DB));
    else
       d:= 0;
    fi;
    theta:= CartanDecomposition(L).CartanInv;
    Print(ForAll( Basis(LieCentre(B)), x -> theta(x) in LieCentre(B) ),"\n");
    d:= d+Dimension(Intersection( CartanDecomposition(L).P, LieCentre(B) ));
    return d;

end;

#############################################################################


corelg.rankT_e:= function( L, L0, c, e )

    local B, DB, d, t, theta;

    t:= SL2Triple(L,e); 
    Print(t[1] in Subspace(L,c.gn[1]),"\n");
    B:= Intersection(L0,LieCentralizer(L,Subalgebra(L,t)));
    if Dimension(B) = 0 then return 0; fi;
    B:= corelg.betterbasis(L,B);
    DB:= LieDerivedSubalgebra(B);
    if Dimension(DB)>0 then
       d:= Dimension(CartanSubspace(DB));
    else
       d:= 0;
    fi;
    theta:= CartanDecomposition(L).CartanInv;
    Print(ForAll( Basis(LieCentre(B)), x -> theta(x) in LieCentre(B) ),"\n");
    d:= d+Dimension(Intersection( CartanDecomposition(L).P, LieCentre(B) ));
    return d;

end;

#############################################################################


corelg.desZT:= function( L, L0, c )

    local e, B, DB, d, theta, t;

    e:= corelg.finde(L,L0,c); t:= SL2Triple(L,e); 
    Print(t[1] in Subspace(L,c.gn[1]),"\n");
    B:= Intersection(L0,LieCentralizer(L,Subalgebra(L,t)));
    B:= corelg.betterbasis(L,B);
    DB:= LieDerivedSubalgebra(B);
    if Dimension(DB)>0 then
       Display( SatakeDiagram(DB) );
    fi;
    theta:= CartanDecomposition(L).CartanInv;
    Print(ForAll( Basis(LieCentre(B)), x -> theta(x) in LieCentre(B) ),"\n");
    Print( Dimension(Intersection( CartanDecomposition(L).K, LieCentre(B) )),
  "  ",Dimension(Intersection( CartanDecomposition(L).P, LieCentre(B) )),"\n");

end;

#############################################################################

corelg.diagmat:= function( L, H ) # H torus, find a matrix C such that 
                           # C^-1*ad h*C is diagonal for all h in H


    local spaces, i, h, sp0, s, sp, A, fct, f, V, bas, C, F, one,
          fct0, num, fam, cf, b, c, r;

    F:= LeftActingDomain(L);
    one := One(F);
    spaces:= [ ShallowCopy( BasisVectors( Basis(L) ) ) ];
    for i in [1..Dimension(H)] do
        h:= Basis(H)[i];
        sp0:= [ ];
        for s in spaces do
            sp:= Basis( Subspace( L, s ), s );
            A:= List( s, u -> Coefficients( sp, h*u ) );
            fct:= Set( Factors( CharacteristicPolynomial(A) ) );

            num  := IndeterminateNumberOfUnivariateLaurentPolynomial(fct[1]);
            fam  := FamilyObj( fct[1] );

            fct0:= [ ];
            for f in fct do
               if Degree(f) = 1 then
                  Add( fct0, f );
               elif Degree(f) = 2 then # we just take square roots...
                  cf := CoefficientsOfUnivariatePolynomial(f);
                  b  := cf[2]; 
                  c  := cf[1]; 
                  r  := (-b+Sqrt(b^2-4*c))/2;
                  if not r in F then Error("cannot do this over ",F); fi;
                  Add( fct0, PolynomialByExtRep( fam, 
                                            [ [], -r, [num,1], one] ) );
                  r  := (-b-Sqrt(b^2-4*c))/2;
                  if not r in F then Error("cannot do this over ",F); fi;
                  Add( fct0, PolynomialByExtRep( fam, 
                           [ [], -r, [num,1], one] ) );
               else
                  Error("not split!");
                  return fail;
               fi;
            od;

            for f in fct0 do
                V:= NullspaceMat( Value( f, A ) );
                Add( sp0, List( V, x -> LinearCombination( s, x ) ) );
            od;
        od;
        spaces:= sp0;
    od;
   
    bas:= Concatenation( spaces );
    C:= List( bas, x -> Coefficients( Basis(L), x ) );
    return TransposedMat(C); 

end;

#############################################################################

corelg.toruseqns:= function( L, H ) # H CSA, find equations for the algebraic 
                             # group with Lie alg ad H.

    local C, ad, add, diag, lat, n, A, d, r, S, P, Q, names, inds, i, j, R, 
          eqns, indet, lhs, rhs, a;

    n:= Dimension(L);

    C:= corelg.diagmat(L,H);
    ad:= List( Basis(H), u -> AdjointMatrix( Basis(L), u ) );
    add:= List( ad, u -> C^-1*u*C );

    diag:= List( add, a -> List( [1..n], i -> a[i][i] ) );
    lat:= NullspaceMat( TransposedMat(diag) );

    # saturation...
    A:= SqrtFieldMakeRational(lat);
    d:= Lcm( List( Filtered( Flat(A), x -> not IsZero(x) ), DenominatorRat ) );
    A:= d*A;
    r:=SmithNormalFormIntegerMatTransforms(A);
    Q:= r.coltrans; P:= r.rowtrans;
    S:= List( r.normal, ShallowCopy );
    for i in [1..Length(S)] do 
        if S[i][i] > 0 then S[i][i]:= 1; fi;
    od;
    lat:= P^-1*S*Q^-1;

    names:= [ ];
    inds:= [ ];
    for i in [1..n] do
        for j in [1..n] do  
            Add( names, Concatenation( "a", String(i), "_", String(j) ) );
            Add( inds, [i,j] ); 
        od;
    od;
    R:= PolynomialRing( SqrtField, names );
    
    eqns:= [ ];
    indet:= IndeterminatesOfPolynomialRing( R );
    for i in [1..n] do
        for j in [1..n] do
            if i <> j then
               Add( eqns, indet[ Position(inds,[i,j]) ] );
            fi;
        od;
    od;

    for i in [1..Length(lat)] do
        lhs:= One(R); rhs:= One(R);
        for j in [1..n] do
            if not IsZero(lat[i][j]) then
               if lat[i][j] > 0 then
                  lhs:= lhs*indet[ Position(inds,[j,j]) ]^lat[i][j];
               else
                  rhs:= rhs*indet[ Position(inds,[j,j]) ]^-lat[i][j];
               fi;
            fi;
        od;
        Add( eqns, lhs-rhs );
    od;
    
    A:= List( [1..n], i -> List( [1..n], j -> indet[ Position(inds,[i,j]) ] ) );
Display(A);
    A:= C^-1*A*C;

    a:= [ ];
    for i in [1..n] do for j in [1..n] do Add( a, A[i][j] ); od; od;
    return List( eqns, p -> Value( p, indet, a ) );
    

end;

#############################################################################

corelg.torusparam:= function( L, H ) # H CSA, find a parametrisation for the algebraic 
                              # group with Lie alg ad H.

    local C, ad, add, diag, lat, n, A, d, r, S, P, Q, names, i, j, R, 
          indet, mats, m;

    n:= Dimension(L);

    C:= corelg.diagmat(L,H);
    ad:= List( Basis(H), u -> AdjointMatrix( Basis(L), u ) );
    add:= List( ad, u -> C^-1*u*C );

    diag:= List( add, a -> List( [1..n], i -> a[i][i] ) );
    lat:= NullspaceMat( TransposedMat(diag) );

    # saturation...
    A:= SqrtFieldMakeRational(lat);
    d:= Lcm( List( Filtered( Flat(A), x -> not IsZero(x) ), DenominatorRat ) );
    A:= d*A;
    r:=SmithNormalFormIntegerMatTransforms(A);
    Q:= r.coltrans; P:= r.rowtrans;
    S:= List( r.normal, ShallowCopy );
    for i in [1..Length(S)] do 
        if S[i][i] > 0 then S[i][i]:= 1; fi;
    od;
    lat:= P^-1*S*Q^-1;

    lat:= NullspaceMat( TransposedMat(lat) );
    # another saturation...
    A:= lat;
    d:= Lcm( List( Filtered( Flat(A), x -> not IsZero(x) ), DenominatorRat ) );
    A:= d*A;
    r:=SmithNormalFormIntegerMatTransforms(A);
    Q:= r.coltrans; P:= r.rowtrans;
    S:= List( r.normal, ShallowCopy );
    for i in [1..Length(S)] do 
        if S[i][i] > 0 then S[i][i]:= 1; fi;
    od;
    lat:= P^-1*S*Q^-1;

    names:= [ ];
    for i in [1..Length(lat)] do
        Add( names, Concatenation( "a", String(i) ) );
    od;
    R:= PolynomialRing( SqrtField, names );
    
    indet:= IndeterminatesOfPolynomialRing( R );
    mats:=  [ ];
    for i in [1..Length(lat)] do
        m:= NullMat( n, n, R );
        for j in [1..n] do
            m[j][j]:= indet[i]^lat[i][j];
        od;
        Add( mats, C*m*C^-1 );
    od;

    return mats;
    

end;

#############################################################################

corelg.resmat:= function( L, T, c ) # T list of mats (parametrised torus),
                             # c carrier algebra,
                             # compute a normalised basis of c_1,
                             # and the restriction of the elements of T...,

   local vv, cfv, hdv, i, k, perm, mats, A, m, c0, cc, names, P, g0, sp, cf, b0;

   vv:= CanonicalBasis( Subspace( L, c.gp[1] ) );
   cfv:= List( vv, u -> Coefficients( Basis(L), u ) );
   hdv:= [ ];
   for i in [1..Length(cfv)] do
       k:= 1;
       while IsZero(cfv[i][k]) do k:= k+1; od;
       Add( hdv, k );
   od;

   perm:= Sortex( hdv );
   cfv:= Permuted( cfv, perm );
   vv:= Permuted( BasisVectors( vv ), perm );

   mats:= [ ];
   for A in T do
       m:= [ ];
       for i in [1..Length(cfv)] do
           c0:= A*cfv[i];
           cc:= [ ];
           for k in [1..Length(hdv)] do
               cf:= c0[hdv[k]]/cfv[k][hdv[k]];
               c0:= c0 -cf*cfv[k];
               Add( cc, cf );
           od;
           Add( m, cc );
       od;
       if not IsZero( c0 ) then
          Add( mats, fail );
       else
          Add( mats, m );
       fi;
   od;

   names:= List( [1..Length(vv)], i -> Concatenation("x",String(i)));
   P:= PolynomialRing( SqrtField, names );
   b0:= Concatenation( List( c.gp[1], x -> List( c.gn[1], y -> x*y ) ) );
   Append( b0, Concatenation( List( c.g0, x -> List( c.g0, y -> x*y ) ) ) );
   g0:= Subalgebra( L, b0 );
   sp:= Basis( Subspace( L, vv ), vv );
   m:= List( CanonicalBasis(g0), x -> TransposedMat(
                      List( vv, v -> Coefficients(sp,x*v)) ) );
   m:= List( m, u -> u*IndeterminatesOfPolynomialRing(P));
   
   return rec( bas:= vv, mats:= mats, densep:= DeterminantMat(m), P:= P );

end;

#############################################################################

corelg.resmat_mat:= function( L, T, c ) # T list of mats (parametrised torus),
                             # c carrier algebra,
                             # compute a normalised basis of c_1,
                             # and the restriction of the elements of T...,

   local vv, cfv, hdv, i, k, perm, mats, A, m, c0, cc, names, P, g0, sp, cf, b0;

   vv:= CanonicalBasis( Subspace( L, c.gp[1] ) );
   cfv:= List( vv, u -> Coefficients( Basis(L), u ) );
   hdv:= [ ];
   for i in [1..Length(cfv)] do
       k:= 1;
       while IsZero(cfv[i][k]) do k:= k+1; od;
       Add( hdv, k );
   od;

   perm:= Sortex( hdv );
   cfv:= Permuted( cfv, perm );
   vv:= Permuted( BasisVectors( vv ), perm );

   mats:= [ ];
   for A in T do
       m:= [ ];
       for i in [1..Length(cfv)] do
           c0:= A*cfv[i];
           cc:= [ ];
           for k in [1..Length(hdv)] do
               cf:= c0[hdv[k]]/cfv[k][hdv[k]];
               c0:= c0 -cf*cfv[k];
               Add( cc, cf );
           od;
           Add( m, cc );
       od;
       if not IsZero( c0 ) then
          Add( mats, fail );
       else
          Add( mats, m );
       fi;
   od;

   names:= List( [1..Length(vv)], i -> Concatenation("x",String(i)));
   P:= PolynomialRing( SqrtField, names );
   b0:= Concatenation( List( c.gp[1], x -> List( c.gn[1], y -> x*y ) ) );
   Append( b0, Concatenation( List( c.g0, x -> List( c.g0, y -> x*y ) ) ) );
   g0:= Subalgebra( L, b0 );
   sp:= Basis( Subspace( L, vv ), vv );
   m:= List( CanonicalBasis(g0), x -> TransposedMat(
                      List( vv, v -> Coefficients(sp,x*v)) ) );
   m:= List( m, u -> u*IndeterminatesOfPolynomialRing(P));
   
   return rec( bas:= vv, mats:= mats, densep:= m, P:= P );

end;

#############################################################################


corelg.expmat:= function( L, u, c1 )

   local R, s, B, m, ex, k, M;

   R:= PolynomialRing( LeftActingDomain(L), ["s"] );
   s:= IndeterminatesOfPolynomialRing(R)[1];
   B:= Basis( Subspace( L, c1 ), c1 );
   m:= TransposedMat( List( B, v -> Coefficients(B,u*v) ) );
   
   ex:=  m^0;
   M:= s*m;
   k:= 1;
   while not IsZero(M) do
         ex:= ex + M/(Factorial(k)*One(LeftActingDomain(L)));
         M:= M*(s*m);
         k:= k+1;
   od;

   return ex; 

end; 

#############################################################################

corelg.expmat0:= function( L, u, c1, s )

   local B, m, ex, k, M;

   B:= Basis( Subspace( L, c1 ), c1 );
   m:= TransposedMat( List( B, v -> Coefficients(B,u*v) ) );
   
   ex:=  m^0;
   M:= s*m;
   k:= 1;
   while not IsZero(M) do
         ex:= ex + M/(Factorial(k)*One(LeftActingDomain(L)));
         M:= M*(s*m);
         k:= k+1;
   od;

   return ex; 

end; 


#############################################################################

corelg.Cayleytriples:= function( L, L0, c ) 
      # L the ambient Lie algebra,
      # L0 zero component of L,
      # c carrier algebra.

   local t, by, gam, i, j, k, u, e, K, sl2, h, bh, nu, del, names, P, eqns, eqn, cf;

   t:= Length(c.gp[1]);
   by:= Basis( Subspace( L, c.gn[1] ), c.gn[1] );
   gam:= [];
   for i in [1..t] do
       gam[i]:= Coefficients(by,CartanDecomposition(L).CartanInv(c.gp[1][i]));
   od;
   gam:= SqrtFieldMakeRational(gam);

   # get the defining elt...
   e:= corelg.finde(L,L0,c);
   K:= Subalgebra( L, Concatenation( c.g0, Flat(c.gp), Flat(c.gn) ) );
   sl2:= SL2Triple( K, e );
   Print( Coefficients( by, sl2[1] ),"\n");
   h:= sl2[2];

   bh:= Basis( Subspace( L, c.g0 ), c.g0 );
   nu:= SqrtFieldMakeRational( Coefficients( bh, h ) ); 

   del:= List( [1..t], x -> [ ] );
   for i in [1..t] do
       for j in [1..t] do
           del[i][j]:= SqrtFieldMakeRational( Coefficients( bh, c.gp[1][i]*c.gn[1][j] ) );
       od;
   od;

   names:= List( [1..t], i -> Concatenation("e",String(i)));
   P:= PolynomialRing( Rationals, names );   
   e:= IndeterminatesOfPolynomialRing( P );
   eqns:= [ ];
   for u in [1..t] do
       eqn:= Zero(P);
       for k in [1..t] do
           for i in [1..t] do
               cf:= 0;
               for j in [1..t] do
                   cf:= cf+gam[i][j]*del[k][j][u];
               od;
               eqn:= eqn + cf*e[k]*e[i];
           od;
       od;
       Add( eqns, eqn+nu[u] );
   od;

   return eqns;

end;
