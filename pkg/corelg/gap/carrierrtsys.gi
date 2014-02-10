###########################################################################
# functions for constructing certain carrier algebras and root systems
#
#
# the functions contained in this file are
# 
#  RootSystem
#  RootsystemOfCartanSubalgebra
#  LieAlgebraIsomorphismByCanonicalGenerators
#  ChevalleyBasis
#  RootSystemOfZGradedLieAlgebra
#  RegularCarrierAlgebraOfSL2Triple:
#  
#  corelg.WDD   
#  corelg.CartanMatrixOfCanonicalGeneratingSet
#  corelg.myChevalleyBasis
#  corelg.carrZm
#  corelg.rtsys_withgrad
#  corelg.gradedSubalgebraByCharacteristic
#  corelg.carrierAlgebraBySL2Triple
#  corelg.nil_orbs_outer
#

###########################################################################
#Input:  lie algebra L over Gaussian rationals, a characteristic h, 
#        a signature table T
#Output: the wdd's of h
#
corelg.WDD:= function( L, h, T )
local sl2, K, H, ch, t, cc, tableM, adh, possibles, adh2, fac2, myfactors,
      k, dim, rank, poss0, l, c1, c2, fac, f, inds1, inds2, ev, pos,
      myMatrixOfAction;
   
   myMatrixOfAction := function(L,bas,h)
   local cf, tmp, pos, bas2;
      cf  := Coefficients(Basis(L),h);
      pos := Filtered([1..Length(cf)],x->not cf[x]=0*cf[x]);
      bas2 := Basis(L){pos};
      cf   := cf{pos};
      tmp := List(bas2,x->TransposedMatDestructive( 
                   List(bas,b-> Coefficients( bas, x^b ) )));
      return tmp*cf;
   end;

   if T.tipo = "notD" then
      T    := T.tab;
      adh  := TransposedMat( AdjointMatrix( Basis(L), h ) );
      adh  := List(adh,x->List(x,SqrtFieldEltByCyclotomic));
      possibles:= [ ];
      dim  := Dimension(L);
      rank := RankMat(adh);
      for k in [1..Length(T)] do
         if T[k][1][1] = dim-rank then
            Add( possibles, T[k] );
         fi;
      od;
   
      k:= 1;
      while Length(possibles) > 1 do 
         rank  := RankMat( adh-k*adh^0 );
         poss0 := [ ];
         for l in [1..Length(possibles)] do
            if possibles[l][1][k+1] = dim-rank then
               Add( poss0, possibles[l] );
            fi;
          od;
          possibles:= poss0;
          k:= k+1;
      od; 

      return possibles[1][2];         

   else

      myfactors := function(mm)
      local d, n, x, fx, r, i,m, tmp;
         m  := List(mm,x->List(x,SqrtFieldEltByCyclotomic));
         d  := Length(m);
         n  := 0;
         x  := Indeterminate(SqrtField,"x");
         fx := [];
         while Length(fx)<d do
            tmp := List(m,x->List(x,ShallowCopy));
            for i in [1..Length(tmp)] do tmp[i][i] := tmp[i][i]-n; od;
            r := d-RankMat(tmp);  
            if r>0 then 
               for i in [1..r] do 
                  Add(fx,x-(n*One(SqrtField)));
                  if n > 0 then Add(fx,x+(n*One(SqrtField))); fi;
               od; 
            fi;
            n := n+1;
            if n>1000 then return Error("mhm...n greater 1000"); fi;
         od;
         return List(fx,SqrtFieldPolynomialToRationalPolynomial);
      end;

      adh := myMatrixOfAction(L, Basis( T.V1 ), h );
      c1  := [ ];
     #adh := List(adh,x->List(x,SqrtFieldEltByCyclotomic));
     #fac := Factors( SqrtFieldPolynomialToRationalPolynomial(
     #                       CharacteristicPolynomial( adh ) ) );
      fac := myfactors(adh);
      for f in fac do
         ev:= ExtRepPolynomialRatFun( f );
         if ev[1] = [] then 
            ev:= -ev[2];
         else
            ev:= 0;
         fi;
         pos:= PositionProperty( c1, x -> x[1]=ev );
         if pos = fail then
            Add( c1, [ev,1] );
         else
            c1[pos][2]:= c1[pos][2]+1;
         fi;
     od;
     Sort( c1, function(a,b) return a[1]<b[1]; end );

     adh := myMatrixOfAction(L, Basis( T.V2 ), h );
     c2  := [];
  
     
    
    #adh2 := List(adh,x->List(x,SqrtFieldEltByCyclotomic));
    #fac2 := Factors( SqrtFieldPolynomialToRationalPolynomial(
    #                        CharacteristicPolynomial( adh2 ) ) );

     fac  := myfactors(adh);
    #if not AsSet(fac)=AsSet(fac2) then Error("mhm..."); fi;
     for f in fac do
        ev := ExtRepPolynomialRatFun( f );
        if ev[1] = [] then 
           ev := -ev[2];
        else
           ev := 0;
        fi;
        pos:= PositionProperty( c2, x -> x[1]=ev );
        if pos = fail then
           Add( c2, [ev,1] );
        else
           c2[pos][2]:= c2[pos][2]+1;
        fi;
     od;
     Sort( c2, function(a,b) return a[1]<b[1]; end );

     inds1:= Filtered( [1..Length(T.tab1)], x -> T.tab1[x][1]=c1 );
     inds2:= Filtered( [1..Length(T.tab2)], x -> T.tab2[x][1]=c2 );
     return T.tab1[ Intersection( inds1, inds2 )[1] ][2];

   fi;
end;




##############################################################################
##
#F   RootsystemOfCartanSubalgebra( <L> )
#F   RootsystemOfCartanSubalgebra( <L>, <H> )
##
##   <L> is a semisimple lie algebra over Gaussian rationals or SqrtField;
##   this function returns a rootsystem of <L> with respect to <H>, and
##   with repect to CartanSubalgebra(<L>) if <H> is not provided
##
InstallGlobalFunction( RootsystemOfCartanSubalgebra, function( arg ) 

    local F,          # coefficients domain of `L'
          BL,         # basis of `L'
          H,          # A Cartan subalgebra of `L'
          basH,       # A basis of `H'
          sp,         # A vector space
          B,          # A list of bases of subspaces of `L' whose direct sum
                      # is equal to `L'
          newB,       # A new version of `B' being constructed
          i,j,l,      # Loop variables
          facs,       # List of the factors of `p'
          V,          # A basis of a subspace of `L'
          M,          # A matrix
          cf,         # A scalar
          a,          # A root vector
          ind,        # An index
          basR,       # A basis of the root system
          h,          # An element of `H'
          posR,       # A list of the positive roots
          fundR,      # A list of the fundamental roots
          issum,      # A boolean
          CartInt,    # The function that calculates the Cartan integer of
                      # two roots
          C,          # The Cartan matrix
          S,          # A list of the root vectors
          zero,       # zero of `F'
          hts,        # A list of the heights of the root vectors
          sorh,       # The set `Set( hts )'
          sorR,       # The soreted set of roots
          R,          # The root system.
          Rvecs,      # The root vectors.
          x,y,        # Canonical generators.
          noPosR,     # Number of positive roots.
          facs0, num, fam, f, b, c, r, F0, Mold, one, t1, t2, t3, L; 

    # Let a and b be two roots of the rootsystem R.
    # Let s and t be the largest integers such that a-s*b and a+t*b
    # are roots.
    # Then the Cartan integer of a and b is s-t.

    L := arg[1];
    if Length(arg)=2 then H := arg[2]; else H := CartanSubalgebra(L); fi;

    if HasRootSystem(H) then
       return RootSystem(H);
    fi;
   
    CartInt := function( R, a, b )
       local s,t,rt;
       s:=0; t:=0;
       rt:=a-b;
       while (rt in R) or (rt=0*R[1]) do
         rt:=rt-b;
         s:=s+1;
       od;
       rt:=a+b;
       while (rt in R) or (rt=0*R[1]) do
         rt:=rt+b;
         t:=t+1;
       od;
       return s-t;
    end;

    F   := LeftActingDomain( L );
    one := One(F);

    # removed, to speed things up a bit:
    #if Determinant( KillingMatrix( Basis( L ) ) ) = Zero( F ) then
    #  Error("the Killing form of <L> is degenerate" );
    #  return fail;
    #fi;


    # First we compute the common eigenvectors of the adjoint action of a
    # Cartan subalgebra H. Here B will be a list of bases of subspaces
    # of L such that H maps each element of B into itself.
    # Furthermore, B has maximal length w.r.t. this property.

    BL   := Basis( L );
    B    := [ ShallowCopy( BasisVectors( BL ) ) ];
    basH := BasisVectors( Basis( H ) );

    for i in basH do
     #Print("now ",Position(basH,i)," of basH\n");

      newB := [ ];
      for j in B do
          #Print("  now ",Position(B,j)," of B\n");

        if Length(j) = 1 then
           Add( newB, j ); 
        else
           V    := Basis( VectorSpace( F, j, "basis" ), j );
           Mold := List( j, x -> Coefficients( V, i*x ) );

           if fail in Flat(Mold) then
              Info(InfoCorelg,1,"Extension of base field would be necessary; have to return fail");
              return fail;
           fi;

           if IsSqrtField(F) then
              M    := SqrtFieldMakeRational(Mold);
              if M = false then 
                #Error("matrix we want to compute char pol of cannot be made rationals");
                #Print(" matrix cannot be made rations; use CharPol for SqrtField\n");
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
                  r  := (-b+Sqrt(b^2-4*c))/2;  #have Sqrt method for rat in SqrtField!
                  if not r in F then Error("cannot do this over ",F); fi;
                  Add( facs0, PolynomialByExtRep( fam, [ [], -r, [num,1], one] ) );
                  r  := (-b-Sqrt(b^2-4*c))/2;
                  if not r in F then Error("cannot do this over ",F); fi;
                  Add( facs0, PolynomialByExtRep( fam, [ [], -r, [num,1], one] ) );

               else
                  Error("not split");
                  return fail;
               fi;
           od;

           for l in facs0 do
             #t1 := Runtime();
              V := NullspaceMat( Value( l, Mold ) );
             #t2 := Runtime();
              Add( newB, List( V, x -> LinearCombination( j, x ) ) );
             #t3 := Runtime();
             #Print("ns, lc",t2-t1," ", t3-t2,"\n");
           od;
        fi;

      od;
      B:= newB;

   od;

  # Now we throw away the subspace H.
    B:= Filtered( B, x -> ( not x[1] in H ) );

  # If an element of B is not one dimensional then H does not split
  # completely, and hence we cannot compute the root system.

   for i in [ 1 .. Length(B) ] do
      if Length( B[i] ) <> 1 then
         Error("the Cartan subalgebra of <L> in not split" );
         return fail;
      fi;
   od;

  # Now we compute the set of roots S.
  # A root is just the list of eigenvalues of the basis elements of H
  # on an element of B.

   S    := [];


   zero := Zero( F );
   for i in [ 1 .. Length(B) ] do
      a   := [ ];
      ind := 0;
      cf  := zero;
      while cf = zero do
         ind := ind+1;
         cf  := Coefficients( BL, B[i][1] )[ ind ];
      od;
      for j in [1..Length(basH)] do
         Add( a, Coefficients( BL, basH[j]*B[i][1] )[ind] / cf );
      od;
      Add( S, a );
   od;

   Rvecs := List( B, x -> x[1] );

  # A set of roots basR is calculated such that the set
  # { [ x_r, x_{-r} ] | r\in R } is a basis of H.

   basH := [ ];
   basR := [ ];
   sp   := MutableBasis( F, [], Zero(L) );
   i    :=1;
   while Length( basH ) < Dimension( H ) do
      a:= S[i];
      j:= Position( S, -a );
      h:= B[i][1]*B[j][1];
      if not IsContainedInSpan( sp, h ) then
      #if not corelg.eltInSubspace(L,BasisVectors(sp), h) then
         CloseMutableBasis( sp, h );
         Add( basR, a );
         Add( basH, h );
      fi;
      i:=i+1;
   od;

  # A root a is said to be positive if the first nonzero element of
  # [ CartInt( S, a, basR[j] ) ] is positive.
  # We calculate the set of positive roots.

   posR:= [ ];
   i:=1;
   while Length( posR ) < Length( S )/2 do
      a:= S[i];
      if (not a in posR) and (not -a in posR) then
         cf := 0;
         j  := 0;
         while cf = 0 do
            j  := j+1;
            cf := CartInt( S, a, basR[j] );
         od;
         if 0 < cf then
            Add( posR, a );
         else
            Add( posR, -a );
         fi;
      fi;
      i:=i+1;
   od;

  # A positive root is called simple if it is not the sum of two other
  # positive roots.
  # We calculate the set of simple roots fundR.

    fundR:= [ ];
   for a in posR do
      issum:= false;
      for i in [1..Length(posR)] do
         for j in [i+1..Length(posR)] do
            if a = posR[i]+posR[j] then
               issum:=true;
            fi;
         od;
      od;
      if not issum then
         Add( fundR, a );
      fi;
   od;

  # Now we calculate the Cartan matrix C of the root system.

   C:= List( fundR, i -> List( fundR, j -> CartInt( S, i, j ) ) );

  # Every root can be written as a sum of the simple roots.
  # The height of a root is the sum of the coefficients appearing
  # in that expression.
  # We order the roots according to increasing height.

   V    := BasisNC( VectorSpace( F, fundR ), fundR );
   hts  := List( posR, r -> Sum( Coefficients( V, r ) ) );
   sorh := Set( hts );

   sorR:= [ ];
   for i in [1..Length(sorh)] do
      Append( sorR, Filtered( posR, r -> hts[Position(posR,r)] = sorh[i] ) );
   od;
   Append( sorR, -1*sorR );
   Rvecs:= List( sorR, r -> Rvecs[ Position(S,r) ] );
    
  # We calculate a set of canonical generators of L. Those are elements
  # x_i, y_i, h_i such that h_i=x_i*y_i, h_i*x_j = c_{ij} x_j,
  # h_i*y_j = -c_{ij} y_j for i \in {1..rank}
    
   x:= Rvecs{[1..Length(C)]};
   noPosR:= Length( Rvecs )/2;
   y:= Rvecs{[1+noPosR..Length(C)+noPosR]};
   for i in [1..Length(x)] do
      V:= VectorSpace( LeftActingDomain(L), [ x[i] ] );
      B:= Basis( V, [x[i]] );
      y[i]:= y[i]*2/Coefficients( B, (x[i]*y[i])*x[i] )[1];
   od;
    
   h:= List([1..Length(C)], j -> x[j]*y[j] );
    
  # Now we construct the root system, and install as many attributes
  # as possible. The roots are represented als lists [ \alpha(h_1),....
  # ,\alpha(h_l)], where the h_i form the Cartan part of the canonical
  # generators.
    
   R:= Objectify( NewType( NewFamily( "RootSystemFam", IsObject ),
               IsAttributeStoringRep and IsRootSystemFromLieAlgebra ), 
               rec() );
   SetCanonicalGenerators( R, [ x, y, h ] );
   SetUnderlyingLieAlgebra( R, L );
   SetPositiveRootVectors( R, Rvecs{[1..noPosR]});
   SetNegativeRootVectors( R, Rvecs{[noPosR+1..2*noPosR]} );
   SetCartanMatrix( R, C );
    
   posR:= [ ];
   for i in [1..noPosR] do
      B:= Basis( VectorSpace( F, [ Rvecs[i] ] ), [ Rvecs[i] ] );
      posR[i]:= List( h, hj ->  Coefficients( B, hj*Rvecs[i] )[1] );
   od;

  #roots are rationals
   if IsSqrtField(F) then
      posR := List(posR, x-> List(x, SqrtFieldEltToCyclotomic));
   fi;
 
   SetPositiveRoots( R, posR );
   SetNegativeRoots( R, -posR ); 
   SetSimpleSystem( R, posR{[1..Length(C)]} );
   SetRootSystem(H,R);
   return R;
end);



#####################################################################
corelg.CartanMatrixOfCanonicalGeneratingSet := function(L,R)

   local mat, i, j, u, k;

   mat:= List( R[3], x -> [] );
   for i in [1..Length(R[3])] do
       for j in [1..Length(R[3])] do
           if i = j then 
              mat[i][j]:= 2;
           else
              u:= R[3][j]*R[1][i];
              k:= -3;
              while not u = k*R[1][i] do k:= k+1; od;
              mat[i][j]:= k;
           fi;
        od;
   od;
   return mat;
  #return List(R[1],e-> List(R[3],h->  
  #           Coefficients(Basis(Subspace(L,[e]),[e]),h*e)[1]));
end;



######################################################################
InstallMethod( RootSystem,
   "for Lie algebras",
   true,
   [ IsLieAlgebra ], 0, function(L)

return RootsystemOfCartanSubalgebra( L, CartanSubalgebra(L) );

end );


##############################################################################
##
#F LieAlgebraIsomorphismByCanonicalGenerators( <L1>, <R1>, <L2>, <R2> )
##
## <L1> and <L2> both are semisimple lie algebras over Gaussian rationals or SqrtField
## and either <R1> and <R2> both are canonical generators of <L1> and <L2> defining 
## the same Cartan Matrix, or <R1> and <R2> are rootsystems or Cartan subalgebras of 
## <L1> and <L2>, respectively; this functions constructs an isomorphism from
## <L1> to <L2> by mapping canonical generators onto canonical generators.
## Attention: This function does not check whether the map actually is a
##            Lie isomorphisms! 
##
InstallGlobalFunction(LieAlgebraIsomorphismByCanonicalGenerators, function( L1, R1, L2, R2 )
local b1, b2, c1, c2, t, tp, en, i, cm1, cm2, tmp;

 
   #R1 and R2 are canonical generators wrt same ordering
    if IsList(R1) and IsList(R2) then       
      #check if can gens really define the same Cartan Matrix
       cm1 := corelg.CartanMatrixOfCanonicalGeneratingSet(L1,R1);
       cm2 := corelg.CartanMatrixOfCanonicalGeneratingSet(L2,R2);
       if not cm1=cm2 then 
          Error("Cartan Matrices of canonical gen sets are different"); 
       fi;
       b1:= SLAfcts.canbas( L1, R1 );
       b2:= SLAfcts.canbas( L2, R2 );
       tmp := AlgebraHomomorphismByImagesNC( L1, L2, corelg.myflat(b1), corelg.myflat(b2) );
       SetIsIsomorphismOfLieAlgebras(tmp,true);
       return tmp;
    fi;

   #otherwise, R1 and R2 are CSA or Rootsystems
   #construct canonical generators and then isomorphism
    if not IsRootSystem(R1) then R1 := RootsystemOfCartanSubalgebra(L1,R1); fi;
    if not IsRootSystem(R2) then R2 := RootsystemOfCartanSubalgebra(L2,R2); fi;

    t  := CartanType( CartanMatrix(R1) );
    tp := ShallowCopy( t.types );
    en := ShallowCopy( t.enumeration );
    SortParallel( tp, en );
    c1 := [ [], [], [] ];
    for i in [1..Length(en)] do
        Append( c1[1], CanonicalGenerators(R1)[1]{en[i]} );
        Append( c1[2], CanonicalGenerators(R1)[2]{en[i]} );
        Append( c1[3], CanonicalGenerators(R1)[3]{en[i]} );
    od;

    t  := CartanType( CartanMatrix(R2) );
    tp := ShallowCopy( t.types );
    en := ShallowCopy( t.enumeration );
    SortParallel( tp, en );
    c2 := [ [], [], [] ];
    for i in [1..Length(en)] do
        Append( c2[1], CanonicalGenerators(R2)[1]{en[i]} );
        Append( c2[2], CanonicalGenerators(R2)[2]{en[i]} );
        Append( c2[3], CanonicalGenerators(R2)[3]{en[i]} );
    od;

    return LieAlgebraIsomorphismByCanonicalGenerators(L1,c1,L2,c2);
end);


###############################################
#this is print:
InstallMethod( ViewObj,
   "for IsIsomorphismOfLieAlgebras",
   true,
   [ IsIsomorphismOfLieAlgebras ], 100,
   function( o )
   local r,t,m,i,minus,signs, tmp;
   Print(Concatenation(["<Lie algebra isomorphism between Lie algebras of dimension ",
         String(Dimension(Source(o)))," over ",String(LeftActingDomain(Source(o))),">"]));
end );


##############################################################################
corelg.myChevalleyBasis := function(LL,R)

     local tp, K, f, cc, h, pr, xx, yy, i, sp, rt, pos;

     if HasChevalleyBasis(R) then return ChevalleyBasis(R); fi;
     tp := CartanType( CartanMatrix(R) );
     K  := DirectSumOfAlgebras( List( tp.types, x -> 
                SimpleLieAlgebra( x[1], x[2], LeftActingDomain(LL) ) ) );
     f  := LieAlgebraIsomorphismByCanonicalGenerators( K, RootSystem(K), LL, R );
     cc := List( ChevalleyBasis(K), x -> List( x, y -> Image( f, y ) ) );
     h  := CanonicalGenerators(R)[3];
     pr := PositiveRoots(R);
     xx := [ ]; yy:= [ ];
     for i in [1..Length(cc[1])] do
         sp  := BasisNC( SubspaceNC( LL, [cc[1][i]],"basis" ), [ cc[1][i] ] );
         rt  := List( h, u -> Coefficients( sp, u*cc[1][i] )[1] );
         pos := Position( pr, rt );
        #if pos <> i then Print(i,"  ",pos,"\n"); fi;
         xx[pos]:= cc[1][i]; yy[pos]:= cc[2][i];
     od; 

     cc := [ xx, yy, h ];
     SetChevalleyBasis( R, cc );
     return cc;

end;



##############################################################################
##
#F   ChevalleyBasis( <R> );
##
##   <L> is a semisimple lie algebra over Gaussian rationals or SqrtField
##   and <R> is  a rootsystem of <L> (with respect to some Cartan subalgebra);
##   this function returns a Chevalley basis
##   of this rootsystem / the rootsystem of the Cartan subalgebra
##
#InstallGlobalFunction( ChevalleyBasisOfRootsystem, function( L, R)

InstallMethod( ChevalleyBasis,
   "for a root system",
   true,
   [ IsRootSystem ], 0,

function( R )
   local L, R1, cb1, iso, perm, pr, pr1, cb, tmp, pos, i, cg1;
   L    := UnderlyingLieAlgebra(R);
   return corelg.myChevalleyBasis( L, R );
end);




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
InstallGlobalFunction( RootSystemOfZGradedLieAlgebra, function( arg )

    # g a grading in carrier form, ie a record with components g0, gp, gn.

    local F,          # coefficients domain of L
          BL,         # basis of L
          H,          # A Cartan subalgebra of L
          basH,       # A basis of H
          sp,         # A vector space
          B,          # A list of bases of subspaces of L whose direct sum
                      # is equal to L
          newB,       # A new version of B being constructed
          i,j,l,      # Loop variables
          facs,       # List of the factors of p
          V,          # A basis of a subspace of L
          M,          # A matrix
          cf,         # A scalar
          a,          # A root vector
          ind,        # An index
          basR,       # A basis of the root system
          h,          # An element of H
          posR,       # A list of the positive roots
          fundR,      # A list of the fundamental roots
          issum,      # A boolean
          CartInt,    # The function that calculates the Cartan integer of
                      # two roots
          C,          # The Cartan matrix
          S,          # A list of the root vectors
          zero,       # zero of F
          hts,        # A list of the heights of the root vectors
          sorh,       # The set Set( hts )
          sorR,       # The soreted set of roots
          R,          # The root system.
          Rvecs,      # The root vectors.
          x,y,        # Canonical generators.
          noPosR,     # Number of positive roots.
          facs0, num, fam, f, b, c, r, possp, g0, F0, Mold, one, L, g; 

    # Let a and b be two roots of the rootsystem R.
    # Let s and t be the largest integers such that a-s*b and a+t*b
    # are roots.
    # Then the Cartan integer of a and b is s-t.

    CartInt := function( R, a, b )
       local s,t,rt;
       s:=0; t:=0;
       rt:=a-b;
       while (rt in R) or (rt=0*R[1]) do
         rt:=rt-b;
         s:=s+1;
       od;
       rt:=a+b;
       while (rt in R) or (rt=0*R[1]) do
         rt:=rt+b;
         t:=t+1;
       od;
       return s-t;
    end;

    L:= arg[1]; g:= arg[2];
    if Length(arg) = 3 then
       H:= arg[3];
    else
       H:= CartanSubalgebra(L);
    fi;
       

    F   := LeftActingDomain( L );
    one := One(F);
    if DeterminantMat( KillingMatrix( Basis( L ) ) ) = Zero( F ) then
      Error("the Killing form of <L> is degenerate" );
      return fail;
    fi;


    # First we compute the common eigenvectors of the adjoint action of a
    # Cartan subalgebra H. Here B will be a list of bases of subspaces
    # of L such that H maps each element of B into itself.
    # Furthermore, B has maximal length w.r.t. this property.

    BL:= Basis( L );
    B:= [ ShallowCopy( BasisVectors( BL ) ) ];
    basH:= BasisVectors( Basis( H ) );

   
    for i in basH do

      newB:= [ ];
      for j in B do

         if Length(j) = 1 then
           Add( newB, j ); 
         else
           V    := Basis( VectorSpace( F, j, "basis" ), j );
           Mold := List( j, x -> Coefficients( V, i*x ) );
           if IsSqrtField(F) then
              M    := List(Mold,x->List(x,SqrtFieldEltToCyclotomic));
              f    := CharacteristicPolynomial( M );
              facs := Set(Factors( f ));
              f    := SqrtFieldRationalPolynomialToSqrtFieldPolynomial(f);
              facs := Set(List(facs,SqrtFieldRationalPolynomialToSqrtFieldPolynomial));
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
                   if not r in F then Error("cannot do this over ",F); fi;
                  Add( facs0, PolynomialByExtRep( fam, [ [], -r, [num,1], one] ) );
                  r  := (-b-Sqrt(b^2-4*c))/2;
                  if not r in F then Error("cannot do this over ",F); fi;
                  Add( facs0, PolynomialByExtRep( fam, [ [], -r, [num,1], one] ) );

               else
                  Error("not split!");
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

    # Now we throw away the subspace H.

    B := Filtered( B, x -> ( not x[1] in H ) );
    #B:= Filtered( B, x -> ( not corelg.eltInSubspace(L,BasisVectors(Basis(H)),x[1])));
    

    # If an element of B is not one dimensional then H does not split
    # completely, and hence we cannot compute the root system.

    for i in [ 1 .. Length(B) ] do
      if Length( B[i] ) <> 1 then
        Error("the Cartan subalgebra of <L> in not split" );
        return fail;
      fi;
    od;

    # Now we compute the set of roots S.
    # A root is just the list of eigenvalues of the basis elements of H
    # on an element of B.

    S:= [];
    zero:= Zero( F );
    for i in [ 1 .. Length(B) ] do
      a:= [ ];
      ind:= 0;
      cf:= zero;
      while cf = zero do
        ind:= ind+1;
        cf:= Coefficients( BL, B[i][1] )[ ind ];
      od;
      for j in [1..Length(basH)] do
        Add( a, Coefficients( BL, basH[j]*B[i][1] )[ind] / cf );
      od;
      Add( S, a );
    od;

    Rvecs:= List( B, x -> x[1] );

    # A set of roots basR is calculated such that the set
    # { [ x_r, x_{-r} ] | r\in R } is a basis of H.

    basH:= [ ];
    basR:= [ ];
    sp:= MutableBasis( F, [], Zero(L) );
    i:=1;
    while Length( basH ) < Dimension( H ) do
      a:= S[i];
      j:= Position( S, -a );
      h:= B[i][1]*B[j][1];
      #if not corelg.eltInSubspace(L,BasisVectors(sp),h ) then
      if not IsContainedInSpan( sp, h ) then
        CloseMutableBasis( sp, h );
        Add( basR, a );
        Add( basH, h );
      fi;
      i:=i+1;
    od;



    # A root a is said to be positive if the corr root space lies in g_k with k>0,
    # or in g_k with k=0 and the first nonzero element of
    # [ CartInt( S, a, basR[j] ) ] is positive.
    # We calculate the set of positive roots.

    posR:= [ ];
    i:=1;
    possp:= SubspaceNC( L, Concatenation( g.gp ),"basis" );
    g0:= SubspaceNC( L, g.g0,"basis" );
    while Length( posR ) < Length( S )/2 do
      a:= S[i];
      if (not a in posR) and (not -a in posR) then

        if B[i][1] in possp then
        #if corelg.eltInSubspace(L,Basis(possp),B[i][1] ) then
           Add( posR, a );
        elif B[i][1] in g0 then
           cf:= 0;
           j:= 0;
           while cf = 0 do
             j:= j+1;
             cf:= CartInt( S, a, basR[j] );
           od;
           if 0 < cf then
             Add( posR, a );
           else
             Add( posR, -a );
           fi;
        else
           Add( posR, -a );
        fi;
      fi;
      i:=i+1;
    od;

    # A positive root is called simple if it is not the sum of two other
    # positive roots.
    # We calculate the set of simple roots fundR.

    fundR:= [ ];
    for a in posR do
      issum:= false;
      for i in [1..Length(posR)] do
        for j in [i+1..Length(posR)] do
          if a = posR[i]+posR[j] then
            issum:=true;
          fi;
        od;
      od;
      if not issum then
        Add( fundR, a );
      fi;
    od;

    # Now we calculate the Cartan matrix C of the root system.

    C:= List( fundR, i -> List( fundR, j -> CartInt( S, i, j ) ) );

    # Every root can be written as a sum of the simple roots.
    # The height of a root is the sum of the coefficients appearing
    # in that expression.
    # We order the roots according to increasing height.

    V:= BasisNC( VectorSpace( F, fundR ), fundR );
    hts:= List( posR, r -> Sum( Coefficients( V, r ) ) );
    sorh:= Set( hts );

    sorR:= [ ];
    for i in [1..Length(sorh)] do
      Append( sorR, Filtered( posR, r -> hts[Position(posR,r)] = sorh[i] ) );
    od;
    Append( sorR, -1*sorR );
    Rvecs:= List( sorR, r -> Rvecs[ Position(S,r) ] );
    
    # We calculate a set of canonical generators of L. Those are elements
    # x_i, y_i, h_i such that h_i=x_i*y_i, h_i*x_j = c_{ij} x_j,
    # h_i*y_j = -c_{ij} y_j for i \in {1..rank}
    
    x:= Rvecs{[1..Length(C)]};
    noPosR:= Length( Rvecs )/2;
    y:= Rvecs{[1+noPosR..Length(C)+noPosR]};
    for i in [1..Length(x)] do
        V:= VectorSpace( LeftActingDomain(L), [ x[i] ] );
        B:= Basis( V, [x[i]] );
        y[i]:= y[i]*2/Coefficients( B, (x[i]*y[i])*x[i] )[1];
    od;
    
    h:= List([1..Length(C)], j -> x[j]*y[j] );
    
    # Now we construct the root system, and install as many attributes
    # as possible. The roots are represented als lists [ \alpha(h_1),....
    # ,\alpha(h_l)], where the h_i form the Cartan part of the canonical
    # generators.
    
    R:= Objectify( NewType( NewFamily( "RootSystemFam", IsObject ),
                IsAttributeStoringRep and IsRootSystemFromLieAlgebra ), 
                rec() );
    SetCanonicalGenerators( R, [ x, y, h ] );
    SetUnderlyingLieAlgebra( R, L );
    SetPositiveRootVectors( R, Rvecs{[1..noPosR]});
    SetNegativeRootVectors( R, Rvecs{[noPosR+1..2*noPosR]} );
    SetCartanMatrix( R, C );
    
    posR:= [ ];
    for i in [1..noPosR] do
        B:= Basis( VectorSpace( F, [ Rvecs[i] ] ), [ Rvecs[i] ] );
        posR[i]:= List( h, hj ->  Coefficients( B, hj*Rvecs[i] )[1] );
    od;
    
   
   #roots are rationals
    if IsSqrtField(F) then
       posR := List(posR, x-> List(x, SqrtFieldEltToCyclotomic));
    fi;

    SetPositiveRoots( R, posR );
    SetNegativeRoots( R, -posR ); 
    SetSimpleSystem( R, posR{[1..Length(C)]} );

    return R;
    
end);






########################################################################################
########################################################################################
#
# the following functions define RegularCarrierAlgebraOfSL2Triple:
#    - corelg.carrZm
#    - corelg.rtsys_withgrad
#    - corelg.gradedSubalgebraByCharacteristic
#    - corelg.carrierAlgebraBySL2Triple
########################################################################################
########################################################################################


########################################################################################
corelg.carrZm:= function( L, gr, e )

   local h, lams, sp, i, gp, gn, eigensp, g0, g1, gm, m, K, k, dim,t0;

   sp:= SubalgebraNC( L, gr[1] );
   h:= BasisVectors(CanonicalBasis( CartanSubalgebra( Intersection( sp, LieNormalizer(L,
                            SubalgebraNC(L,[e]))))));
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


###############################################################################

corelg.rtsys_withgrad :=    function( L, rvecs, H, g )

    # g a grading in carrier form, ie a record with components g0, gp, gn.
    # rvecs is a list of the root vectors

    local F,          # coefficients domain of L
          BL,         # basis of L
          basH,       # A basis of H
          sp,         # A vector space
          B,          # A list of bases of subspaces of L whose direct sum
                      # is equal to L
          newB,       # A new version of B being constructed
          i,j,l,      # Loop variables
          facs,       # List of the factors of p
          V,          # A basis of a subspace of L
          M,          # A matrix
          cf,         # A scalar
          a,          # A root vector
          ind,        # An index
          basR,       # A basis of the root system
          h,          # An element of H
          posR,       # A list of the positive roots
          fundR,      # A list of the fundamental roots
          issum,      # A boolean
          CartInt,    # The function that calculates the Cartan integer of
                      # two roots
          C,          # The Cartan matrix
          S,          # A list of the root vectors
          zero,       # zero of F
          hts,        # A list of the heights of the root vectors
          sorh,       # The set Set( hts )
          sorR,       # The soreted set of roots
          R,          # The root system.
          Rvecs,      # The root vectors.
          x,y,        # Canonical generators.
          noPosR,     # Number of positive roots.
          facs0, num, fam, f, b, c, r, possp, g0, v, fct; 

    # Let a and b be two roots of the rootsystem R.
    # Let s and t be the largest integers such that a-s*b and a+t*b
    # are roots.
    # Then the Cartan integer of a and b is s-t.

    CartInt := function( R, a, b )
       local s,t,rt;
       s:=0; t:=0;
       rt:=a-b;
       while (rt in R) or (rt=0*R[1]) do
         rt:=rt-b;
         s:=s+1;
       od;
       rt:=a+b;
       while (rt in R) or (rt=0*R[1]) do
         rt:=rt+b;
         t:=t+1;
       od;
       return s-t;
    end;

    F:= LeftActingDomain(L);
    BL:= Basis(L);
    basH:= Basis(H);
    B:= List( rvecs, x -> [x] );

    # Now we compute the set of roots S.
    # A root is just the list of eigenvalues of the basis elements of H
    # on an element of B.

    S:= [];
    zero:= Zero( F );
    for i in [ 1 .. Length(B) ] do
      a:= [ ];
      ind:= 0;
      cf:= zero;
      while cf = zero do
        ind:= ind+1;
        cf:= Coefficients( BL, B[i][1] )[ ind ];
      od;
      for j in [1..Length(basH)] do
        Add( a, Coefficients( BL, basH[j]*B[i][1] )[ind] / cf );
      od;
      Add( S, a );
    od;

    Rvecs:= List( B, x -> x[1] );

    # A set of roots basR is calculated such that the set
    # { [ x_r, x_{-r} ] | r\in R } is a basis of H.

    basH:= [ ];
    basR:= [ ];
    sp:= MutableBasis( F, [], Zero(L) );
    i:=1;
    while Length( basH ) < Dimension( H ) do
      a:= S[i];
      j:= Position( S, -a );
      h:= B[i][1]*B[j][1];
      #if not corelg.eltInSubspace(L,BasisVectors(sp), h) then
      if not IsContainedInSpan( sp, h ) then
        CloseMutableBasis( sp, h );
        Add( basR, a );
        Add( basH, h );
      fi;
      i:=i+1;
    od;

    # A root a is said to be positive if the corr root space lies in g_k with k>0,
    # or in g_k with k=0 and the first nonzero element of
    # [ CartInt( S, a, basR[j] ) ] is positive.
    # We calculate the set of positive roots.

    posR:= [ ];
    i:=1;
    possp:= SubspaceNC( L, Concatenation( g.gp ),"basis" );
    g0:= SubspaceNC( L, g.g0,"basis" );
    while Length( posR ) < Length( S )/2 do
      a:= S[i];
      if (not a in posR) and (not -a in posR) then

        if B[i][1] in possp then
           Add( posR, a );
        elif B[i][1] in g0 then
           cf:= zero;
           j:= 0;
           while cf = zero do
             j:= j+1;
             cf:= CartInt( S, a, basR[j] );
           od;
           if 0 < cf then
             Add( posR, a );
           else
             Add( posR, -a );
           fi;
        else
           Add( posR, -a );
        fi;
      fi;
      i:=i+1;
    od;

    # A positive root is called simple if it is not the sum of two other
    # positive roots.
    # We calculate the set of simple roots fundR.

    fundR:= [ ];
    for a in posR do
      issum:= false;
      for i in [1..Length(posR)] do
        for j in [i+1..Length(posR)] do
          if a = posR[i]+posR[j] then
            issum:=true;
          fi;
        od;
      od;
      if not issum then
        Add( fundR, a );
      fi;
    od;

    # Now we calculate the Cartan matrix C of the root system.

    C:= List( fundR, i -> List( fundR, j -> CartInt( S, i, j ) ) );

    # Every root can be written as a sum of the simple roots.
    # The height of a root is the sum of the coefficients appearing
    # in that expression.
    # We order the roots according to increasing height.

    V:= BasisNC( VectorSpace( F, fundR ), fundR );
    hts:= List( posR, r -> Sum( Coefficients( V, r ) ) );
    sorh:= Set( hts );

    sorR:= [ ];
    for i in [1..Length(sorh)] do
      Append( sorR, Filtered( posR, r -> hts[Position(posR,r)] = sorh[i] ) );
    od;
    Append( sorR, -1*sorR );
    Rvecs:= List( sorR, r -> Rvecs[ Position(S,r) ] );
    
    # We calculate a set of canonical generators of L. Those are elements
    # x_i, y_i, h_i such that h_i=x_i*y_i, h_i*x_j = c_{ij} x_j,
    # h_i*y_j = -c_{ij} y_j for i \in {1..rank}
    
    x:= Rvecs{[1..Length(C)]};
    noPosR:= Length( Rvecs )/2;
    y:= Rvecs{[1+noPosR..Length(C)+noPosR]};
    for i in [1..Length(x)] do
        V:= VectorSpace( LeftActingDomain(L), [ x[i] ] );
        B:= Basis( V, [x[i]] );
        y[i]:= y[i]*2/Coefficients( B, (x[i]*y[i])*x[i] )[1];
    od;
    
    h:= List([1..Length(C)], j -> x[j]*y[j] );
    
    # Now we construct the root system, and install as many attributes
    # as possible. The roots are represented als lists [ \alpha(h_1),....
    # ,\alpha(h_l)], where the h_i form the Cartan part of the canonical
    # generators.
    
    R:= Objectify( NewType( NewFamily( "RootSystemFam", IsObject ),
                IsAttributeStoringRep and IsRootSystemFromLieAlgebra ), 
                rec() );
    SetCanonicalGenerators( R, [ x, y, h ] );
    SetUnderlyingLieAlgebra( R, L );
    SetPositiveRootVectors( R, Rvecs{[1..noPosR]});
    SetNegativeRootVectors( R, Rvecs{[noPosR+1..2*noPosR]} );
    SetCartanMatrix( R, C );

    fct:= function(x) 
       if IsGaussRat(x) then return x;  else return x![1][1][1]; fi;
    end;
    
    posR:= [ ];
    for i in [1..noPosR] do
        B:= Basis( VectorSpace( F, [ Rvecs[i] ] ), [ Rvecs[i] ] );
        v:= List( h, hj ->  Coefficients( B, hj*Rvecs[i] )[1] );
        posR[i]:= List( v, x -> fct(x) );
    od;
    
    SetPositiveRoots( R, posR );
    SetNegativeRoots( R, -posR ); 
    SetSimpleSystem( R, posR{[1..Length(C)]} );

    return R;
    
end;

########################################################################

corelg.gradedSubalgebraByCharacteristic:= function( L, gr, h )
   
    # here L is a Z/2-graded Lie algebra, grading in gr, two element list...
    # h nuetral elt of sl2 triple. We get the Z-graded subalgebra such that
    # g_k = { x\in L \mid x in gr[k mod 2], [h,x] = 2*k*x}

    local adh, id, g0, g1, grad, gp, gn, k, done, cf, sp;

    adh:= TransposedMat( AdjointMatrix( Basis(L), h ) );
    id:= adh^0;
    g0:= SubspaceNC( L, gr[1],"basis" );
    g1:= SubspaceNC( L, gr[2],"basis" );
    grad:= [g0,g1];
    gp:= [ ];
    k:= 1;
    done:= false;
    while not done do
       cf:= NullspaceMat( adh-2*k*id );
       if cf <> [] then
          sp:= Intersection( grad[(k mod 2)+1], SubspaceNC( L, List( cf, c -> c*Basis(L) ) ) );
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
          sp:= Intersection( grad[(k mod 2)+1], SubspaceNC( L, List( cf, c -> c*Basis(L) ) ) );
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


corelg.carrierAlgebraBySL2Triple:= function( L, grad, sl2 )

   local R, B, ch, posR, N, rts, rr, pi, r1, zero, stack, res, r, 
         start, rrr, ips, i, vv, u, h, C, CT, pi_0, pi_1, t, s, pos,
         ct, eqns, rhs, eqn, j, sol, h0, psi0, psi1, good, x, y, es, fs, 
         valmat, val, chars, u0, v, done, gr1, gr2, g2, h_mats1, h_mats2, 
         mat, sl2s, id1, id2, Omega, V, e, ff, found, co, k, sp, extended,
         zz, bas, sim, Bw, W0, types, weights, wrts, tp, a, c, comb, hZ, hs,
         info, posRv, negRv, g0, g1, gm, CM, rr0, l0, l1, gr, deg, R0, gs, grading,
         cardat, U, gsp, grr, r0, gp, gn, K0, rvs, F, fct, rsp;

   gs:= corelg.gradedSubalgebraByCharacteristic( L, grad, sl2[2] );

   F:= LeftActingDomain(L);

   K0:= SubalgebraNC( L, Concatenation( gs.g0, corelg.myflat( gs.gp ), corelg.myflat( gs.gn ) ) );
   K0:= LieDerivedSubalgebra( K0 );
   gs.g0:= BasisVectors( Basis( Intersection( K0, SubspaceNC( L, gs.g0,"basis" ) ) ) );

   rvs:= [ ];
   for v in  PositiveRootVectors(RootSystem(L)) do
       if v in K0 then Add( rvs, v ); fi;
       #if corelg.eltInSubspace(L,BasisVectors(Basis(K0)),v) then Add( rvs, v); fi;
   od;
   for v in  NegativeRootVectors(RootSystem(L)) do
       if v in K0 then Add( rvs, v ); fi;
       #if corelg.eltInSubspace(L,BasisVectors(Basis(K0)),v) then Add( rvs, v); fi;
   od;

   SetCartanSubalgebra( K0, Intersection( CartanSubalgebra(L), K0 ) ); #!!!
## SetCartanSubalgebra( K0, Intersection( MaximallyCompactCartanSubalgebra(L), K0 ) );
   if Length(rvs)+Dimension( CartanSubalgebra(K0) ) <> Dimension(K0) then
      R0:= RootsystemOfCartanSubalgebra(K0);
      rvs:= Concatenation( PositiveRootVectors(R0), NegativeRootVectors(R0) );
   fi;
  
   R0:= corelg.rtsys_withgrad( K0, rvs, Intersection( CartanSubalgebra(L), K0 ), gs ); #!!!
 ##R0:= corelg.rtsys_withgrad( K0, rvs, Intersection( MaximallyCompactCartanSubalgebra(L), K0 ), gs );
   
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
                          #not corelg.eltInSubspace(rsp,BasisVectors(Basis(r.sp)),r1[i]) then
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
             U:= SubalgebraNC( L, Concatenation( r0.g0, corelg.myflat(r0.gp), corelg.myflat(r0.gn) ), "basis" );
             U:= LieDerivedSubalgebra(U);
             r0.g0:= BasisVectors( Basis( Intersection( U, SubspaceNC( L, r0.g0, "basis" ) ) ) );

             return r0;
          fi;
       fi;
   od;

   return "not found!!";

end;


##############################################################################
##
#F RegularCarrierAlgebraOfSL2Triple( <L>, <sl2> )
##
## <L> is a semisimple lie algebra over Gaussian rationals or SqrtField
## and <sl2> is an SL2-triple in <L> of the form [f,h,e] with ef=h, he=2e, hf=-2f;
## this function returns the Z-graded carrier algebra of <sl2> normalised by
## CartanSubalgebra(<L>).
##
InstallGlobalFunction( RegularCarrierAlgebraOfSL2Triple, function( L, sl2 )
local cr, K0, H0, grading;
 
   if not HasCartanDecomposition(L) then 
      Error("no Cartan decomposition attached"); 
   fi;
   grading := [Basis(CartanDecomposition(L).K),Basis(CartanDecomposition(L).P)];
   cr      := corelg.carrZm( L, grading, sl2[3] );
   K0      := SubalgebraNC( L, cr.g0 );
   H0      := Intersection( CartanSubalgebra(L), K0 );
   if LieNormalizer(K0,H0) = H0 then
      return cr;
   else
      return corelg.carrierAlgebraBySL2Triple( L, grading, sl2 );
   fi;

end);




##################################################################################
##################################################################################








##################################################################################
#
#  Added nil orbs for outer aut... modification of SLA function

corelg.nil_orbs_outer:= function( L, gr0, gr1, gr2 )

     # Here L is a simple graded Lie algebra; gr0 a basis of the
     # elts of degree 0, gr1 of degree 1, and gr2 of degree -1.
     # We find the nilpotent G_0-orbits in g_1.
     # We *do not* assume that the given CSA of L is also a CSA of g_0.

     local F, g0, s, r, HL, Hs, R, Ci, hL, hl, C, rank, posRv_L, posR_L,
           posR, i, j, sums, fundR, inds, tr, h_candidates, BH, W, h, 
           c_h, ph, stb, v, w, is_rep, h0, wr, Omega, good_h, g1, g2, h_mats1,
           h_mats2, mat, sl2s, id1, id2, V, e, f, bb, ef, found, good, co, x, 
           C_h0, sp, sp0, y, b, bas, c, Cs, B, Rs, nas, b0, ranks, in_weylch,
           charact, k, sol, info;

     F:= LeftActingDomain(L);

     g0:= SubalgebraNC( L, gr0, "basis" );

     s:= LieDerivedSubalgebra( g0 );
     r:= LieCentre(g0);

     HL:= CartanSubalgebra(L);
     Hs:= Intersection( s, HL );
     SetCartanSubalgebra( s, Hs );

     R:= RootSystem(L);
     Ci:= CartanMatrix( R )^-1;
     hL:= ChevalleyBasis(L)[3];

     hl:= List( NilpotentOrbits(L), x -> Ci*WeightedDynkinDiagram(x) );
     for i in [1..Length(hl)] do
         if hl[i] = 0*hl[i] then
            Unbind( hl[i] );
         fi;
     od;
     hl:= Filtered( hl, x -> IsBound(x) );

     C:= CartanMatrix( R );
     rank:= Length(C);

     Rs:= RootsystemOfCartanSubalgebra(s);
     Cs:= CartanMatrix( Rs );
     ranks:= Length( Cs );

     bas:= ShallowCopy( CanonicalGenerators(Rs)[3] );
     Append( bas, BasisVectors( Basis(r) ) );
     b0:= Basis( VectorSpace( F, bas ), bas );

     in_weylch:= function( h )

          local cf, u;

          u:= h*hL;
          if not u in g0 then return false; fi;
          #if not corelg.eltInSubspace(L,BasisVectors(Basis(g0)),u) then return false; fi;
          cf:= Coefficients( b0, u ){[1..ranks]};
          if ForAll( Cs*cf, x -> x >= 0 ) then
             return true;
          else
             return false;
          fi;

     end;

     charact:= function( h )

          local cf;

          cf:= Coefficients( b0, h ){[1..ranks]};
          return Cs*cf;

     end;

     h_candidates:= SLAfcts.loop_W( C, hl, in_weylch );
     
     info:= "Constructed ";
     Append( info, String(Length(h_candidates)) );
     Append( info, " Cartan elements to be checked.");

     Info(InfoSLA,2,info);

     # now we need to compute sl_2 triples wrt the h-s found...

     Omega:= [0..Dimension(L)];
     good_h:= [ ];

     g1:= Basis( SubspaceNC( L, gr1 ), gr1 );
     g2:= Basis( SubspaceNC( L, gr2 ), gr2 );

     # the matrices of hL[i] acting on g1
     h_mats1:= [ ];
     for h0 in bas do
         mat:= [ ];
         for i in [1..Length(g1)] do
             Add( mat, Coefficients( g1, h0*g1[i] ) );
         od;
         Add( h_mats1, mat );
     od;

     # those of wrt g2...
     h_mats2:= [ ];
     for h0 in bas do
         mat:= [ ];
         for i in [1..Length(g1)] do
             Add( mat, Coefficients( g2, h0*g2[i] ) );
         od;
         Add( h_mats2, mat );
     od;

     sl2s:= [ ];
     id1:= IdentityMat( Length(g1) );
     id2:= IdentityMat( Length(g2) );
     for h in h_candidates do

         c_h:= Coefficients( b0, h*hL );

         mat:= c_h*h_mats1;
         mat:= mat - 2*id1;
         V:= NullspaceMat( mat );
         e:= List( V, v -> v*gr1 );

         mat:= c_h*h_mats2;
         mat:= mat + 2*id2;
         V:= NullspaceMat( mat );
         f:= List( V, v -> v*gr2 );

         # check whether h0 in [e,f]....
         bb:= [ ];
         for x in e do
             for y in f do
                 Add( bb, x*y );
             od;
         od;
         ef:= SubspaceNC( L, bb );

         h0:= h*hL;

         if h0 in ef then  #otherwise we can just discard h...
         #if corelg.eltInSubspace(L,BasisVectors(Basis(ef)),h0) then
            found:= false;
            good:= false;
            while not found do

                co:= List( e, x -> Random(Omega) );
                x:= co*e;
                sp:= SubspaceNC( L, List( f, y -> x*y) );

                if Dimension(sp) = Length(e) and h0 in sp then
                #if Dimension(sp) = Length(e) and 
                #   corelg.eltInSubspace(L,BasisVectors(Basis(sp)),h0) then 

                   # look for a nice one...
                   for i in [1..Length(co)] do
                       k:= 0;
                       found:= false;
                       while not found do
                           co[i]:= k;
                           x:= co*e;
                           sp:= SubspaceNC( L, List( f, y -> x*y) );

                           if Dimension(sp) = Length(e) and h0 in sp then
                           #if Dimension(sp) = Length(e) and 
                           #   corelg.eltInSubspace(L,BasisVectors(Basis(sp)),h0) then 
                              found:= true;
                           else
                              k:= k+1;
                           fi;
                       od;
                   od;

                   mat:= List( f, u -> Coefficients( Basis(sp), x*u ) );
                   sol:= SolutionMat( mat, Coefficients( Basis(sp), h0 ) );

                   Add( good_h, h0 );
                   Add( sl2s, [sol*f,h0,x] );

                   found:= true;

                else
                   C_h0:= LieCentralizer( g0, SubalgebraNC( g0, [h0] ) );
                   sp0:= SubspaceNC( L, List( Basis(C_h0), y -> y*x ) );
                   if Dimension(sp0) = Length(e) then
                      found:= true;
                      good:= false;
                   fi;
                fi;
      
            od;

         fi;
     od;

     return rec( hs:= good_h, sl2:= sl2s, chars:= List( good_h, charact ) );

end;





################################################################################









