###########################################################################
#
# this file contains the entry corelg.SuperLie, which is used in forms.gi
#
#
#


# the next function is from the library; we correct the order of simple
# roots of F4...

corelg.simliealg:= function( type, n, F )

    local T,               # The table of the Lie algebra constructed.
          i,j,k,l,         # Loop variables.
          lst,             # A list.
          R,               # Positive roots
          cc,              # List of coefficients.
          lenR,            # length of 'R'
          Rij,             # The sum of two roots from 'R'.
          eps,             # The so-called "epsilon"-function.
          epsmat,          # A matrix used to calculate the eps-function.
          dim,             # The dimension of the Lie algebra.
          C,               # Cartan matrix
          L,               # Lie algebra, result
          vectors,         # vectors spanning a Cartan subalgebra
          CSA,             # List of indices of the basis vectors of a Cartan
                           # subalgebra.
          e,
          inds,            # List of indices.
          r,r1,r2,         # Roots.
          roots,           # List of roots.
          primes,          # List of lists of corresponding roots.
          B,               # Basis of a vector space.
          cfs,             # List of coefficient lists.
          d,               # Order of the diagram automorphism.
          found,           # Boolean.
          a,
          q,
          perm,            # Permutation representing the diagram automorphism.
          shorts,
          posR,            # Positive roots.
          CartanMatrixToPositiveRoots; # Function for determining the
                                       # positive roots.


    CartanMatrixToPositiveRoots:= function( C )

        local   rank,  posr,  ready,  ind,  le,  i,  a,  j,  ej,  r,  b,
                q;

        rank:= Length( C );

        # `posr' will be a list of the positive roots. We start with the
        # simple roots, which are simply unit vectors.

        posr:= IdentityMat( rank );

        ready:= false;
        ind:= 1;
        le:= rank;
        while ind <= le  do

            # We loop over those elements of `posR' that have been found in
            # the previous round, i.e., those at positions ranging from
            # `ind' to `le'.

            le:= Length( posr );
            for i in [ind..le] do
                a:= posr[i];

                # We determine whether a+ej is a root (where ej is the j-th
                # simple root.
                for j in [1..rank] do
                    ej:= posr[j];

                    # We determine the maximum number `r' such that a-r*ej is
                    # a root.
                    r:= -1;
                    b:= ShallowCopy( a );
                    while b in posr do
                        b:= b-ej;
                        r:=r+1;
                    od;
                    q:= r-LinearCombination( TransposedMat( C )[j], a );
                    if q>0 and (not a+ej in posr ) then
                        Add( posr, a+ej );
                    fi;
                od;
            od;
            ind:= le+1;
            le:= Length( posr );
        od;

        return posr;
    end;


    # The following function is the so-called epsilon function.
    eps:= function( a, b, epm )
        local rk;

        rk:= Length( epm );
        return Product( [1..rk],i ->
                       Product( [1..rk], j ->
                               epm[i][j] ^ ( a[i]*b[j] ) ) );
    end;

    if type in [ "A", "D", "E" ] then

        # We are in the simply-laced case. Here we construct the root
        # system and the matrix of the epsilon function. Then we can
        # fill the multiplication table directly.

        C:= 2*IdentityMat( n );
        if type = "A" then
            for i in [1..n-1] do
                C[i][i+1]:= -1;
                C[i+1][i]:= -1;
            od;
        elif type = "D" then
            if n < 4 then
                Error("<n> must be >= 4");
            fi;
            for i in [1..n-2] do
                C[i][i+1]:= -1;
                C[i+1][i]:= -1;
            od;
            C[n-2][n]:=-1;
            C[n][n-2]:= -1;
        else

            C:= [
                 [ 2, 0, -1, 0, 0, 0, 0, 0 ], [ 0, 2, 0, -1, 0, 0, 0, 0 ],
                 [ -1, 0, 2, -1, 0, 0, 0, 0 ], [ 0, -1, -1, 2, -1, 0, 0, 0 ],
                 [ 0, 0, 0, -1, 2, -1, 0, 0 ], [ 0, 0, 0, 0, -1, 2, -1, 0 ],
                 [ 0, 0, 0, 0, 0, -1, 2, -1 ], [ 0, 0, 0, 0, 0, 0, -1, 2 ] ];

            if n = 6 then
                C:= C{ [ 1 .. 6 ] }{ [ 1 .. 6 ] };
            elif n = 7 then
                C:= C{ [ 1 .. 7 ] }{ [ 1 .. 7 ] };
            elif n < 6 or 8 < n then
                Error( "<n> must be one of 6, 7, 8" );
            fi;
        fi;
        R:= CartanMatrixToPositiveRoots( C );


        # We conctruct `epsmat', which satisfies
        #                  /
        #                 |-1 if i=j,
        #  epsmat[i][j] = |-1 if i and j are connected, and i>j
        #                 | 1 if i and j are not connected or i<j.
        #                  \
        # (where `connected' means connected in the Dynkin diagram.

        epsmat:= [];
        for i in [ 1 .. n ] do
            epsmat[i]:= [];
            for j in [ 1 .. i-1 ] do
                epsmat[i][j]:= 1;
            od;
            epsmat[i][i]:= -1;
            for j in [ i+1 .. n ] do
                epsmat[i][j]:= (-1)^C[i][j];
            od;
        od;

        lenR:= Length( R );
        dim:= 2*lenR + n;

        posR:= List( R, r -> Zero(F)*r );

        # Initialize the s.c. table
        T:= EmptySCTable( dim, Zero(F), "antisymmetric" );

        lst:= [ 1 .. n ] + 2 * lenR;

        for i in [1..lenR] do
            for j in [i..lenR] do
                Rij:= R[i]+R[j];
                if Rij in R then
                    k:= Position(R,Rij);
                    e:= eps(R[i],R[j],epsmat)*One(F);
                    SetEntrySCTable( T, i, j, [ e, k ] );
                    SetEntrySCTable( T, i+lenR, j+lenR, [ -e, k+lenR ] );
                fi;
                if i = j and T[i][j+lenR] = [[],[]] then
                    # We form the product x_{\alpha_i}*x_{-\alpha_i}, which
                    # will be an element of the Cartan subalgebra.

                    inds:= Filtered( [1..n], x -> R[i][x] <> 0 );
                    T[i][j+lenR]:= [ lst{inds}, R[i]{inds}*One(F) ];
                    T[j+lenR][i]:= [ lst{inds}, -R[i]{inds}*One(F) ];
                fi;
            od;
        od;
        for i in [1..lenR] do
            for j in [1..lenR] do
                Rij:= R[i]-R[j];
                if Rij in R then
                    k:= Position(R,Rij);
                    SetEntrySCTable( T, i, j+lenR,
                            [-One(F)*eps(R[i],-R[j],epsmat),k] );
                elif -Rij in R then
                    k:= Position(R,-Rij);
                    SetEntrySCTable( T, i, j+lenR,
                            [One(F)*eps(R[i],-R[j],epsmat),k+lenR] );
                fi;
            od;
            for j in [1..n] do

                # We take care of the comutation relations of the form
                # [h_j,x_{\beta_i}]= < \beta_i, \alpha_j > x_{\beta_i}.
                cc:= LinearCombination( R[i], C[j] );
                if cc <> 0*cc then

                    posR[i][j]:= One(F)*cc;

                    T[2*lenR+j][i]:=[[i],[One(F)*cc]];
                    T[i][2*lenR+j]:=[[i],[-One(F)*cc]];
                    T[2*lenR+j][i+lenR]:=[[i+lenR],[-One(F)*cc]];
                    T[i+lenR][2*lenR+j]:=[[i+lenR],[One(F)*cc]];
                fi;
            od;
        od;

        L:= LieAlgebraByStructureConstants( F, T );

        # A Cartan subalgebra is spanned by the last 'n' basis elements.
        CSA:= [ dim-n+1 .. dim ];
        vectors:= BasisVectors( CanonicalBasis( L ) ){ CSA };
        SetCartanSubalgebra( L, SubalgebraNC( L, vectors, "basis" ) );
        SetIsRestrictedLieAlgebra( L, Characteristic( F ) > 0 );

    elif type in [ "B", "C", "F", "G" ] then

        # Now we are in the non simply laced case. In each case we construct
        # a simply laced root system, which has a diagram automorphism.
        # We take an epsilon function which is invariant under the diagram
        # automorphism. Furthermore, the permutation `perm' will represent
        # the diagram aotomorphism as acting on the roots (so that
        # Permuted( r, perm ) is the result of applying the diagram
        # automorphism to the root r).

        if type = "B" then

            # In this case we construct D_{n+1}.
            if n <= 1 then
                Error( "<n> must be >= 2");
            fi;
            C:= 2*IdentityMat( n+1 );
            for i in [1..n-1] do
                C[i][i+1]:= -1;
                C[i+1][i]:= -1;
            od;
            C[n-1][n+1]:=-1;
            C[n+1][n-1]:= -1;
            R:= CartanMatrixToPositiveRoots( C );

            epsmat:= NullMat( n+1, n+1 ) + 1;
            for i in [ 1 .. n-1 ] do
                epsmat[i+1][i]:= -1;
                epsmat[i][i]:= -1;
            od;
            epsmat[n+1][n-1]:= -1;
            epsmat[n][n]:= -1;
            epsmat[n+1][n+1]:= -1;

            perm:= (n,n+1);
            d:= 2;

        elif type = "C" then

            # In this case we construct A_{2n-1}.
            if n < 2 then
                Error( "<n> must be >= 3");
            fi;
            C:= 2*IdentityMat( 2*n-1 );
            for i in [1..2*n-2] do
                C[i][i+1]:= -1;
                C[i+1][i]:= -1;
            od;
            R:= CartanMatrixToPositiveRoots( C );

            epsmat:= NullMat( 2*n-1, 2*n-1 ) + 1;
            for i in [ 1 .. n-1 ] do
                epsmat[i][i+1]:= -1;
                epsmat[i][i]:= -1;
            od;
            for i in [n..2*n-2] do
                epsmat[i+1][i]:= -1;
                epsmat[i][i]:= -1;
            od;
            epsmat[2*n-1][2*n-1]:= -1;

            perm:= ();
            for i in [1..n-1] do
                perm:= perm*(i,2*n-i);
            od;
            d:= 2;

        elif type = "F" then

            # In this case we construct E_6.
            if n <> 4 then
                Error( "<n> must be equal to 4");
            fi;

            C:= IdentityMat( 6 );
            C[1][3]:=-1; C[2][4]:=-1; C[3][4]:=-1; C[4][5]:=-1; C[5][6]:=-1;
            C:= C+TransposedMat( C );
            R:= CartanMatrixToPositiveRoots( C );

            epsmat:= NullMat( 6, 6 ) + 1;
            for i in [1..6] do epsmat[i][i]:= -1; od;
            epsmat[1][3]:=-1; epsmat[3][4]:=-1; epsmat[5][4]:=-1;
            epsmat[6][5]:=-1; epsmat[2][4]:=-1;

            perm:= (1,6)*(3,5);
            d:= 2;

        elif type = "G" then

            # In this case we conctruct D_4.
            if n <> 2 then
                Error( "<n> must be equal to 2");
            fi;

            C:= IdentityMat( 4 );
            C[1][2]:=-1; C[2][3]:=-1; C[2][4]:=-1;
            C:= C+TransposedMat( C );
            R:= CartanMatrixToPositiveRoots( C );

            epsmat:= NullMat( 4, 4 ) + 1;
            for i in [1..4] do epsmat[i][i]:= -1; od;
            epsmat[1][2]:=-1; epsmat[4][2]:=-1; epsmat[3][2]:=-1;

            perm:= (1,3,4);
            d:= 3;

        fi;

        # Now `roots' will be the list of positive roots of the resulting Lie
        # algebra. They are formed from the roots in `R' by applying the
        # diagram automorphism. If a r\in R is invariant under the
        # automorphism, then it is added to `roots' (and its prime is
        # the root itself). Otherwise we add \frac{1}{d}(r+\phi(r)+\cdots
        # + \phi^{d-1}(r)), where \phi is the diagram automorphism.
        # In this case the prime of the root are all \phi^i(r).

        if d = 2 then

            roots:= [ ];
            primes:= [ ];
            for r in R do
                r1:= Permuted( r, perm );
                if r = r1 then
                    Add( roots, r );
                    Add( primes, [ r ] );
                else
                    if not (r+r1)/2 in roots then
                        Add( roots, (r+r1)/2 );
                        Add( primes, [ r, r1 ] );
                    fi;
                fi;
            od;


           if type = "F" then
              roots:= Permuted( roots, (1,4,2) );
              primes:= Permuted( primes, (1,4,2) );
              C:= Permuted( C, (1,4,2) );
           fi;

            B:= Basis( VectorSpace( Rationals, roots{[1..n]} ),roots{[1..n]});
            cfs:= List( roots, x -> Coefficients( B, x ) );

        elif d = 3 then
            roots:= [ ];
            primes:= [ ];
            for r in R do
                r1:= Permuted( r, perm );
                if r = r1 then
                    Add( roots, r );
                    Add( primes, [ r ] );
                else
                    r2:= (r+r1+Permuted(r1,perm))/3;
                    if not r2 in roots then
                        Add( roots, r2 );
                        Add( primes, [ r, r1, Permuted( r1, perm ) ] );
                    fi;
                fi;
            od;

            B:= Basis( VectorSpace( Rationals, roots{[1..n]} ),roots{[1..n]});
            cfs:= List( roots, x -> Coefficients( B, x ) );
        fi;

        # `shorts' will be a list of indices indicating where the
        # short simple roots are. The coefficients on those places
        # in `cfs' need to be divided by `d'.

        shorts:= Filtered( [1..n], ii -> Length( primes[ii] ) > 1 );
        for i in [1..Length(cfs)] do
            for j in shorts do
                cfs[i][j]:= cfs[i][j]/d;
            od;
        od;

        Append( R, -R );
        lenR:= Length( roots );
        dim:= 2*lenR + n;

        posR:= List( [1..lenR], ii -> List( [1..n], jj -> Zero( F ) ) );

        # Initialize the s.c. table
        T:= EmptySCTable( dim, Zero(F), "antisymmetric" );

        lst:= [ 1 .. n ] + 2 * lenR;

        for i in [1..lenR] do
            for j in [i..lenR] do
                Rij:= roots[i]+roots[j];
                if Rij in roots then

                    # We look for `r' in `primes[i]' and `r1' in `primes[j]'
                    # such that `r+r1' lies in `R'.
                    found:= false;
                    for k in [1..Length(primes[i])] do
                        if found then break; fi;
                        r:= primes[i][k];
                        for l in [1..Length(primes[j])] do
                            r1:= primes[j][l];
                            if r+r1 in R then
                                found := true; break;
                            fi;
                        od;
                    od;

                    # `q' will be the maximal integer such that `roots[i]-
                    # roots[j]' is a root.

                    k:= Position( roots, Rij );
                    q:=0; a:= roots[i] - roots[j];
                    while a in roots or -a in roots do
                        q:=q+1;
                        a:= a-roots[j];
                    od;

                    e:= eps(r,r1,epsmat)*(q+1)*One(F);
                    SetEntrySCTable( T, i, j, [ e, k ] );
                    SetEntrySCTable( T, i+lenR, j+lenR, [ -e, k+lenR ] );
                fi;
                if i = j and T[i][j+lenR] = [[],[]] then
                    # We form the product x_{\alpha_i}*x_{-\alpha_i}, which
                    # will be an element of the Cartan subalgebra.

                    inds:= Filtered( [1..n], x -> cfs[i][x] <> 0 );
                    if Length( primes[i] ) = 1 then
                        T[i][j+lenR]:= [ lst{inds}, cfs[i]{inds}*One(F) ];
                        T[j+lenR][i]:= [ lst{inds}, -cfs[i]{inds}*One(F) ];
                    else
                        T[i][j+lenR]:= [ lst{inds}, cfs[i]{inds}*d*One(F) ];
                        T[j+lenR][i]:= [ lst{inds}, -cfs[i]{inds}*d*One(F) ];
                    fi;
                fi;
            od;
        od;
        for i in [1..lenR] do
            for j in [1..lenR] do
                Rij:= roots[i]-roots[j];
                if Rij in roots then

                    found:= false;
                    for k in [1..Length(primes[i])] do
                        if found then break; fi;
                        r:= primes[i][k];
                        for l in [1..Length(primes[j])] do
                            r1:= primes[j][l];
                            if r-r1 in R then
                                found := true; break;
                            fi;
                        od;
                    od;

                    k:= Position( roots, Rij );
                    q:=0; a:= roots[i] + roots[j];
                    while a in roots or -a in roots do
                        q:=q+1;
                        a:= a+roots[j];
                    od;

                    SetEntrySCTable( T, i, j+lenR,
                            [-One(F)*(q+1)*eps(r,-r1,epsmat),k] );

                elif -Rij in roots then

                    found:= false;
                    for k in [1..Length(primes[i])] do
                        if found then break; fi;
                        r:= primes[i][k];
                        for l in [1..Length(primes[j])] do
                            r1:= primes[j][l];
                            if r-r1 in R then
                                found := true; break;
                            fi;
                        od;
                    od;

                    k:= Position( roots, -Rij );
                    q:=0; a:= roots[i] + roots[j];
                    while a in roots or -a in roots do
                        q:=q+1;
                        a:= a+roots[j];
                    od;
                    SetEntrySCTable( T, i, j+lenR,
                            [One(F)*(q+1)*eps(r,-r1,epsmat),k+lenR] );
                fi;
            od;
            for j in [1..n] do

                # Now we take care of the relations [h,x_{\beta}]....

                cc:= LinearCombination( roots[i], C[j] );
                if Length( primes[j] ) > 1 then
                    # i.e., `roots[j]' is "short".
                    cc:= d*cc;
                fi;

                if cc <> 0*cc then

                    posR[i][j]:= One(F)*cc;

                    T[2*lenR+j][i]:=[[i],[One(F)*cc]];
                    T[i][2*lenR+j]:=[[i],[-One(F)*cc]];
                    T[2*lenR+j][i+lenR]:=[[i+lenR],[-One(F)*cc]];
                    T[i+lenR][2*lenR+j]:=[[i+lenR],[One(F)*cc]];
                fi;
            od;
        od;

        L:= LieAlgebraByStructureConstants( F, T );

        # A Cartan subalgebra is spanned by the last 'n' basis elements.
        CSA:= [ dim-n+1 .. dim ];
        vectors:= BasisVectors( CanonicalBasis( L ) ){ CSA };
        SetCartanSubalgebra( L, SubalgebraNC( L, vectors, "basis" ) );
        SetIsRestrictedLieAlgebra( L, Characteristic( F ) > 0 );

    fi;

    R:= Objectify( NewType( NewFamily( "RootSystemFam", IsObject ),
                IsAttributeStoringRep and IsRootSystemFromLieAlgebra ),
                rec() );
    SetUnderlyingLieAlgebra( R, L );
    SetPositiveRoots( R, posR );
    SetNegativeRoots( R, -posR );
    SetSimpleSystem( R, posR{[1..n]} );
    SetCanonicalGenerators( R, [ CanonicalBasis( L ){[1..n]},
                                 CanonicalBasis( L ){[lenR+1..lenR+n]},
                                 vectors ] );
    SetPositiveRootVectors( R, CanonicalBasis(L){[1..lenR]} );
    SetNegativeRootVectors( R, CanonicalBasis(L){[lenR+1..2*lenR]} );
    SetChevalleyBasis( L, [ PositiveRootVectors( R ),
                            NegativeRootVectors( R ),
                            vectors ] );

    if not ( Characteristic( F ) in [ 2, 3 ] ) then

        C:= 2*IdentityMat( n );
        for i in [1..n] do
            for j in [1..n] do
                if i <> j then
                    q:= 0;
                    r:= posR[i]+posR[j];
                    while r in posR do
                        q:=q+1;
                        r:= r+posR[j];
                    od;
                    C[i][j]:= -q;
                fi;
            od;
        od;

        SetCartanMatrix( R, C );

        SetSemiSimpleType( L, Concatenation( type, String( n ) ) );
    fi;

    SetRootSystem( L, R );

    if Characteristic( F ) = 0 then
       SetIsSimpleAlgebra( L, true );
    fi;

    return L;


end;


corelg.SuperLie:=function(arg)


local T, L, Lo, To, R, P, S, s, V, B, p, TT, i, j, g, K, f, g0, g1, k, l, D, U,  KK1, KK2, KK3, KK4, BB, bb, N, CC, a, b, c, d, ww, K0, P0, H0, t, pos, makeCartInv, sp, rts, v, q, n, t0, h, F0, CartInt, fundr, allrts; 

q:= arg[1];
n:= arg[2];
t0:= arg[3];
h:= arg[4];
if Length(arg)=5 then
   F0:= arg[5];
else
   F0:= CF(4);
fi;

t:= ShallowCopy(t0);

Lo:=corelg.simliealg(q,n,Rationals);
To:= StructureConstantsTable(Basis(Lo));;
#T:=Sub3( a, n);
#L:=LieAlgebraByStructureConstants( Rationals,T);
R:= RootSystem(Lo);
P:= PositiveRoots(R);
S:= SimpleSystem(R);
#C:= ChevalleyBasis(L);
V:= VectorSpace(Rationals,  S);
B:= Basis(V, S);
s:= Length(S);;
p:= Length(P);;
D:=[];; #set of alpha such that X=x-x
U:=[];; #set of alpha such that X=x+x
TT:=EmptySCTable( 2*p+s,  Zero(F0),  "antisymmetric" );;
a:=0;
b:=0;
c:=0;
d:=0;

#############################################################################################################################################

#cerchiamo di sistemare il vettore dei\Phi(alpha), ovvero tenere a mente le permutazioni
#cioè in questo modo creo il vettore N con le radici permutate
#mentre il vettore P ha quelle normali.

K:=[];

for i in [1..s] do
	Add(K, P[Position(P,P[OnPoints(i,h)])] );
od;
#ho appena sistemato il vettore K con le permutazioni definite da h sul simple system

############################################################################################################################################

#resta da sistemare il resto delle radici positive

for i in [s+1..p] do 
	N:=P[1]*0;
	CC:=Coefficients(B,P[i]); #mi scrivo i coefficienti di P[i] rispetto al simple System
	for j in [1..s] do
		N:=N+CC[j]*P[Position(P,K[j])]; #Calcolo la permutaziona associata a questa radice
	od;
	Add(K, P[Position(P,N)] );
od;

#Vettore K con le immagini della prmutazione è a posto.

###############################################################################################################################################

#Sistemiamo il vettore dei segni, per ora l ho solo definito sul Simple System

for i in [s+1..p] do
	for j in [1..s] do
		if (P[i]-P[j] in P ) then
			f:=Position(P, P[i]-P[j]); #calcolo la posizione in P di P[i]-P[j]
			g0:=Position(P,K[j]); #calcolo la posizione in P della permutazione di P[j]
			g1:=Position(P, K[f]); #calcolo la posizione in P della permutazione di P[i]-P[j]
		Add(t, t[j]*t[f]*To[g0][g1][2][1]*(To[j][f][2][1]^(-1)));
		break;
		fi;
	od;
od;


##############################################################################################################################################

#sistemiamo i prodotti [H,X] e [H,Y]

for i in [1..s] do
	for j in [1..p] do
		if (Position(P,K[i]) > i) then #permutazione solo sugli h, quindi mi bastano i primi s
			if ( Position(P,K[j]) > Position(P,P[j])) then
a:=0;
b:=0;
c:=0;
d:=0;

if not ( IsEmpty(To[2*p+Position(P,P[i])][j][2]) =true) then  
	if not ( IsEmpty(To[2*p+Position(P,K[i])][j][2])=true ) then 

a:=a+To[2*p+Position(P,P[i])][j][2][1];
b:=b+To[2*p+Position(P,K[i])][j][2][1];
c:=c+To[2*p+Position(P,P[i])][Position(P,K[j])][2][1];
d:=d+To[2*p+Position(P,K[i])][Position(P,K[j])][2][1];

SetEntrySCTable( TT, 2*p+Position(P,P[i]), 	  	    j,[    (a+b+c+d)/2, 		  p+j ]);
SetEntrySCTable( TT, 2*p+Position(P,P[i]),   Position(P,K[j]),[    (a+b+c+d)/2, p+Position(P,K[j])]);
SetEntrySCTable( TT, 2*p+Position(P,P[i]), 		  p+j,[ -1*(a+b+c+d)/2, 		    j] );
SetEntrySCTable( TT, 2*p+Position(P,P[i]), p+Position(P,K[j]),[ -1*(a+b+c+d)/2, Position(P,K[j])]);

SetEntrySCTable( TT, 2*p+Position(P,K[i]), 	  	    j, [-1*(a-b-c+d)/2, p+Position(P,K[j]) ]);
SetEntrySCTable( TT, 2*p+Position(P,K[i]),   Position(P,K[j]), [   (a-b-c+d)/2, 		       p+j]);
SetEntrySCTable( TT, 2*p+Position(P,K[i]), 		  p+j, [   (a-b-c+d)/2, 	  Position(P,K[j]) ]);
SetEntrySCTable( TT, 2*p+Position(P,K[i]), p+Position(P,K[j]), [-1*(a-b-c+d)/2, 	   	       j ]);

	else

a:=a+To[2*p+Position(P,P[i])][j][2][1];
#b:=b+To[2*p+Position(P,K[i])][j][2][1];
#c:=c+To[2*p+Position(P,P[i])][Position(P,K[j])][2][1];
d:=d+To[2*p+Position(P,K[i])][Position(P,K[j])][2][1];

SetEntrySCTable( TT, 2*p+Position(P,P[i]), 	  	    j,[    (a+d)/2, 		  p+j ]);
SetEntrySCTable( TT, 2*p+Position(P,P[i]),   Position(P,K[j]),[    (a+d)/2, p+Position(P,K[j])]);
SetEntrySCTable( TT, 2*p+Position(P,P[i]), 		  p+j,[ -1*(a+d)/2, 		    j] );
SetEntrySCTable( TT, 2*p+Position(P,P[i]), p+Position(P,K[j]),[ -1*(a+d)/2, Position(P,K[j])]);

SetEntrySCTable( TT, 2*p+Position(P,K[i]), 	  	    j,[ -1*(a+d)/2, p+Position(P,K[j]) ]);
SetEntrySCTable( TT, 2*p+Position(P,K[i]),   Position(P,K[j]),[    (a+d)/2, 		       p+j]);
SetEntrySCTable( TT, 2*p+Position(P,K[i]), 		  p+j,[    (a+d)/2, 	  Position(P,K[j] )]);
SetEntrySCTable( TT, 2*p+Position(P,K[i]), p+Position(P,K[j]),[ -1*(a+d)/2, 	   	       j]) ;

	fi;
else
	if not ( IsEmpty(To[2*p+Position(P,K[i])][j][2])=true ) then 

#a:=a+To[2*p+Position(P,P[i])][j][2][1];
b:=b+To[2*p+Position(P,K[i])][j][2][1];
c:=c+To[2*p+Position(P,P[i])][Position(P,K[j])][2][1];
#d:=d+To[2*p+Position(P,K[i])][Position(P,K[j])][2][1];

SetEntrySCTable( TT, 2*p+Position(P,P[i]), 	  	    j,[    (b+c)/2, 		  p+j ]);
SetEntrySCTable( TT, 2*p+Position(P,P[i]),   Position(P,K[j]),[    (b+c)/2, p+Position(P,K[j])]);
SetEntrySCTable( TT, 2*p+Position(P,P[i]), 		  p+j,[ -1*(b+c)/2, 		    j ]);
SetEntrySCTable( TT, 2*p+Position(P,P[i]), p+Position(P,K[j]),[ -1*(b+c)/2, Position(P,K[j])]);

SetEntrySCTable( TT, 2*p+Position(P,K[i]), 	  	    j,[ -1*(-b-c)/2, p+Position(P,K[j]) ]);
SetEntrySCTable( TT, 2*p+Position(P,K[i]),   Position(P,K[j]),[    (-b-c)/2, 		       p+j]);
SetEntrySCTable( TT, 2*p+Position(P,K[i]), 		  p+j,[    (-b-c)/2, 	  Position(P,K[j] )]);
SetEntrySCTable( TT, 2*p+Position(P,K[i]), p+Position(P,K[j]),[ -1*(-b-c)/2, 	   	       j]) ;

	fi;

fi;

			else
				if ( Position(P,K[j]) = Position(P,P[j])) then

a:=0;
b:=0;
c:=0;
d:=0;

if not ( IsEmpty(To[2*p+Position(P,P[i])][j][2]) =true) then  
	if not ( IsEmpty(To[2*p+Position(P,K[i])][j][2])=true ) then 

a:=a+To[2*p+Position(P,P[i])][j][2][1];
c:=c+To[2*p+Position(P,K[i])][j][2][1];

SetEntrySCTable( TT, 2*p+Position(P,P[i]), 	  	    j,[    (a+c),   		   p+j]);
SetEntrySCTable( TT, 2*p+Position(P,P[i]), 		  p+j,[ -1*(a+c), 		    j ]);
	else
a:=a+To[2*p+Position(P,P[i])][j][2][1];

SetEntrySCTable( TT, 2*p+Position(P,P[i]), 	  	    j,[    a,   		   p+j]);
SetEntrySCTable( TT, 2*p+Position(P,P[i]), 		  p+j,[ -1*a, 		    j ]);
	fi;
else
	if not ( IsEmpty(To[2*p+Position(P,K[i])][j][2])=true ) then 

c:=c+To[2*p+Position(P,K[i])][j][2][1];

SetEntrySCTable( TT, 2*p+Position(P,P[i]), 	  	    j,[    c,   		   p+j]);
SetEntrySCTable( TT, 2*p+Position(P,P[i]), 		  p+j,[ -1*c, 		    j ]);
	fi;
fi;

				fi;
			fi;
		else
			if (Position(P,K[i]) = i) then #permutazione solo sugli h, quindi mi bastano i primi s
				if ( Position(P,K[j]) > Position(P,P[j])) then
a:=0;
b:=0;
c:=0;
d:=0;

if not ( IsEmpty(To[2*p+Position(P,P[i])][j][2]) =true) then  
	if not ( IsEmpty(To[2*p+Position(P,P[i])][Position(P,K[j])][2])=true ) then 

a:=a+To[2*p+Position(P,P[i])][j][2][1];
c:=c+To[2*p+Position(P,P[i])][Position(P,K[j])][2][1];

SetEntrySCTable( TT, 2*p+Position(P,P[i]), 	  	    j,[    (a+c),   		   p+j]);
SetEntrySCTable( TT, 2*p+Position(P,P[i]),   Position(P,K[j]),[    (a+c),   p+Position(P,K[j])]);
SetEntrySCTable( TT, 2*p+Position(P,P[i]), 		  p+j,[ -1*(a+c), 		    j ]);
SetEntrySCTable( TT, 2*p+Position(P,P[i]), p+Position(P,K[j]),[ -1*(a+c),     Position(P,K[j])]);

	else

a:=a+To[2*p+Position(P,P[i])][j][2][1];

SetEntrySCTable( TT, 2*p+Position(P,P[i]), 	  	    j,[    a, 		  p+j ]);
SetEntrySCTable( TT, 2*p+Position(P,P[i]),   Position(P,K[j]),[    a, p+Position(P,K[j])]);
SetEntrySCTable( TT, 2*p+Position(P,P[i]), 		  p+j,[ -1*a, 		    j ]);
SetEntrySCTable( TT, 2*p+Position(P,P[i]), p+Position(P,K[j]),[ -1*a, Position(P,K[j])]);
	fi;

else
	if not ( IsEmpty(To[2*p+Position(P,P[i])][Position(P,K[j])][2]) =true ) then

c:=c+To[2*p+Position(P,P[i])][Position(P,K[j])][2][1];

SetEntrySCTable( TT, 2*p+Position(P,P[i]), 	  	    j,[    c, 		  p+j ]);
SetEntrySCTable( TT, 2*p+Position(P,P[i]),   Position(P,K[j]),[    c, p+Position(P,K[j])]);
SetEntrySCTable( TT, 2*p+Position(P,P[i]), 		  p+j,[ -1*c, 		    j ]);
SetEntrySCTable( TT, 2*p+Position(P,P[i]), p+Position(P,K[j]),[ -1*c, Position(P,K[j])]);

	fi;
fi;

				else
					if ( Position(P,K[j]) = Position(P,P[j])) then
a:=0;

if not ( IsEmpty(To[2*p+Position(P,P[i])][j][2]) =true) then  

a:=a+To[2*p+Position(P,P[i])][j][2][1];

SetEntrySCTable( TT, 2*p+Position(P,P[i]),	j,[2*a, p+j] );
SetEntrySCTable( TT, 2*p+Position(P,P[i]),	p+j,[-2*a, j] );
fi;

					fi;
				fi;
			fi;
		fi;
	od;
od;


#######################################################################################################################################
#####################################################################################################################################

#sistemiamo i prodotti [X,Y] con indici uguali


for i in [1..p] do #alpha varia in P
	if ( Position(P,K[i]) > Position(P, P[i]) ) then #alpha <phi_alpha

KK1:=[]; #++ nessun problema
KK2:=[]; #+- nessun problema
KK3:=[]; #-+ occhio ai prodotti misti
KK4:=[]; #-- occhio ai prodotti misti

if (t[i]=1) then

if (P[i]+K[i] in P) then
if not ( IsEmpty( To[i][Position(P,K[i])][2] )= true ) then
if ( t[Position(P, P[i]+K[i])] = 1 ) then

if not ( (+1*To[i][Position(P, K[i])][2][1]+1*To[Position(P, K[i])][i][2][1]) = 0 ) then
Add(KK1, (+1*To[i][Position(P, K[i])][2][1]+1*To[Position(P, K[i])][i][2][1])/2 ); #aggiungo il valore del prodotto relativo all'indice j
Add(KK1, p+Position(P,P[i]+K[i])); #aggiungo la relativa posizione  di j nella nuova base
Add(KK4, +1*(+1*To[i][Position(P, K[i])][2][1]+1*To[Position(P, K[i])][i][2][1])/2 ); #aggiungo il valore del prodotto relativo all'indice j
Add(KK4, p+Position(P,P[i]+K[i])); 
fi;

else

if not ( (-1*To[i][Position(P, K[i])][2][1]+1*To[Position(P, K[i])][i][2][1]) = 0 ) then
Add(KK2, (-1*To[i][Position(P, K[i])][2][1]+1*To[Position(P, K[i])][i][2][1])/2 ); #aggiungo il valore del prodotto relativo all'indice j
Add(KK2, p+Position(P,P[i]+K[i])); #aggiungo la relativa posizione  di j nella nuova base
Add(KK3, (+1*To[i][Position(P, K[i])][2][1]-1*To[Position(P, K[i])][i][2][1])/2 ); #aggiungo il valore del prodotto relativo all'indice j
Add(KK3, p+Position(P,P[i]+K[i]));  
fi;

fi; fi; fi;

if not (IsEmpty( To[i][p+i][2]) = true) then

for j in [1..s] do 
	for k in [1..Length(To[i][p+i][2])] do
		if (2*p+j = To[i][p+i][1][k] ) then

			if ( OnPoints(j,h) > j ) then
				if (2*p+OnPoints(j,h) in To[i][p+i][1] ) then

ww:=Position(To[i][p+i][1],2*p+OnPoints(j,h));

if not (To[i][p+i][2][k] = To[i][p+i][2][ww] ) then
Add(KK1, +2*(To[i][p+i][2][k] + To[i][p+i][2][ww] )); #aggiungo il valore del prodotto relativo all'indice j
Add(KK1, 2*p+j); #aggiungo la relativa posizione  di j nella nuova base
Add(KK2, -2*(To[i][p+i][2][k] - To[i][p+i][2][ww] )); #aggiungo il valore del prodotto relativo all'indice j
Add(KK2, 2*p+Position(P,K[j])); #aggiungo la relativa posizione  di j nella nuova base
Add(KK3, -2*(To[i][p+i][2][k] - To[i][p+i][2][ww] )); #aggiungo il valore del prodotto relativo all'indice j
Add(KK3, 2*p+Position(P,K[j])); 
Add(KK4, -2*(To[i][p+i][2][k] + To[i][p+i][2][ww] )); #aggiungo il valore del prodotto relativo all'indice j
Add(KK4, 2*p+j);
else
Add(KK1, +2*(To[i][p+i][2][k] + To[i][p+i][2][ww] )); #aggiungo il valore del prodotto relativo all'indice j
Add(KK1, 2*p+j); #aggiungo la relativa posizione  di j nella nuova base
Add(KK4, -2*(To[i][p+i][2][k] + To[i][p+i][2][ww] )); #aggiungo il valore del prodotto relativo all'indice j
Add(KK4, 2*p+j);
fi; 	
				else
Add(KK1, 2*To[i][p+i][2][k] ); #aggiungo il valore del prodotto relativo all'indice j
Add(KK1, 2*p+j); #aggiungo la relativa posizione  di j nella nuova base
Add(KK2, -2*To[i][p+i][2][k] ); #aggiungo il valore del prodotto relativo all'indice j
Add(KK2, 2*p+Position(P,K[j])); #aggiungo la relativa posizione  di j nella nuova base
Add(KK3, -2*To[i][p+i][2][k] ); #aggiungo il valore del prodotto relativo all'indice j
Add(KK3, 2*p+Position(P,K[j])); 
Add(KK4, -2*To[i][p+i][2][k] ); #aggiungo il valore del prodotto relativo all'indice j
Add(KK4, 2*p+j); 				
				fi;					
			else
					if ( OnPoints(j,h) = j ) then
Add(KK1, +2*To[i][p+i][2][k] ); #aggiungo il valore del prodotto relativo all'indice j
Add(KK1, 2*p+j); #aggiungo la relativa posizione  di j nella nuova base
Add(KK4, -2*To[i][p+i][2][k] ); #aggiungo il valore del prodotto relativo all'indice j
Add(KK4, 2*p+j); 
					fi;
			fi;

fi; od; od; fi;


else #t[i]=-1

if not (IsEmpty( To[i][p+i][2]) = true) then

for j in [1..s] do 
	for k in [1..Length(To[i][p+i][2])] do
		if (2*p+j = To[i][p+i][1][k] ) then

			if ( OnPoints(j,h) > j ) then
				if (2*p+OnPoints(j,h) in To[i][p+i][1] ) then
ww:=Position(To[i][p+i][1],2*p+OnPoints(j,h));

if not (To[i][p+i][2][k] = To[i][p+i][2][ww] ) then
Add(KK1, +2*(To[i][p+i][2][k] + To[i][p+i][2][ww] )); #aggiungo il valore del prodotto relativo all'indice j
Add(KK1, 2*p+j); #aggiungo la relativa posizione  di j nella nuova base
Add(KK2, -2*(To[i][p+i][2][k] - To[i][p+i][2][ww] )); #aggiungo il valore del prodotto relativo all'indice j
Add(KK2, 2*p+Position(P,K[j])); #aggiungo la relativa posizione  di j nella nuova base
Add(KK3, -2*(To[i][p+i][2][k] - To[i][p+i][2][ww] )); #aggiungo il valore del prodotto relativo all'indice j
Add(KK3, 2*p+Position(P,K[j])); 
Add(KK4, -2*(To[i][p+i][2][k] + To[i][p+i][2][ww] )); #aggiungo il valore del prodotto relativo all'indice j
Add(KK4, 2*p+j);
else
Add(KK1, +2*(To[i][p+i][2][k] + To[i][p+i][2][ww] )); #aggiungo il valore del prodotto relativo all'indice j
Add(KK1, 2*p+j); #aggiungo la relativa posizione  di j nella nuova base
Add(KK4, -2*(To[i][p+i][2][k] + To[i][p+i][2][ww] )); #aggiungo il valore del prodotto relativo all'indice j
Add(KK4, 2*p+j);
fi; 		
				else
Add(KK1, 2*To[i][p+i][2][k] ); #aggiungo il valore del prodotto relativo all'indice j
Add(KK1, 2*p+j); #aggiungo la relativa posizione  di j nella nuova base
Add(KK2, -2*To[i][p+i][2][k] ); #aggiungo il valore del prodotto relativo all'indice j
Add(KK2, 2*p+Position(P,K[j])); #aggiungo la relativa posizione  di j nella nuova base
Add(KK3, -2*To[i][p+i][2][k] ); #aggiungo il valore del prodotto relativo all'indice j
Add(KK3, 2*p+Position(P,K[j])); 
Add(KK4, -2*To[i][p+i][2][k] ); #aggiungo il valore del prodotto relativo all'indice j
Add(KK4, 2*p+j); 				
				fi;					
			else
					if ( OnPoints(j,h) = j ) then
Add(KK1, +2*To[i][p+i][2][k] ); #aggiungo il valore del prodotto relativo all'indice j
Add(KK1, 2*p+j); #aggiungo la relativa posizione  di j nella nuova base
Add(KK4, -2*To[i][p+i][2][k] ); #aggiungo il valore del prodotto relativo all'indice j
Add(KK4, 2*p+j); 
					fi;
			fi;

fi; od; od; fi;

if (P[i]+K[i] in P) then
if not ( IsEmpty(To[i][Position(P,K[i])][2])= true ) then
if ( t[Position(P, P[i]+K[i])] = 1 ) then

if not ( (-1*To[i][Position(P, K[i])][2][1]-1*To[Position(P, K[i])][i][2][1]) = 0 ) then
Add(KK1, (-1*To[i][Position(P, K[i])][2][1]-1*To[Position(P, K[i])][i][2][1])/2 ); #aggiungo il valore del prodotto relativo all'indice j
Add(KK1, p+Position(P,P[i]+K[i])); #aggiungo la relativa posizione  di j nella nuova base
Add(KK4, (-1*To[i][Position(P, K[i])][2][1]-1*To[Position(P, K[i])][i][2][1])/2 ); #aggiungo il valore del prodotto relativo all'indice j
Add(KK4, p+Position(P,P[i]+K[i]));  
fi;

else
if not ( (+1*To[i][Position(P, K[i])][2][1]-1*To[Position(P, K[i])][i][2][1]) = 0 ) then
Add(KK2, (+1*To[i][Position(P, K[i])][2][1]-1*To[Position(P, K[i])][i][2][1])/2 ); #aggiungo il valore del prodotto relativo all'indice j
Add(KK2, p+Position(P,P[i]+K[i])); #aggiungo la relativa posizione  di j nella nuova base
Add(KK3, (-1*To[i][Position(P, K[i])][2][1]+1*To[Position(P, K[i])][i][2][1])/2 ); #aggiungo il valore del prodotto relativo all'indice j
Add(KK3, p+Position(P,P[i]+K[i])); 
fi;

fi; fi; fi;

fi; 

SetEntrySCTable( TT, i, p+i , KK1);
SetEntrySCTable( TT, i, p+Position(P,K[i]) , KK2);
SetEntrySCTable( TT, Position(P, K[i]), p+i , KK3);
SetEntrySCTable( TT, Position(P, K[i]), p+Position(P, K[i]) , KK4);

	else
		if ( Position(P,K[i]) = Position(P, P[i]) ) then #alpha <phi_alpha

KK1:=[]; #++ nessun problema

for j in [1..s] do 
	for k in [1..Length(To[i][p+i][2])] do
		if (2*p+j = To[i][p+i][1][k] ) then

if (t[i]=1) then

			if ( OnPoints(j,h) > j ) then
				if (2*p+OnPoints(j,h) in To[i][p+i][1] ) then
ww:=Position(To[i][p+i][1],2*p+OnPoints(j,h));
Add(KK1, +4*(To[i][p+i][2][k] + To[i][p+i][2][ww] ) ); #aggiungo il valore del prodotto relativo all'indice j
Add(KK1, 2*p+j); #aggiungo la relativa posizione  di j nella nuova base
				else
Add(KK1, +4*To[i][p+i][2][k] ); #aggiungo il valore del prodotto relativo all'indice j
Add(KK1, 2*p+j); #aggiungo la relativa posizione  di j nella nuova base
				fi;					
			else
					if ( OnPoints(j,h) = j ) then
Add(KK1, +4*To[i][p+i][2][k] ); #aggiungo il valore del prodotto relativo all'indice j
Add(KK1, 2*p+j); #aggiungo la relativa posizione  di j nella nuova base
					fi;
			fi;

else #t[i]=-1

#KK1:=[]; #-- occhio ai prodotti misti
			if ( OnPoints(j,h) > j ) then
				if (2*p+OnPoints(j,h) in To[i][p+i][1] ) then
ww:=Position(To[i][p+i][1],2*p+OnPoints(j,h));
Add(KK1, -4*(To[i][p+i][2][k] + To[i][p+i][2][ww] ) ); #aggiungo il valore del prodotto relativo all'indice j
Add(KK1, 2*p+j); #aggiungo la relativa posizione  di j nella nuova base	
				else
Add(KK1, -4*To[i][p+i][2][k] ); #aggiungo il valore del prodotto relativo all'indice j
Add(KK1, 2*p+j); 				
				fi;					
			else
					if ( OnPoints(j,h) = j ) then
Add(KK1, -4*To[i][p+i][2][k] ); #aggiungo il valore del prodotto relativo all'indice j
Add(KK1, 2*p+j); 
					fi;
			fi;
fi;

		fi;
	od;
od;

SetEntrySCTable( TT, i, p+i , KK1);

		fi;
	fi;
od;
#####################################################################################################################



##############################################################################################################################################

#Prodotti [X,X] con indici differenti


for i in [1..p] do #variazione alpha per gli X
	for j in [i+1..p] do #variazione beta >alpha per gli Y
		if ( Position(P, K[i]) > Position(P,P[i]) ) then #alpha< phi_alpha (ragiono con gli X)
			if ( p+Position(P, K[j]) > p+Position(P,P[j]) ) then # beta< phi_beta (ragiono con gli Y)

KK1:=[]; # i p+j
KK2:=[]; # i p+jbar
KK3:=[]; # ibar p+j
KK4:=[]; # ibar p+jbar

if (t[i]=1) then
if (t[j]=1) then

if (P[i]+P[j] in P ) then #qui so che aplha+beta < phi_aplha+phi_beta quindi sono nel caso sicuro
if not (IsEmpty(To[Position(P, P[i])][Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]+K[i]) > p+Position(P,P[j]+P[i])) then

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, Position(P, P[i]+P[j]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, Position(P, K[i]+K[j]) );
Add(KK3, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, Position(P, K[i]+K[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, Position(P, P[i]+P[j]) );
else
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, Position(P, P[i]+P[j]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, Position(P, K[i]+K[j]) );
Add(KK3, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, Position(P, K[i]+K[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, Position(P, P[i]+P[j]) );
fi; 
else
	if (p+Position(P,K[j]+K[i]) = p+Position(P,P[j]+P[i])) then

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, Position(P, P[i]+P[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, Position(P, P[i]+P[j]) );
else
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, Position(P, K[i]+K[j]) );
Add(KK3, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, Position(P, K[i]+K[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, Position(P, K[j]+K[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, Position(P, P[i]+P[j]) );
Add(KK3, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, Position(P, P[i]+P[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, Position(P,K[j]+K[i]) );
else
Add(KK1, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, Position(P, K[j]+K[i]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, Position(P, P[i]+P[j]) );
Add(KK3, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, Position(P, P[i]+P[j]) );
Add(KK4, +1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, Position(P, K[j]+K[i]) );
fi; 
	fi; 

fi; 
fi; fi;

#######
############
###############

if (P[j]-P[i] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][p+Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]-K[i]) > p+Position(P,P[j]-P[i])) then

if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK1, Position(P, P[j]-P[i]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK2, Position(P, K[j]-K[i]) );
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK3, Position(P, K[j]-K[i]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, Position(P, P[j]-P[i]) );
else
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK1, Position(P, P[j]-P[i]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK2, Position(P, K[j]-K[i]) );
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK3, Position(P, K[j]-K[i]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, Position(P, P[j]-P[i]) );
fi; 

else
	if (p+Position(P,K[j]-K[i]) = p+Position(P,P[j]-P[i])) then

if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK1, Position(P, P[j]-P[i]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, Position(P, P[j]-P[i]) );
else
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK2, Position(P, K[j]-K[i]) );
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK3, Position(P, K[j]-K[i]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))

if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK1, Position(P, K[j]-K[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK2, Position(P, P[j]-P[i]) );
Add(KK3, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK3, Position(P, P[j]-P[i]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, Position(P, K[j]-K[i]) );
else
Add(KK1, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK1, Position(P, K[j]-K[i]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK2, Position(P, P[j]-P[i]) );
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK3, Position(P, P[j]-P[i]) );
Add(KK4, +1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, Position(P, K[j]-K[i]) );
fi; 
	fi; 

fi; 
fi; fi;


######
##########
#############

if (P[i]+K[j] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][Position(P,K[j])][2])=true) then
if (p+Position(P,K[i]+P[j]) > p+Position(P,P[i]+K[j])) then

if (t[Position(P,P[i]+K[j] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, Position(P, P[i]+K[j]) );
Add(KK2, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, Position(P, K[i]+P[j]) );
Add(KK3, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, Position(P, K[i]+P[j]) );
Add(KK4, -1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, Position(P, P[i]+K[j]) );
else
Add(KK1, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, Position(P, P[i]+K[j]) );
Add(KK2, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, Position(P, K[i]+P[j]) );
Add(KK3, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, Position(P, K[i]+P[j]) );
Add(KK4, -1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, Position(P, P[i]+K[j]) );
fi; 
else
	if (p+Position(P,K[i]+P[j]) = p+Position(P,P[i]+K[j])) then

if (t[Position(P,K[j]+P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, Position(P, P[i]+K[j]) );
Add(KK4, -1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, Position(P, P[i]+K[j]) );
else
Add(KK2, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, Position(P, K[i]+P[j]) );
Add(KK3, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, Position(P, K[i]+P[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))

if (t[Position(P,P[i]+K[j] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, Position(P, K[i]+P[j]) );
Add(KK2, -1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, Position(P, P[i]+K[j]) );
Add(KK3, -1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, Position(P, P[i]+K[j]) );
Add(KK4, -1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, Position(P, K[i]+P[j]) );
else
Add(KK1, -1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, Position(P, K[i]+P[j]) );
Add(KK2, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, Position(P, P[i]+K[j]) );
Add(KK3, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, Position(P, P[i]+K[j]) );
Add(KK4, +1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, Position(P, K[i]+P[j]) );
fi; 
	fi; 

fi; 
fi; fi;

######
##########
#############


if (K[j]-P[i] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][p+Position(P,K[j])][2])=true) then
if (p+Position(P,P[j]-K[i]) > p+Position(P,K[j]-P[i])) then

if (t[Position(P,K[j]-P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK1, Position(P, K[j]-P[i]) );
Add(KK2, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK2, Position(P, P[j]-K[i]) );
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK3, Position(P, P[j]-K[i]) );
Add(KK4, -1*(-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, Position(P, K[j]-P[i]) );
else
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK1, Position(P, K[j]-P[i]) );
Add(KK2, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK2, Position(P, P[j]-K[i]) );
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK3, Position(P, P[j]-K[i]) );
Add(KK4, -1*(-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, Position(P, K[j]-P[i]) );
fi; 

else
	if (p+Position(P,P[j]-K[i]) = p+Position(P,K[j]-P[i])) then

if (t[Position(P,K[j]-P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK1, Position(P, K[j]-P[i]) );
Add(KK4, -1*(-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, Position(P, K[j]-P[i]) );
else
Add(KK2, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK2, Position(P, P[j]-K[i]) );
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK3, Position(P, P[j]-K[i]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))

if (t[Position(P,K[j]-P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK1, Position(P, P[j]-K[i]) );
Add(KK2, -1*(-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK2, Position(P, K[j]-P[i]) );
Add(KK3, -1*(+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK3, Position(P, K[j]-P[i]) );
Add(KK4, -1*(-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, Position(P, P[j]-K[i]) );
else
Add(KK1, -1*(+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK1, Position(P, P[j]-K[i]) );
Add(KK2, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK2, Position(P, K[j]-P[i]) );
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK3, Position(P, K[j]-P[i]) );
Add(KK4, +1*(-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, Position(P, P[j]-K[i]) );
fi; 
	fi; 

fi; 
fi; fi;

###########
##################
####################



else # ovvero t[j]=-1

###############
#####################
#########################

if (P[i]+P[j] in P ) then #qui so che aplha+beta < phi_aplha+phi_beta quindi sono nel caso sicuro
if not (IsEmpty(To[Position(P, P[i])][Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]+K[i]) > p+Position(P,P[j]+P[i])) then

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, Position(P, P[i]+P[j]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, Position(P, K[i]+K[j]) );
Add(KK3, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, Position(P, K[i]+K[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, Position(P, P[i]+P[j]) );
else
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, Position(P, P[i]+P[j]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, Position(P, K[i]+K[j]) );
Add(KK3, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, Position(P, K[i]+K[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, Position(P, P[i]+P[j]) );
fi; 
else
	if (p+Position(P,K[j]+K[i]) = p+Position(P,P[j]+P[i])) then

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, Position(P, P[i]+P[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, Position(P, P[i]+P[j]) );
else
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, Position(P, K[i]+K[j]) );
Add(KK3, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, Position(P, K[i]+K[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, Position(P, K[j]+K[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, Position(P, P[i]+P[j]) );
Add(KK3, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, Position(P, P[i]+P[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, Position(P,K[j]+K[i]) );
else
Add(KK1, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, Position(P, K[j]+K[i]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, Position(P, P[i]+P[j]) );
Add(KK3, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, Position(P, P[i]+P[j]) );
Add(KK4, +1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, Position(P, K[j]+K[i]) );
fi; 
	fi; 

fi; 
fi; fi;

#######
############
###############

if (P[j]-P[i] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][p+Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]-K[i]) > p+Position(P,P[j]-P[i])) then

if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK1, Position(P, P[j]-P[i]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK2, Position(P, K[j]-K[i]) );
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK3, Position(P, K[j]-K[i]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, Position(P, P[j]-P[i]) );
else
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK1, Position(P, P[j]-P[i]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK2, Position(P, K[j]-K[i]) );
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK3, Position(P, K[j]-K[i]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, Position(P, P[j]-P[i]) );
fi; 
else
	if (p+Position(P,K[j]-K[i]) = p+Position(P,P[j]-P[i])) then

if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK1, Position(P, P[j]-P[i]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, Position(P, P[j]-P[i]) );
else
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK2, Position(P, K[j]-K[i]) );
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK3, Position(P, K[j]-K[i]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))

if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK1, Position(P, K[j]-K[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK2, Position(P, P[j]-P[i]) );
Add(KK3, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK3, Position(P, P[j]-P[i]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, Position(P, K[j]-K[i]) );
else
Add(KK1, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK1, Position(P, K[j]-K[i]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK2, Position(P, P[j]-P[i]) );
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK3, Position(P, P[j]-P[i]) );
Add(KK4, +1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, Position(P, K[j]-K[i]) );
fi; 
	fi; 

fi; 
fi; fi;

######
##########
#############

if (P[i]+K[j] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][Position(P,K[j])][2])=true) then
if (p+Position(P,K[i]+P[j]) > p+Position(P,P[i]+K[j])) then

if (t[Position(P,P[i]+K[j] )] = 1) then
Add(KK1, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, Position(P, P[i]+K[j]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, Position(P, K[i]+P[j]) );
Add(KK3, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, Position(P, K[i]+P[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, Position(P, P[i]+K[j]) );
else
Add(KK1, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, Position(P, P[i]+K[j]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, Position(P, K[i]+P[j]) );
Add(KK3, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, Position(P, K[i]+P[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, Position(P, P[i]+K[j]) );
fi; 
else
	if (p+Position(P,K[i]+P[j]) = p+Position(P,P[i]+K[j])) then

if (t[Position(P,K[j]+P[i] )] = 1) then
Add(KK1, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, Position(P, P[i]+K[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, Position(P, P[i]+K[j]) );
else
Add(KK2, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, Position(P, K[i]+P[j]) );
Add(KK3, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, Position(P, K[i]+P[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))

if (t[Position(P,P[i]+K[j] )] = 1) then
Add(KK1, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, Position(P, K[i]+P[j]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, Position(P, P[i]+K[j]) );
Add(KK3, -1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, Position(P, P[i]+K[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, Position(P, K[i]+P[j]) );
else
Add(KK1, -1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, Position(P, K[i]+P[j]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, Position(P, P[i]+K[j]) );
Add(KK3, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, Position(P, P[i]+K[j]) );
Add(KK4, +1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, Position(P, K[i]+P[j]) );
fi; 
	fi; 

fi; 
fi; fi;

######
##########
#############


if (K[j]-P[i] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][p+Position(P,K[j])][2])=true) then
if (p+Position(P,P[j]-K[i]) > p+Position(P,K[j]-P[i])) then

if (t[Position(P,K[j]-P[i] )] = 1) then
Add(KK1, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK1, Position(P, K[j]-P[i]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK2, Position(P, P[j]-K[i]) );
Add(KK3, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK3, Position(P, P[j]-K[i]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, Position(P, K[j]-P[i]) );
else
Add(KK1, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK1, Position(P, K[j]-P[i]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK2, Position(P, P[j]-K[i]) );
Add(KK3, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK3, Position(P, P[j]-K[i]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, Position(P, K[j]-P[i]) );
fi; 
else
	if (p+Position(P,P[j]-K[i]) = p+Position(P,K[j]-P[i])) then

if (t[Position(P,K[j]-P[i] )] = 1) then
Add(KK1, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK1, Position(P, K[j]-P[i]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, Position(P, K[j]-P[i]) );
else
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK2, Position(P, P[j]-K[i]) );
Add(KK3, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK3, Position(P, P[j]-K[i]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))

if (t[Position(P,K[j]-P[i] )] = 1) then
Add(KK1, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK1, Position(P, P[j]-K[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK2, Position(P, K[j]-P[i]) );
Add(KK3, -1*(-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK3, Position(P, K[j]-P[i]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, Position(P, P[j]-K[i]) );

else

Add(KK1, -1*(-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK1, Position(P, P[j]-K[i]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK2, Position(P, K[j]-P[i]) );
Add(KK3, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK3, Position(P, K[j]-P[i]) );
Add(KK4, +1*(+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, Position(P, P[j]-K[i]) );

fi; 
	fi; 

fi; 
fi; fi;




fi;

##############
#####################
########################à

else #t[i]=-1

if (t[j]=1) then


###############
#####################
#########################

if (P[i]+P[j] in P ) then #qui so che aplha+beta < phi_aplha+phi_beta quindi sono nel caso sicuro
if not (IsEmpty(To[Position(P, P[i])][Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]+K[i]) > p+Position(P,P[j]+P[i])) then

if (t[Position(P,P[j]+P[i] )] = 1) then

Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, Position(P, P[i]+P[j]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, Position(P, K[i]+K[j]) );
Add(KK3, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, Position(P, K[i]+K[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, Position(P, P[i]+P[j]) );

else

Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, Position(P, P[i]+P[j]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, Position(P, K[i]+K[j]) );
Add(KK3, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, Position(P, K[i]+K[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, Position(P, P[i]+P[j]) );

fi; 

else
	if (p+Position(P,K[j]+K[i]) = p+Position(P,P[j]+P[i])) then

if (t[Position(P,P[j]+P[i] )] = 1) then

Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, Position(P, P[i]+P[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, Position(P, P[i]+P[j]) );

else

Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, Position(P, K[i]+K[j]) );
Add(KK3, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, Position(P, K[i]+K[j]) );

fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))

if (t[Position(P,P[j]+P[i] )] = 1) then

Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, Position(P, K[j]+K[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, Position(P, P[i]+P[j]) );
Add(KK3, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, Position(P, P[i]+P[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, Position(P,K[j]+K[i]) );

else

Add(KK1, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, Position(P, K[j]+K[i]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, Position(P, P[i]+P[j]) );
Add(KK3, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, Position(P, P[i]+P[j]) );
Add(KK4, +1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, Position(P, K[j]+K[i]) );

fi; 
	fi; 

fi; 
fi; fi;





#######
############
###############





if (P[j]-P[i] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][p+Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]-K[i]) > p+Position(P,P[j]-P[i])) then

if (t[Position(P,P[j]-P[i] )] = 1) then

Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK1, Position(P, P[j]-P[i]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK2, Position(P, K[j]-K[i]) );
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK3, Position(P, K[j]-K[i]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, Position(P, P[j]-P[i]) );

else

Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK1, Position(P, P[j]-P[i]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK2, Position(P, K[j]-K[i]) );
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK3, Position(P, K[j]-K[i]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, Position(P, P[j]-P[i]) );

fi; 

else
	if (p+Position(P,K[j]-K[i]) = p+Position(P,P[j]-P[i])) then

if (t[Position(P,P[j]-P[i] )] = 1) then

Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK1, Position(P, P[j]-P[i]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, Position(P, P[j]-P[i]) );

else

Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK2, Position(P, K[j]-K[i]) );
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK3, Position(P, K[j]-K[i]) );

fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))

if (t[Position(P,P[j]-P[i] )] = 1) then

Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK1, Position(P, K[j]-K[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK2, Position(P, P[j]-P[i]) );
Add(KK3, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK3, Position(P, P[j]-P[i]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, Position(P, K[j]-K[i]) );

else

Add(KK1, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK1, Position(P, K[j]-K[i]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK2, Position(P, P[j]-P[i]) );
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK3, Position(P, P[j]-P[i]) );
Add(KK4, +1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, Position(P, K[j]-K[i]) );

fi; 
	fi; 

fi; 
fi; fi;


######
##########
#############


if (P[i]+K[j] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][Position(P,K[j])][2])=true) then
if (p+Position(P,K[i]+P[j]) > p+Position(P,P[i]+K[j])) then

if (t[Position(P,P[i]+K[j] )] = 1) then

Add(KK1, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, Position(P, P[i]+K[j]) );
Add(KK2, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, Position(P, K[i]+P[j]) );
Add(KK3, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, Position(P, K[i]+P[j]) );
Add(KK4, -1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, Position(P, P[i]+K[j]) );

else

Add(KK1, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, Position(P, P[i]+K[j]) );
Add(KK2, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, Position(P, K[i]+P[j]) );
Add(KK3, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, Position(P, K[i]+P[j]) );
Add(KK4, -1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, Position(P, P[i]+K[j]) );

fi; 

else
	if (p+Position(P,K[i]+P[j]) = p+Position(P,P[i]+K[j])) then

if (t[Position(P,K[j]+P[i] )] = 1) then

Add(KK1, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, Position(P, P[i]+K[j]) );
Add(KK4, -1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, Position(P, P[i]+K[j]) );

else

Add(KK2, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, Position(P, K[i]+P[j]) );
Add(KK3, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, Position(P, K[i]+P[j]) );

fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))

if (t[Position(P,P[i]+K[j] )] = 1) then

Add(KK1, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, Position(P, K[i]+P[j]) );
Add(KK2, -1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, Position(P, P[i]+K[j]) );
Add(KK3, -1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, Position(P, P[i]+K[j]) );
Add(KK4, -1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, Position(P, K[i]+P[j]) );

else

Add(KK1, -1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, Position(P, K[i]+P[j]) );
Add(KK2, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, Position(P, P[i]+K[j]) );
Add(KK3, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, Position(P, P[i]+K[j]) );
Add(KK4, +1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, Position(P, K[i]+P[j]) );

fi; 
	fi; 

fi; 
fi; fi;

######
##########
#############


if (K[j]-P[i] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][p+Position(P,K[j])][2])=true) then
if (p+Position(P,P[j]-K[i]) > p+Position(P,K[j]-P[i])) then

if (t[Position(P,K[j]-P[i] )] = 1) then

Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK1, Position(P, K[j]-P[i]) );
Add(KK2, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK2, Position(P, P[j]-K[i]) );
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK3, Position(P, P[j]-K[i]) );
Add(KK4, -1*(-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, Position(P, K[j]-P[i]) );

else

Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK1, Position(P, K[j]-P[i]) );
Add(KK2, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK2, Position(P, P[j]-K[i]) );
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK3, Position(P, P[j]-K[i]) );
Add(KK4, -1*(-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, Position(P, K[j]-P[i]) );

fi; 

else
	if (p+Position(P,P[j]-K[i]) = p+Position(P,K[j]-P[i])) then

if (t[Position(P,K[j]-P[i] )] = 1) then

Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK1, Position(P, K[j]-P[i]) );
Add(KK4, -1*(-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, Position(P, K[j]-P[i]) );

else

Add(KK2, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK2, Position(P, P[j]-K[i]) );
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK3, Position(P, P[j]-K[i]) );

fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))

if (t[Position(P,K[j]-P[i] )] = 1) then

Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK1, Position(P, P[j]-K[i]) );
Add(KK2, -1*(-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK2, Position(P, K[j]-P[i]) );
Add(KK3, -1*(+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK3, Position(P, K[j]-P[i]) );
Add(KK4, -1*(-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, Position(P, P[j]-K[i]) );

else

Add(KK1, -1*(+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK1, Position(P, P[j]-K[i]) );
Add(KK2, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK2, Position(P, K[j]-P[i]) );
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK3, Position(P, K[j]-P[i]) );
Add(KK4, +1*(-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, Position(P, P[j]-K[i]) );

fi; 
	fi; 

fi; 
fi; fi;



else #t[j]=-1


###############
#####################
#########################

if (P[i]+P[j] in P ) then #qui so che aplha+beta < phi_aplha+phi_beta quindi sono nel caso sicuro
if not (IsEmpty(To[Position(P, P[i])][Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]+K[i]) > p+Position(P,P[j]+P[i])) then

if (t[Position(P,P[j]+P[i] )] = 1) then

Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, Position(P, P[i]+P[j]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, Position(P, K[i]+K[j]) );
Add(KK3, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, Position(P, K[i]+K[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, Position(P, P[i]+P[j]) );

else

Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, Position(P, P[i]+P[j]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, Position(P, K[i]+K[j]) );
Add(KK3, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, Position(P, K[i]+K[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, Position(P, P[i]+P[j]) );

fi; 

else
	if (p+Position(P,K[j]+K[i]) = p+Position(P,P[j]+P[i])) then

if (t[Position(P,P[j]+P[i] )] = 1) then

Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, Position(P, P[i]+P[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, Position(P, P[i]+P[j]) );

else

Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, Position(P, K[i]+K[j]) );
Add(KK3, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, Position(P, K[i]+K[j]) );

fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))

if (t[Position(P,P[j]+P[i] )] = 1) then

Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, Position(P, K[j]+K[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, Position(P, P[i]+P[j]) );
Add(KK3, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, Position(P, P[i]+P[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, Position(P,K[j]+K[i]) );

else

Add(KK1, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, Position(P, K[j]+K[i]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, Position(P, P[i]+P[j]) );
Add(KK3, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, Position(P, P[i]+P[j]) );
Add(KK4, +1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, Position(P, K[j]+K[i]) );

fi; 
	fi; 

fi; 
fi; fi;





#######
############
###############


if (P[j]-P[i] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][p+Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]-K[i]) > p+Position(P,P[j]-P[i])) then

if (t[Position(P,P[j]-P[i] )] = 1) then

Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK1, Position(P, P[j]-P[i]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK2, Position(P, K[j]-K[i]) );
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK3, Position(P, K[j]-K[i]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, Position(P, P[j]-P[i]) );

else

Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK1, Position(P, P[j]-P[i]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK2, Position(P, K[j]-K[i]) );
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK3, Position(P, K[j]-K[i]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, Position(P, P[j]-P[i]) );

fi; 

else
	if (p+Position(P,K[j]-K[i]) = p+Position(P,P[j]-P[i])) then

if (t[Position(P,P[j]-P[i] )] = 1) then

Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK1, Position(P, P[j]-P[i]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, Position(P, P[j]-P[i]) );

else

Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK2, Position(P, K[j]-K[i]) );
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK3, Position(P, K[j]-K[i]) );

fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))

if (t[Position(P,P[j]-P[i] )] = 1) then

Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK1, Position(P, K[j]-K[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK2, Position(P, P[j]-P[i]) );
Add(KK3, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK3, Position(P, P[j]-P[i]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, Position(P, K[j]-K[i]) );

else

Add(KK1, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK1, Position(P, K[j]-K[i]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK2, Position(P, P[j]-P[i]) );
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK3, Position(P, P[j]-P[i]) );
Add(KK4, +1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, Position(P, K[j]-K[i]) );

fi; 
	fi; 

fi; 
fi; fi;


######
##########
#############


if (P[i]+K[j] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][Position(P,K[j])][2])=true) then
if (p+Position(P,K[i]+P[j]) > p+Position(P,P[i]+K[j])) then

if (t[Position(P,P[i]+K[j] )] = 1) then

Add(KK1, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, Position(P, P[i]+K[j]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, Position(P, K[i]+P[j]) );
Add(KK3, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, Position(P, K[i]+P[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, Position(P, P[i]+K[j]) );

else

Add(KK1, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, Position(P, P[i]+K[j]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, Position(P, K[i]+P[j]) );
Add(KK3, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, Position(P, K[i]+P[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, Position(P, P[i]+K[j]) );

fi; 

else
	if (p+Position(P,K[i]+P[j]) = p+Position(P,P[i]+K[j])) then

if (t[Position(P,K[j]+P[i] )] = 1) then

Add(KK1, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, Position(P, P[i]+K[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, Position(P, P[i]+K[j]) );

else

Add(KK2, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, Position(P, K[i]+P[j]) );
Add(KK3, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, Position(P, K[i]+P[j]) );

fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))

if (t[Position(P,P[i]+K[j] )] = 1) then

Add(KK1, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, Position(P, K[i]+P[j]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, Position(P, P[i]+K[j]) );
Add(KK3, -1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, Position(P, P[i]+K[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, Position(P, K[i]+P[j]) );

else

Add(KK1, -1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, Position(P, K[i]+P[j]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, Position(P, P[i]+K[j]) );
Add(KK3, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, Position(P, P[i]+K[j]) );
Add(KK4, +1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, Position(P, K[i]+P[j]) );

fi; 
	fi; 

fi; 
fi; fi;

######
##########
#############


if (K[j]-P[i] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][p+Position(P,K[j])][2])=true) then
if (p+Position(P,P[j]-K[i]) > p+Position(P,K[j]-P[i])) then

if (t[Position(P,K[j]-P[i] )] = 1) then

Add(KK1, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK1, Position(P, K[j]-P[i]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK2, Position(P, P[j]-K[i]) );
Add(KK3, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK3, Position(P, P[j]-K[i]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, Position(P, K[j]-P[i]) );

else

Add(KK1, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK1, Position(P, K[j]-P[i]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK2, Position(P, P[j]-K[i]) );
Add(KK3, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK3, Position(P, P[j]-K[i]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, Position(P, K[j]-P[i]) );

fi; 

else
	if (p+Position(P,P[j]-K[i]) = p+Position(P,K[j]-P[i])) then

if (t[Position(P,K[j]-P[i] )] = 1) then

Add(KK1, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK1, Position(P, K[j]-P[i]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, Position(P, K[j]-P[i]) );

else

Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK2, Position(P, P[j]-K[i]) );
Add(KK3, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK3, Position(P, P[j]-K[i]) );

fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))

if (t[Position(P,K[j]-P[i] )] = 1) then

Add(KK1, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK1, Position(P, P[j]-K[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK2, Position(P, K[j]-P[i]) );
Add(KK3, -1*(-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK3, Position(P, K[j]-P[i]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, Position(P, P[j]-K[i]) );

else

Add(KK1, -1*(-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK1, Position(P, P[j]-K[i]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK2, Position(P, K[j]-P[i]) );
Add(KK3, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK3, Position(P, K[j]-P[i]) );
Add(KK4, +1*(+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, Position(P, P[j]-K[i]) );

fi; 
	fi; 

fi; 
fi; fi;


fi; fi;

SetEntrySCTable( TT, i,  j , KK1);
SetEntrySCTable( TT, i,  Position(P, K[j]), KK2);
SetEntrySCTable( TT, Position(P,K[i]), j, KK3);
SetEntrySCTable( TT, Position(P,K[i]) ,Position(P, K[j]), KK4);

else

	if ( p+Position(P, K[j]) = p+Position(P,P[j]) ) then #beta = phi_beta

KK1:=[];
KK2:=[];

if (t[i]=1) then
if (t[j]=1) then

if (P[i]+P[j] in P ) then #qui so che aplha+beta < phi_aplha+phi_beta quindi sono nel caso sicuro
if not (IsEmpty(To[Position(P, P[i])][Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]+K[i]) > p+Position(P,P[j]+P[i])) then

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, P[i]+P[j]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, K[i]+K[j]) );
else
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, P[i]+P[j]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, K[i]+K[j]) );
fi; 

else
	if (p+Position(P,K[j]+K[i]) = p+Position(P,P[j]+P[i])) then
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, P[i]+P[j]) );
else
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, K[i]+K[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, K[j]+K[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, P[i]+P[j]) );
else
Add(KK1, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, K[j]+K[i]) );
Add(KK2, +1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, P[i]+P[j]) );
fi; 
	fi; 
fi; 
fi; fi;

#######
############
###############

if (P[j]-P[i] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][p+Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]-K[i]) > p+Position(P,P[j]-P[i])) then

if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, Position(P, P[j]-P[i]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, Position(P, K[j]-K[i]) );
else
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, Position(P, P[j]-P[i]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, Position(P, K[j]-K[i]) );
fi; 
else
	if (p+Position(P,K[j]-K[i]) = p+Position(P,P[j]-P[i])) then
if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, Position(P, P[j]-P[i]) );
else
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, Position(P, K[j]-K[i]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, Position(P, K[j]-K[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, Position(P, P[j]-P[i]) );
else
Add(KK1, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, Position(P, K[j]-K[i]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, Position(P, P[j]-P[i]) );
fi; 
	fi; 
fi; 
fi; fi;

###########
##################
####################

else #t[j]=-1

if (P[i]+P[j] in P ) then #qui so che aplha+beta < phi_aplha+phi_beta quindi sono nel caso sicuro
if not (IsEmpty(To[Position(P, P[i])][Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]+K[i]) > p+Position(P,P[j]+P[i])) then
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, K[j]+K[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, P[i]+P[j]) );
else
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, K[j]+K[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, P[i]+P[j]) );
fi; 
else
	if (p+Position(P,K[j]+K[i]) = p+Position(P,P[j]+P[i])) then
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, P[i]+P[j]) );
else
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, K[j]+K[i]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, K[j]+K[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, P[i]+P[j]) );
else
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, P[i]+P[j]) );
Add(KK2, +1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, K[j]+K[i]) );
fi; 
	fi; 
fi; 
fi; fi;

#######
############
###############

if (P[j]-P[i] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][p+Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]-K[i]) > p+Position(P,P[j]-P[i])) then
if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, Position(P, K[j]-K[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, Position(P, P[j]-P[i]) );
else
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, Position(P, K[j]-K[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, Position(P, P[j]-P[i]) );
fi; 
else
	if (p+Position(P,K[j]-K[i]) = p+Position(P,P[j]-P[i])) then
if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, Position(P, P[j]-P[i]) );
else
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, Position(P, K[j]-K[i]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, Position(P, P[j]-P[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, Position(P, K[j]-K[i]) );
else
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, Position(P, P[j]-P[i]) );
Add(KK2, +1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, Position(P, K[j]-K[i]) );
fi; 
	fi; 
fi; 
fi; fi;

fi;

else #t[i]=-1
if (t[j]=1) then

if (P[i]+P[j] in P ) then #qui so che aplha+beta < phi_aplha+phi_beta quindi sono nel caso sicuro
if not (IsEmpty(To[Position(P, P[i])][Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]+K[i]) > p+Position(P,P[j]+P[i])) then
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, P[i]+P[j]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, K[i]+K[j]) );
else
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, P[i]+P[j]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, K[i]+K[j]) );
fi; 
else
	if (p+Position(P,K[j]+K[i]) = p+Position(P,P[j]+P[i])) then
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, P[i]+P[j]) );
else
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, K[i]+K[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, K[j]+K[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, P[i]+P[j]) );
else
Add(KK1, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, K[j]+K[i]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, P[i]+P[j]) );
fi; 
	fi; 
fi; 
fi; fi;

#######
############
###############

if (P[j]-P[i] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][p+Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]-K[i]) > p+Position(P,P[j]-P[i])) then
if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, Position(P, P[j]-P[i]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, Position(P, K[j]-K[i]) );
else
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, Position(P, P[j]-P[i]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, Position(P, K[j]-K[i]) );
fi; 
else
	if (p+Position(P,K[j]-K[i]) = p+Position(P,P[j]-P[i])) then
if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, Position(P, P[j]-P[i]) );
else
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, Position(P, K[j]-K[i]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, Position(P, K[j]-K[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, Position(P, P[j]-P[i]) );
else
Add(KK1, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, Position(P, K[j]-K[i]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, Position(P, P[j]-P[i]) );
fi; 
	fi; 
fi; 
fi; fi;

###########
##################
####################

else #t[i]=t[j]=-1

if (P[i]+P[j] in P ) then #qui so che aplha+beta < phi_aplha+phi_beta quindi sono nel caso sicuro
if not (IsEmpty(To[Position(P, P[i])][Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]+K[i]) > p+Position(P,P[j]+P[i])) then
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, K[i]+K[j]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, P[j]+P[i]) );
else
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, K[i]+K[j]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, P[j]+P[i]) );
fi; 
else
	if (p+Position(P,K[j]+K[i]) = p+Position(P,P[j]+P[i])) then
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, P[j]+P[i]) );
else
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, K[i]+K[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, P[j]+P[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, K[i]+K[j]) );
else
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, P[j]+P[i]) );
Add(KK2, +1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, K[i]+K[j]) );
fi; 
	fi; 
fi; 
fi; fi;

#######
############
###############

if (P[j]-P[i] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][p+Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]-K[i]) > p+Position(P,P[j]-P[i])) then
if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, Position(P, K[j]-K[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, Position(P, P[j]-P[i]) );
else
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, Position(P, K[j]-K[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, Position(P, P[j]-P[i]) );
fi; 
else
	if (p+Position(P,K[j]-K[i]) = p+Position(P,P[j]-P[i])) then
if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, Position(P, P[j]-P[i]) );
else
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, Position(P, K[j]-K[i]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, Position(P, P[j]-P[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, Position(P, K[j]-K[i]) );
else
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, Position(P, P[j]-P[i]) );
Add(KK2, +1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, Position(P, K[j]-K[i]) );
fi; 
	fi; 
fi; 
fi; fi;

fi; fi;

SetEntrySCTable( TT, i,  j , KK1);
SetEntrySCTable( TT, Position(P,K[i]), j, KK2);

fi;
fi;

else
	if ( p+Position(P, K[i]) = p+Position(P,P[i]) ) then #alpha< phi_alpha
			if ( p+Position(P, K[j]) > p+Position(P,P[j]) ) then #beta< phi_beta

KK1:=[];
KK2:=[];

if (t[i]=1) then
if (t[j]=1) then

if (P[i]+P[j] in P ) then #qui so che aplha+beta < phi_aplha+phi_beta quindi sono nel caso sicuro
if not (IsEmpty(To[Position(P, P[i])][Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]+K[i]) > p+Position(P,P[j]+P[i])) then
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, P[i]+P[j]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, K[i]+K[j]) );
else
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, P[i]+P[j]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, K[i]+K[j]) );
fi; 
else
	if (p+Position(P,K[j]+K[i]) = p+Position(P,P[j]+P[i])) then
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, P[i]+P[j]) );
else
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, K[i]+K[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, K[j]+K[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, P[i]+P[j]) );
else
Add(KK1, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, K[j]+K[i]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, P[i]+P[j]) );
fi; 
	fi; 
fi; 
fi; fi;

#######
############
###############

if (P[j]-P[i] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][p+Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]-K[i]) > p+Position(P,P[j]-P[i])) then
if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, Position(P, P[j]-P[i]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, Position(P, K[j]-K[i]) );
else
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, Position(P, P[j]-P[i]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, Position(P, K[j]-K[i]) );
fi; 
else
	if (p+Position(P,K[j]-K[i]) = p+Position(P,P[j]-P[i])) then
if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, Position(P, P[j]-P[i]) );
else
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, Position(P, K[j]-K[i]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, Position(P, K[j]-K[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, Position(P, P[j]-P[i]) );
else
Add(KK1, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, Position(P, K[j]-K[i]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, Position(P, P[j]-P[i]) );
fi; 
	fi; 
fi; 
fi; fi;


###########
##################
####################

else #t[j]=-1

if (P[i]+P[j] in P ) then #qui so che aplha+beta < phi_aplha+phi_beta quindi sono nel caso sicuro
if not (IsEmpty(To[Position(P, P[i])][Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]+K[i]) > p+Position(P,P[j]+P[i])) then
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, P[i]+P[j]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, K[i]+K[j]) );
else
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, P[i]+P[j]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, K[i]+K[j]) );
fi; 
else
	if (p+Position(P,K[j]+K[i]) = p+Position(P,P[j]+P[i])) then
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, P[i]+P[j]) );
else
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, K[i]+K[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, K[j]+K[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, P[i]+P[j]) );
else
Add(KK1, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, K[j]+K[i]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, P[i]+P[j]) );
fi; 
	fi; 
fi; 
fi; fi;

#######
############
###############

if (P[j]-P[i] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][p+Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]-K[i]) > p+Position(P,P[j]-P[i])) then
if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, Position(P, P[j]-P[i]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, Position(P, K[j]-K[i]) );
else
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, Position(P, P[j]-P[i]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, Position(P, K[j]-K[i]) );
fi; 
else
	if (p+Position(P,K[j]-K[i]) = p+Position(P,P[j]-P[i])) then
if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, Position(P, P[j]-P[i]) );
else
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, Position(P, K[j]-K[i]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, Position(P, K[j]-K[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, Position(P, P[j]-P[i]) );
else
Add(KK1, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, Position(P, K[j]-K[i]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, Position(P, P[j]-P[i]) );
fi; 
	fi; 
fi; 
fi; fi;

###########
##################
####################
fi;

else #t[i]=1
if (t[j]=1) then

if (P[i]+P[j] in P ) then #qui so che aplha+beta < phi_aplha+phi_beta quindi sono nel caso sicuro
if not (IsEmpty(To[Position(P, P[i])][Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]+K[i]) > p+Position(P,P[j]+P[i])) then
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, K[j]+K[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, P[i]+P[j]) );
else
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, K[j]+K[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, P[i]+P[j]) );
fi; 
else
	if (p+Position(P,K[j]+K[i]) = p+Position(P,P[j]+P[i])) then
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, P[i]+P[j]) );
else
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, K[j]+K[i]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, K[j]+K[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, P[i]+P[j]) );
else
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, P[i]+P[j]) );
Add(KK2, +1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, K[j]+K[i]) );
fi; 
	fi; 
fi; 
fi; fi;

#######
############
###############

if (P[j]-P[i] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][p+Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]-K[i]) > p+Position(P,P[j]-P[i])) then
if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, Position(P, K[j]-K[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, Position(P, P[j]-P[i]) );
else
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, Position(P, K[j]-K[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, Position(P, P[j]-P[i]) );
fi; 
else
	if (p+Position(P,K[j]-K[i]) = p+Position(P,P[j]-P[i])) then
if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, Position(P, P[j]-P[i]) );
else
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, Position(P, K[j]-K[i]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, Position(P, P[j]-P[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, Position(P, K[j]-K[i]) );
else
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, Position(P, P[j]-P[i]) );
Add(KK2, +1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, Position(P, K[j]-K[i]) );
fi; 
	fi; 
fi; 
fi; fi;

else

if (P[i]+P[j] in P ) then #qui so che aplha+beta < phi_aplha+phi_beta quindi sono nel caso sicuro
if not (IsEmpty(To[Position(P, P[i])][Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]+K[i]) > p+Position(P,P[j]+P[i])) then
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, K[i]+K[j]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, P[j]+P[i]) );
else
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, K[i]+K[j]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, P[j]+P[i]) );
fi; 
else
	if (p+Position(P,K[j]+K[i]) = p+Position(P,P[j]+P[i])) then
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, P[j]+P[i]) );
else
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, K[i]+K[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, P[j]+P[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, K[i]+K[j]) );
else
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, P[j]+P[i]) );
Add(KK2, +1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, K[i]+K[j]) );
fi; 
	fi; 
fi; 
fi; fi;

######
############
###############

if (P[j]-P[i] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][p+Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]-K[i]) > p+Position(P,P[j]-P[i])) then
if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, Position(P, K[j]-K[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, Position(P, P[j]-P[i]) );
else
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, Position(P, K[j]-K[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, Position(P, P[j]-P[i]) );
fi; 
else
	if (p+Position(P,K[j]-K[i]) = p+Position(P,P[j]-P[i])) then
if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, Position(P, P[j]-P[i]) );
else
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, Position(P, K[j]-K[i]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))

if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, Position(P, P[j]-P[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, Position(P, K[j]-K[i]) );
else
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, Position(P, P[j]-P[i]) );
Add(KK2, +1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, Position(P, K[j]-K[i]) );
fi; 
	fi; 
fi; 
fi; fi;
fi;
fi;

SetEntrySCTable( TT, i,  j , KK1);
SetEntrySCTable( TT, i,  Position(P,K[j]), KK2);

##########################
###########################################
#################################################

		else
			if ( p+Position(P, K[j]) = p+Position(P,P[j]) ) then #beta = phi_beta

KK1:=[];

if (t[i]=1) then
if (t[j]=1) then

if (P[i]+P[j] in P ) then #qui so che aplha+beta < phi_aplha + phi_beta quindi sono nel caso sicuro
if not (IsEmpty(To[Position(P, P[i])][Position(P,P[j])][2])=true) then
if (t[Position(P,P[j]+P[i])]=1) then
Add(KK1, +2*To[Position(P, P[i])][Position(P,P[j])][2][1]);
Add(KK1, Position(P, P[i]+P[j]) );
fi; fi; fi;

if (P[j]-P[i] in P ) then
if not (IsEmpty(To[Position(P, P[i])][p+Position(P,P[j])][2])=true) then
if (t[Position(P,P[j]-P[i])]=1) then
if (p+Position(P,K[j]-K[i]) > p+Position(P,P[j]-P[i])) then
Add(KK1, +2*To[Position(P, P[i])][p+Position(P,P[j])][2][1]);
Add(KK1, Position(P, P[j]-P[i]) );
else 
	if (p+Position(P,K[j]-K[i]) = p+Position(P, P[j]-P[i])) then
Add(KK1, +2*To[Position(P, P[i])][p+Position(P,P[j])][2][1]);
Add(KK1, Position(P, P[j]-P[i]) );
	else
Add(KK1, +2*To[Position(P, P[i])][p+Position(P,P[j])][2][1]);
Add(KK1, Position(P, K[j]-K[i]) );
	fi;
fi; fi; fi; fi;

SetEntrySCTable( TT, i,  j , KK1);

else

if (P[i]+P[j] in P ) then #qui so che aplha+beta < phi_aplha+phi_beta quindi sono nel caso sicuro
if not (IsEmpty(To[Position(P, P[i])][Position(P,P[j])][2])=true) then
if (t[Position(P,P[j]+P[i])]=-1) then
Add(KK1, +2*To[Position(P, P[i])][Position(P,P[j])][2][1]);
Add(KK1, Position(P, P[i]+P[j]) );
fi; fi; fi;

if (P[j]-P[i] in P ) then
if not (IsEmpty(To[Position(P, P[i])][p+Position(P,P[j])][2])=true) then
if (t[Position(P,P[j]-P[i])]=-1) then
if (p+Position(P,K[j]-K[i]) > p+Position(P,P[j]-P[i])) then
Add(KK1, +2*To[Position(P, P[i])][p+Position(P,P[j])][2][1]);
Add(KK1, Position(P, P[j]-P[i]) );
else 
	if (p+Position(P,K[j]-K[i]) = p+Position(P, P[j]-P[i])) then
Add(KK1, +2*To[Position(P, P[i])][p+Position(P,P[j])][2][1]);
Add(KK1, Position(P, P[j]-P[i]) );
	else
Add(KK1, +2*To[Position(P, P[i])][p+Position(P,P[j])][2][1]);
Add(KK1, Position(P, K[j]-K[i]) );
	fi;
fi; fi; fi; fi;
SetEntrySCTable( TT, i,  j , KK1);
fi;

else
if (t[j]=1) then

if (P[i]+P[j] in P ) then #qui so che aplha+beta < phi_aplha+phi_beta quindi sono nel caso sicuro
if not (IsEmpty(To[Position(P, P[i])][Position(P,P[j])][2])=true) then
if (t[Position(P,P[j]+P[i])]=-1) then
Add(KK1, +2*To[Position(P, P[i])][Position(P,P[j])][2][1]);
Add(KK1, Position(P, P[i]+P[j]) );
fi; fi; fi;

if (P[j]-P[i] in P ) then
if not (IsEmpty(To[Position(P, P[i])][p+Position(P,P[j])][2])=true) then
if (t[Position(P,P[j]-P[i])]=-1) then
if (p+Position(P,K[j]-K[i]) > p+Position(P,P[j]-P[i])) then
Add(KK1, +2*To[Position(P, P[i])][p+Position(P,P[j])][2][1]);
Add(KK1, Position(P, P[j]-P[i]) );
else 
	if (p+Position(P,K[j]-K[i]) = p+Position(P, P[j]-P[i])) then
Add(KK1, +2*To[Position(P, P[i])][p+Position(P,P[j])][2][1]);
Add(KK1, Position(P, P[j]-P[i]) );
	else
Add(KK1, +2*To[Position(P, P[i])][p+Position(P,P[j])][2][1]);
Add(KK1, Position(P, K[j]-K[i]) );
	fi;
fi; fi; fi; fi;

SetEntrySCTable( TT, i,  j , KK1);


else


if (P[i]+P[j] in P ) then #qui so che aplha+beta < phi_aplha+phi_beta quindi sono nel caso sicuro
if not (IsEmpty(To[Position(P, P[i])][Position(P,P[j])][2])=true) then
if (t[Position(P,P[j]+P[i])]=1) then
Add(KK1, -2*To[Position(P, P[i])][Position(P,P[j])][2][1]);
Add(KK1, Position(P, P[i]+P[j]) );
fi; fi; fi;

if (P[j]-P[i] in P ) then
if not (IsEmpty(To[Position(P, P[i])][p+Position(P,P[j])][2])=true) then
if (t[Position(P,P[j]-P[i])]=1) then
if (p+Position(P,K[j]-K[i]) > p+Position(P,P[j]-P[i])) then
Add(KK1, -2*To[Position(P, P[i])][p+Position(P,P[j])][2][1]);
Add(KK1, Position(P, P[j]-P[i]) );
else 
	if (p+Position(P,K[j]-K[i]) = p+Position(P, P[j]-P[i])) then
Add(KK1, -2*To[Position(P, P[i])][p+Position(P,P[j])][2][1]);
Add(KK1, Position(P, P[j]-P[i]) );
	else
Add(KK1, -2*To[Position(P, P[i])][p+Position(P,P[j])][2][1]);
Add(KK1, Position(P, K[j]-K[i]) );
	fi;
fi; fi; fi; fi;

SetEntrySCTable( TT, i,  j , KK1);

fi;
fi; fi; fi;
fi; fi;
od; od;

for i in [1..p] do #variazione alpha per gli X
		if ( Position(P, K[i]) > Position(P,P[i]) ) then #alpha< phi_alpha (ragiono con gli X)
KK1:=[];

if (t[i] = 1) then

if (P[i]+K[i] in P ) then
if not  (IsEmpty(To[i][Position(P,K[i])][2] ) = true ) then
if (t[Position(P, P[i]+K[i])]=-1) then
Add(KK1, -1*To[Position(P, P[i])][Position(P,K[i])][2][1]);
Add(KK1, Position(P, K[i]+P[i]) );
fi; fi; fi;

else

if (P[i]+K[i] in P ) then
if not  (IsEmpty(To[i][Position(P,K[i])][2] ) = true ) then
if (t[Position(P, P[i]+K[i])]=-1) then
Add(KK1, +1*To[Position(P, P[i])][Position(P,K[i])][2][1]);
Add(KK1, Position(P, K[i]+P[i]) );
fi; fi; fi;

fi;

SetEntrySCTable( TT, i,  Position(P,K[i]) , KK1);

		fi;
od;


#########################################################################################################
###########################################################################################################


#############################################################################################################################################
###################################################################################################################

#Prodotti [Y,Y] con indici differenti X>Y ovvero alpha > beta !!!!!!!


for i in [1..p-1] do #variazione alpha per gli X
	for j in [i+1..p] do #variazione beta < alpha per gli Y
		if ( Position(P, K[i]) > Position(P,P[i]) ) then #alpha< phi_alpha (ragiono con gli X)
			if ( p+Position(P, K[j]) > p+Position(P,P[j]) ) then # beta< phi_beta (ragiono con gli Y)

KK1:=[]; # i p+j
KK2:=[]; # i p+jbar
KK3:=[]; # ibar p+j
KK4:=[]; # ibar p+jbar

if (t[i]=1) then
if (t[j]=1) then

if (P[i]+P[j] in P ) then #qui so che aplha+beta < phi_aplha+phi_beta quindi sono nel caso sicuro
if not (IsEmpty(To[Position(P, P[i])][Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]+K[i]) > p+Position(P,P[j]+P[i])) then

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, Position(P, P[i]+P[j]) );
Add(KK2, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, Position(P, K[i]+K[j]) );
Add(KK3, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, Position(P, K[i]+K[j]) );
Add(KK4, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, Position(P, P[i]+P[j]) );
else
Add(KK1, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, Position(P, P[i]+P[j]) );
Add(KK2, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, Position(P, K[i]+K[j]) );
Add(KK3, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, Position(P, K[i]+K[j]) );
Add(KK4, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, Position(P, P[i]+P[j]) );
fi; 

else
	if (p+Position(P,K[j]+K[i]) = p+Position(P,P[j]+P[i])) then

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, Position(P, P[i]+P[j]) );
Add(KK4, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, Position(P, P[i]+P[j]) );
else
Add(KK2, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, Position(P, K[i]+K[j]) );
Add(KK3, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, Position(P, K[i]+K[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, Position(P, K[j]+K[i]) );
Add(KK2, -1*(-1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, Position(P, P[i]+P[j]) );
Add(KK3, -1*(-1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, Position(P, P[i]+P[j]) );
Add(KK4, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, Position(P,K[j]+K[i]) );
else
Add(KK1, -1*(-1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, Position(P, K[j]+K[i]) );
Add(KK2, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, Position(P, P[i]+P[j]) );
Add(KK3, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, Position(P, P[i]+P[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, Position(P, K[j]+K[i]) );
fi; 
	fi; 
fi; 
fi; fi;

############

if (P[j]-P[i] in P ) then 
if not (IsEmpty(To[p+Position(P, P[i])][Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]-K[i]) > p+Position(P,P[j]-P[i])) then

if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, Position(P, P[j]-P[i]) );
Add(KK2, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, Position(P, K[j]-K[i]) );
Add(KK3, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, Position(P, K[j]-K[i]) );
Add(KK4, (+1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, Position(P, P[j]-P[i]) );
else
Add(KK1, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, Position(P, P[j]-P[i]) );
Add(KK2, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, Position(P, K[j]-K[i]) );
Add(KK3, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, Position(P, K[j]-K[i]) );
Add(KK4, (+1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, Position(P, P[j]-P[i]) );
fi; 

else
	if (p+Position(P,K[j]-K[i]) = p+Position(P,P[j]-P[i])) then

if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, Position(P, P[j]-P[i]) );
Add(KK4, (+1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, Position(P, P[j]-P[i]) );
else
Add(KK2, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, Position(P, K[j]-K[i]) );
Add(KK3, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, Position(P, K[j]-K[i]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))

if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, Position(P, K[j]-K[i]) );
Add(KK2, -1*(-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, Position(P, P[j]-P[i]) );
Add(KK3, -1*(-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, Position(P, P[j]-P[i]) );
Add(KK4, (+1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, Position(P, K[j]-K[i]) );
else
Add(KK1, -1*(-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, Position(P, K[j]-K[i]) );
Add(KK2, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, Position(P, P[j]-P[i]) );
Add(KK3, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, Position(P, P[j]-P[i]) );
Add(KK4, -1*(+1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, Position(P, K[j]-K[i]) );
fi; 
	fi; 
fi; 
fi; fi;


##########

if (P[i]+K[j] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][Position(P,K[j])][2])=true) then
if (p+Position(P,K[i]+P[j]) > p+Position(P,P[i]+K[j])) then

if (t[Position(P,P[i]+K[j] )] = 1) then
Add(KK1, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, Position(P, P[i]+K[j]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, Position(P, K[i]+P[j]) );
Add(KK3, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, Position(P, K[i]+P[j]) );
Add(KK4, +1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, Position(P, P[i]+K[j]) );
else
Add(KK1, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, Position(P, P[i]+K[j]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, Position(P, K[i]+P[j]) );
Add(KK3, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, Position(P, K[i]+P[j]) );
Add(KK4, +1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, Position(P, P[i]+K[j]) );
fi; 
else
	if (p+Position(P,K[i]+P[j]) = p+Position(P,P[i]+K[j])) then

if (t[Position(P,K[j]+P[i] )] = 1) then
Add(KK1, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, Position(P, P[i]+K[j]) );
Add(KK4, +1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, Position(P, P[i]+K[j]) );
else
Add(KK2, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, Position(P, K[i]+P[j]) );
Add(KK3, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, Position(P, K[i]+P[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))

if (t[Position(P,P[i]+K[j] )] = 1) then
Add(KK1, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, Position(P, K[i]+P[j]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, Position(P, P[i]+K[j]) );
Add(KK3, -1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, Position(P, P[i]+K[j]) );
Add(KK4, +1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, Position(P, K[i]+P[j]) );
else
Add(KK1, -1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, Position(P, K[i]+P[j]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, Position(P, P[i]+K[j]) );
Add(KK3, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, Position(P, P[i]+K[j]) );
Add(KK4, -1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, Position(P, K[i]+P[j]) );
fi; 
	fi; 
fi; 
fi; fi;

##########

if (K[j]-P[i] in P ) then 
if not (IsEmpty(To[p+Position(P, P[i])][Position(P,K[j])][2])=true) then
if (p+Position(P,P[j]-K[i]) > p+Position(P,K[j]-P[i])) then

if (t[Position(P,K[j]-P[i] )] = 1) then
Add(KK1, (-1*To[p+Position(P, P[i])][Position(P,K[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, Position(P, K[j]-P[i]) );
Add(KK2, (+1*To[p+Position(P, P[i])][Position(P,K[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, Position(P, P[j]-K[i]) );
Add(KK3, (-1*To[p+Position(P, P[i])][Position(P,K[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, Position(P, P[j]-K[i]) );
Add(KK4, +1*(-1*To[p+Position(P, P[i])][Position(P,K[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, Position(P, K[j]-P[i]) );
else
Add(KK1, (-1*To[p+Position(P, P[i])][Position(P,K[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, Position(P, K[j]-P[i]) );
Add(KK2, (+1*To[p+Position(P, P[i])][Position(P,K[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, Position(P, P[j]-K[i]) );
Add(KK3, (-1*To[p+Position(P, P[i])][Position(P,K[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, Position(P, P[j]-K[i]) );
Add(KK4, +1*(-1*To[p+Position(P, P[i])][Position(P,K[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, Position(P, K[j]-P[i]) );
fi; 
else
	if (p+Position(P,P[j]-K[i]) = p+Position(P,K[j]-P[i])) then

if (t[Position(P,K[j]-P[i] )] = 1) then
Add(KK1, (-1*To[p+Position(P, P[i])][Position(P,K[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, Position(P, K[j]-P[i]) );
Add(KK4, +1*(-1*To[p+Position(P, P[i])][Position(P,K[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, Position(P, K[j]-P[i]) );
else
Add(KK2, (+1*To[p+Position(P, P[i])][Position(P,K[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, Position(P, P[j]-K[i]) );
Add(KK3, (-1*To[p+Position(P, P[i])][Position(P,K[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, Position(P, P[j]-K[i]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))

if (t[Position(P,K[j]-P[i] )] = 1) then
Add(KK1, (-1*To[p+Position(P, P[i])][Position(P,K[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, Position(P, P[j]-K[i]) );
Add(KK2, -1*(+1*To[p+Position(P, P[i])][Position(P,K[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, Position(P, K[j]-P[i]) );
Add(KK3, -1*(-1*To[p+Position(P, P[i])][Position(P,K[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, Position(P, K[j]-P[i]) );
Add(KK4, +1*(-1*To[p+Position(P, P[i])][Position(P,K[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, Position(P, P[j]-K[i]) );
else
Add(KK1, -1*(-1*To[p+Position(P, P[i])][Position(P,K[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, Position(P, P[j]-K[i]) );
Add(KK2, (+1*To[p+Position(P, P[i])][Position(P,K[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, Position(P, K[j]-P[i]) );
Add(KK3, (-1*To[p+Position(P, P[i])][Position(P,K[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, Position(P, K[j]-P[i]) );
Add(KK4, -1*(-1*To[p+Position(P, P[i])][Position(P,K[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, Position(P, P[j]-K[i]) );
fi; 
	fi; 
fi; 
fi; fi;

####################

else # ovvero t[j]=-1

#########################

if (P[i]+P[j] in P ) then #qui so che aplha+beta < phi_aplha+phi_beta quindi sono nel caso sicuro
if not (IsEmpty(To[Position(P, P[i])][Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]+K[i]) > p+Position(P,P[j]+P[i])) then

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, Position(P, P[i]+P[j]) );
Add(KK2, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, Position(P, K[i]+K[j]) );
Add(KK3, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, Position(P, K[i]+K[j]) );
Add(KK4, +1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, Position(P, P[i]+P[j]) );
else
Add(KK1, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, Position(P, P[i]+P[j]) );
Add(KK2, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, Position(P, K[i]+K[j]) );
Add(KK3, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, Position(P, K[i]+K[j]) );
Add(KK4, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, Position(P, P[i]+P[j]) );
fi; 
else
	if (p+Position(P,K[j]+K[i]) = p+Position(P,P[j]+P[i])) then

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, Position(P, P[i]+P[j]) );
Add(KK4, +1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, Position(P, P[i]+P[j]) );
else
Add(KK2, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, Position(P, K[i]+K[j]) );
Add(KK3, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, Position(P, K[i]+K[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, Position(P, K[j]+K[i]) );
Add(KK2, -1*(-1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, Position(P, P[i]+P[j]) );
Add(KK3, -1*(-1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, Position(P, P[i]+P[j]) );
Add(KK4, +1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, Position(P,K[j]+K[i]) );
else
Add(KK1, -1*(-1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, Position(P, K[j]+K[i]) );
Add(KK2, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, Position(P, P[i]+P[j]) );
Add(KK3, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, Position(P, P[i]+P[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, Position(P, K[j]+K[i]) );
fi; 
	fi; 
fi; 
fi; fi;

###############

if (P[j]-P[i] in P ) then 
if not (IsEmpty(To[p+Position(P, P[i])][Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]-K[i]) > p+Position(P,P[j]-P[i])) then

if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, Position(P, P[j]-P[i]) );
Add(KK2, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, Position(P, K[j]-K[i]) );
Add(KK3, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, Position(P, K[j]-K[i]) );
Add(KK4, +1*(+1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, Position(P, P[j]-P[i]) );
else
Add(KK1, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, Position(P, P[j]-P[i]) );
Add(KK2, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, Position(P, K[j]-K[i]) );
Add(KK3, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, Position(P, K[j]-K[i]) );
Add(KK4, +1*(+1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, Position(P, P[j]-P[i]) );
fi; 
else
	if (p+Position(P,K[j]-K[i]) = p+Position(P,P[j]-P[i])) then

if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, Position(P, P[j]-P[i]) );
Add(KK4, +1*(+1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, Position(P, P[j]-P[i]) );
else
Add(KK2, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, Position(P, K[j]-K[i]) );
Add(KK3, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, Position(P, K[j]-K[i]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))

if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, Position(P, K[j]-K[i]) );
Add(KK2, -1*(-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, Position(P, P[j]-P[i]) );
Add(KK3, -1*(-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, Position(P, P[j]-P[i]) );
Add(KK4, +1*(+1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, Position(P, K[j]-K[i]) );
else
Add(KK1, -1*(-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, Position(P, K[j]-K[i]) );
Add(KK2, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, Position(P, P[j]-P[i]) );
Add(KK3, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, Position(P, P[j]-P[i]) );
Add(KK4, -1*(+1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, Position(P, K[j]-K[i]) );
fi; 
	fi; 
fi; 
fi; fi;

##########

if (P[i]+K[j] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][Position(P,K[j])][2])=true) then
if (p+Position(P,K[i]+P[j]) > p+Position(P,P[i]+K[j])) then

if (t[Position(P,P[i]+K[j] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, Position(P, P[i]+K[j]) );
Add(KK2, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, Position(P, K[i]+P[j]) );
Add(KK3, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, Position(P, K[i]+P[j]) );
Add(KK4, +1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, Position(P, P[i]+K[j]) );
else
Add(KK1, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, Position(P, P[i]+K[j]) );
Add(KK2, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, Position(P, K[i]+P[j]) );
Add(KK3, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, Position(P, K[i]+P[j]) );
Add(KK4, +1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, Position(P, P[i]+K[j]) );
fi; 
else
	if (p+Position(P,K[i]+P[j]) = p+Position(P,P[i]+K[j])) then

if (t[Position(P,K[j]+P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, Position(P, P[i]+K[j]) );
Add(KK4, +1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, Position(P, P[i]+K[j]) );
else
Add(KK2, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, Position(P, K[i]+P[j]) );
Add(KK3, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, Position(P, K[i]+P[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))

if (t[Position(P,P[i]+K[j] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, Position(P, K[i]+P[j]) );
Add(KK2, -1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, Position(P, P[i]+K[j]) );
Add(KK3, -1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, Position(P, P[i]+K[j]) );
Add(KK4, +1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, Position(P, K[i]+P[j]) );
else
Add(KK1, -1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, Position(P, K[i]+P[j]) );
Add(KK2, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, Position(P, P[i]+K[j]) );
Add(KK3, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, Position(P, P[i]+K[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, Position(P, K[i]+P[j]) );
fi; 
	fi; 
fi; 
fi; fi;

#############

if (K[j]-P[i] in P ) then 
if not (IsEmpty(To[p+Position(P, P[i])][Position(P,K[j])][2])=true) then
if (p+Position(P,P[j]-K[i]) > p+Position(P,K[j]-P[i])) then

if (t[Position(P,K[j]-P[i] )] = 1) then
Add(KK1, (+1*To[p+Position(P, P[i])][Position(P,K[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, Position(P, K[j]-P[i]) );
Add(KK2, (-1*To[p+Position(P, P[i])][Position(P,K[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, Position(P, P[j]-K[i]) );
Add(KK3, (+1*To[p+Position(P, P[i])][Position(P,K[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, Position(P, P[j]-K[i]) );
Add(KK4, +1*(+1*To[p+Position(P, P[i])][Position(P,K[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, Position(P, K[j]-P[i]) );
else
Add(KK1, (+1*To[p+Position(P, P[i])][Position(P,K[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, Position(P, K[j]-P[i]) );
Add(KK2, (-1*To[p+Position(P, P[i])][Position(P,K[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, Position(P, P[j]-K[i]) );
Add(KK3, (+1*To[p+Position(P, P[i])][Position(P,K[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, Position(P, P[j]-K[i]) );
Add(KK4, +1*(+1*To[p+Position(P, P[i])][Position(P,K[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, Position(P, K[j]-P[i]) );
fi; 
else
	if (p+Position(P,P[j]-K[i]) = p+Position(P,K[j]-P[i])) then

if (t[Position(P,K[j]-P[i] )] = 1) then
Add(KK1, (+1*To[p+Position(P, P[i])][Position(P,K[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, Position(P, K[j]-P[i]) );
Add(KK4, +1*(+1*To[p+Position(P, P[i])][Position(P,K[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, Position(P, K[j]-P[i]) );
else
Add(KK2, (-1*To[p+Position(P, P[i])][Position(P,K[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, Position(P, P[j]-K[i]) );
Add(KK3, (+1*To[p+Position(P, P[i])][Position(P,K[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, Position(P, P[j]-K[i]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))

if (t[Position(P,K[j]-P[i] )] = 1) then
Add(KK1, (+1*To[p+Position(P, P[i])][Position(P,K[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, Position(P, P[j]-K[i]) );
Add(KK2, -1*(-1*To[p+Position(P, P[i])][Position(P,K[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, Position(P, K[j]-P[i]) );
Add(KK3, -1*(+1*To[p+Position(P, P[i])][Position(P,K[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, Position(P, K[j]-P[i]) );
Add(KK4, +1*(+1*To[p+Position(P, P[i])][Position(P,K[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, Position(P, P[j]-K[i]) );
else
Add(KK1, -1*(+1*To[p+Position(P, P[i])][Position(P,K[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, Position(P, P[j]-K[i]) );
Add(KK2, (-1*To[p+Position(P, P[i])][Position(P,K[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, Position(P, K[j]-P[i]) );
Add(KK3, (+1*To[p+Position(P, P[i])][Position(P,K[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, Position(P, K[j]-P[i]) );
Add(KK4, -1*(+1*To[p+Position(P, P[i])][Position(P,K[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, Position(P, P[j]-K[i]) );
fi; 
	fi; 
fi; 
fi; fi;

fi;

########################à

else #t[i]=-1
if (t[j]=1) then

#########################

if (P[i]+P[j] in P ) then #qui so che aplha+beta < phi_aplha+phi_beta quindi sono nel caso sicuro
if not (IsEmpty(To[Position(P, P[i])][Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]+K[i]) > p+Position(P,P[j]+P[i])) then

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, Position(P, P[i]+P[j]) );
Add(KK2, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, Position(P, K[i]+K[j]) );
Add(KK3, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, Position(P, K[i]+K[j]) );
Add(KK4, +1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, Position(P, P[i]+P[j]) );
else
Add(KK1, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, Position(P, P[i]+P[j]) );
Add(KK2, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, Position(P, K[i]+K[j]) );
Add(KK3, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, Position(P, K[i]+K[j]) );
Add(KK4, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, Position(P, P[i]+P[j]) );
fi; 
else
	if (p+Position(P,K[j]+K[i]) = p+Position(P,P[j]+P[i])) then

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, Position(P, P[i]+P[j]) );
Add(KK4, +1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, Position(P, P[i]+P[j]) );
else
Add(KK2, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, Position(P, K[i]+K[j]) );
Add(KK3, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, Position(P, K[i]+K[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, Position(P, K[j]+K[i]) );
Add(KK2, -1*(-1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, Position(P, P[i]+P[j]) );
Add(KK3, -1*(-1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, Position(P, P[i]+P[j]) );
Add(KK4, +1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, Position(P,K[j]+K[i]) );
else
Add(KK1, -1*(-1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, Position(P, K[j]+K[i]) );
Add(KK2, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, Position(P, P[i]+P[j]) );
Add(KK3, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, Position(P, P[i]+P[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, Position(P, K[j]+K[i]) );
fi; 
	fi; 
fi; 
fi; fi;

###############

if (P[j]-P[i] in P ) then 
if not (IsEmpty(To[p+Position(P, P[i])][Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]-K[i]) > p+Position(P,P[j]-P[i])) then

if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, Position(P, P[j]-P[i]) );
Add(KK2, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, Position(P, K[j]-K[i]) );
Add(KK3, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, Position(P, K[j]-K[i]) );
Add(KK4, +1*(+1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, Position(P, P[j]-P[i]) );
else
Add(KK1, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, Position(P, P[j]-P[i]) );
Add(KK2, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, Position(P, K[j]-K[i]) );
Add(KK3, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, Position(P, K[j]-K[i]) );
Add(KK4, +1*(+1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, Position(P, P[j]-P[i]) );
fi; 

else
	if (p+Position(P,K[j]-K[i]) = p+Position(P,P[j]-P[i])) then

if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, Position(P, P[j]-P[i]) );
Add(KK4, +1*(+1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, Position(P, P[j]-P[i]) );
else
Add(KK2, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, Position(P, K[j]-K[i]) );
Add(KK3, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, Position(P, K[j]-K[i]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))

if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, Position(P, K[j]-K[i]) );
Add(KK2, -1*(-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, Position(P, P[j]-P[i]) );
Add(KK3, -1*(-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, Position(P, P[j]-P[i]) );
Add(KK4, +1*(+1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, Position(P, K[j]-K[i]) );
else
Add(KK1, -1*(-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, Position(P, K[j]-K[i]) );
Add(KK2, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, Position(P, P[j]-P[i]) );
Add(KK3, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, Position(P, P[j]-P[i]) );
Add(KK4, -1*(+1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, Position(P, K[j]-K[i]) );
fi; 
	fi; 
fi; 
fi; fi;

#############

if (P[i]+K[j] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][Position(P,K[j])][2])=true) then
if (p+Position(P,K[i]+P[j]) > p+Position(P,P[i]+K[j])) then

if (t[Position(P,P[i]+K[j] )] = 1) then
Add(KK1, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, Position(P, P[i]+K[j]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, Position(P, K[i]+P[j]) );
Add(KK3, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, Position(P, K[i]+P[j]) );
Add(KK4, +1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, Position(P, P[i]+K[j]) );
else
Add(KK1, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, Position(P, P[i]+K[j]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, Position(P, K[i]+P[j]) );
Add(KK3, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, Position(P, K[i]+P[j]) );
Add(KK4, +1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, Position(P, P[i]+K[j]) );
fi; 
else
	if (p+Position(P,K[i]+P[j]) = p+Position(P,P[i]+K[j])) then

if (t[Position(P,K[j]+P[i] )] = 1) then
Add(KK1, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, Position(P, P[i]+K[j]) );
Add(KK4, +1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, Position(P, P[i]+K[j]) );
else
Add(KK2, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, Position(P, K[i]+P[j]) );
Add(KK3, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, Position(P, K[i]+P[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))

if (t[Position(P,P[i]+K[j] )] = 1) then
Add(KK1, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, Position(P, K[i]+P[j]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, Position(P, P[i]+K[j]) );
Add(KK3, -1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, Position(P, P[i]+K[j]) );
Add(KK4, +1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, Position(P, K[i]+P[j]) );
else
Add(KK1, -1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, Position(P, K[i]+P[j]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, Position(P, P[i]+K[j]) );
Add(KK3, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, Position(P, P[i]+K[j]) );
Add(KK4, -1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, Position(P, K[i]+P[j]) );
fi; 
	fi; 
fi; 
fi; fi;

#############

if (K[j]-P[i] in P ) then 
if not (IsEmpty(To[p+Position(P, P[i])][Position(P,K[j])][2])=true) then
if (p+Position(P,P[j]-K[i]) > p+Position(P,K[j]-P[i])) then

if (t[Position(P,K[j]-P[i] )] = 1) then
Add(KK1, (-1*To[p+Position(P, P[i])][Position(P,K[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, Position(P, K[j]-P[i]) );
Add(KK2, (+1*To[p+Position(P, P[i])][Position(P,K[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, Position(P, P[j]-K[i]) );
Add(KK3, (-1*To[p+Position(P, P[i])][Position(P,K[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, Position(P, P[j]-K[i]) );
Add(KK4, +1*(-1*To[p+Position(P, P[i])][Position(P,K[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, Position(P, K[j]-P[i]) );
else
Add(KK1, (-1*To[p+Position(P, P[i])][Position(P,K[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, Position(P, K[j]-P[i]) );
Add(KK2, (+1*To[p+Position(P, P[i])][Position(P,K[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, Position(P, P[j]-K[i]) );
Add(KK3, (-1*To[p+Position(P, P[i])][Position(P,K[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, Position(P, P[j]-K[i]) );
Add(KK4, +1*(-1*To[p+Position(P, P[i])][Position(P,K[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, Position(P, K[j]-P[i]) );
fi; 
else
	if (p+Position(P,P[j]-K[i]) = p+Position(P,K[j]-P[i])) then

if (t[Position(P,K[j]-P[i] )] = 1) then
Add(KK1, (-1*To[p+Position(P, P[i])][Position(P,K[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, Position(P, K[j]-P[i]) );
Add(KK4, +1*(-1*To[p+Position(P, P[i])][Position(P,K[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, Position(P, K[j]-P[i]) );
else
Add(KK2, (+1*To[p+Position(P, P[i])][Position(P,K[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, Position(P, P[j]-K[i]) );
Add(KK3, (-1*To[p+Position(P, P[i])][Position(P,K[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, Position(P, P[j]-K[i]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))

if (t[Position(P,K[j]-P[i] )] = 1) then
Add(KK1, (-1*To[p+Position(P, P[i])][Position(P,K[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, Position(P, P[j]-K[i]) );
Add(KK2, -1*(+1*To[p+Position(P, P[i])][Position(P,K[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, Position(P, K[j]-P[i]) );
Add(KK3, -1*(-1*To[p+Position(P, P[i])][Position(P,K[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, Position(P, K[j]-P[i]) );
Add(KK4, +1*(-1*To[p+Position(P, P[i])][Position(P,K[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, Position(P, P[j]-K[i]) );
else
Add(KK1, -1*(-1*To[p+Position(P, P[i])][Position(P,K[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, Position(P, P[j]-K[i]) );
Add(KK2, (+1*To[p+Position(P, P[i])][Position(P,K[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, Position(P, K[j]-P[i]) );
Add(KK3, (-1*To[p+Position(P, P[i])][Position(P,K[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, Position(P, K[j]-P[i]) );
Add(KK4, -1*(-1*To[p+Position(P, P[i])][Position(P,K[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, Position(P, P[j]-K[i]) );
fi; 
	fi; 
fi; 
fi; fi;

else #t[j]=-1

#########################

if (P[i]+P[j] in P ) then #qui so che aplha+beta < phi_aplha+phi_beta quindi sono nel caso sicuro
if not (IsEmpty(To[Position(P, P[i])][Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]+K[i]) > p+Position(P,P[j]+P[i])) then

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, Position(P, P[i]+P[j]) );
Add(KK2, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, Position(P, K[i]+K[j]) );
Add(KK3, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, Position(P, K[i]+K[j]) );
Add(KK4, +1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, Position(P, P[i]+P[j]) );
else
Add(KK1, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, Position(P, P[i]+P[j]) );
Add(KK2, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, Position(P, K[i]+K[j]) );
Add(KK3, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, Position(P, K[i]+K[j]) );
Add(KK4, +1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, Position(P, P[i]+P[j]) );
fi; 
else
	if (p+Position(P,K[j]+K[i]) = p+Position(P,P[j]+P[i])) then

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, Position(P, P[i]+P[j]) );
Add(KK4, +1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, Position(P, P[i]+P[j]) );
else
Add(KK2, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, Position(P, K[i]+K[j]) );
Add(KK3, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, Position(P, K[i]+K[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, Position(P, K[j]+K[i]) );
Add(KK2, -1*(-1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, Position(P, P[i]+P[j]) );
Add(KK3, -1*(-1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, Position(P, P[i]+P[j]) );
Add(KK4, +1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, Position(P,K[j]+K[i]) );
else
Add(KK1, -1*(-1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, Position(P, K[j]+K[i]) );
Add(KK2, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, Position(P, P[i]+P[j]) );
Add(KK3, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, Position(P, P[i]+P[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, Position(P, K[j]+K[i]) );
fi; 
	fi; 
fi; 
fi; fi;

###############

if (P[j]-P[i] in P ) then 
if not (IsEmpty(To[p+Position(P, P[i])][Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]-K[i]) > p+Position(P,P[j]-P[i])) then

if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, Position(P, P[j]-P[i]) );
Add(KK2, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, Position(P, K[j]-K[i]) );
Add(KK3, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, Position(P, K[j]-K[i]) );
Add(KK4, +1*(+1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, Position(P, P[j]-P[i]) );
else
Add(KK1, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, Position(P, P[j]-P[i]) );
Add(KK2, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, Position(P, K[j]-K[i]) );
Add(KK3, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, Position(P, K[j]-K[i]) );
Add(KK4, +1*(+1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, Position(P, P[j]-P[i]) );
fi; 
else
	if (p+Position(P,K[j]-K[i]) = p+Position(P,P[j]-P[i])) then

if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, Position(P, P[j]-P[i]) );
Add(KK4, +1*(+1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, Position(P, P[j]-P[i]) );
else
Add(KK2, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, Position(P, K[j]-K[i]) );
Add(KK3, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, Position(P, K[j]-K[i]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))

if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, Position(P, K[j]-K[i]) );
Add(KK2, -1*(-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, Position(P, P[j]-P[i]) );
Add(KK3, -1*(-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, Position(P, P[j]-P[i]) );
Add(KK4, +1*(+1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, Position(P, K[j]-K[i]) );
else
Add(KK1, -1*(-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, Position(P, K[j]-K[i]) );
Add(KK2, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, Position(P, P[j]-P[i]) );
Add(KK3, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, Position(P, P[j]-P[i]) );
Add(KK4, -1*(+1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, Position(P, K[j]-K[i]) );
fi; 
	fi; 
fi; 
fi; fi;

#############

if (P[i]+K[j] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][Position(P,K[j])][2])=true) then
if (p+Position(P,K[i]+P[j]) > p+Position(P,P[i]+K[j])) then

if (t[Position(P,P[i]+K[j] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, Position(P, P[i]+K[j]) );
Add(KK2, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, Position(P, K[i]+P[j]) );
Add(KK3, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, Position(P, K[i]+P[j]) );
Add(KK4, +1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, Position(P, P[i]+K[j]) );
else
Add(KK1, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, Position(P, P[i]+K[j]) );
Add(KK2, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, Position(P, K[i]+P[j]) );
Add(KK3, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, Position(P, K[i]+P[j]) );
Add(KK4, +1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, Position(P, P[i]+K[j]) );
fi; 

else
	if (p+Position(P,K[i]+P[j]) = p+Position(P,P[i]+K[j])) then

if (t[Position(P,K[j]+P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, Position(P, P[i]+K[j]) );
Add(KK4, +1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, Position(P, P[i]+K[j]) );
else
Add(KK2, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, Position(P, K[i]+P[j]) );
Add(KK3, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, Position(P, K[i]+P[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))

if (t[Position(P,P[i]+K[j] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, Position(P, K[i]+P[j]) );
Add(KK2, -1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, Position(P, P[i]+K[j]) );
Add(KK3, -1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, Position(P, P[i]+K[j]) );
Add(KK4, +1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, Position(P, K[i]+P[j]) );
else
Add(KK1, -1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, Position(P, K[i]+P[j]) );
Add(KK2, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, Position(P, P[i]+K[j]) );
Add(KK3, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, Position(P, P[i]+K[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, Position(P, K[i]+P[j]) );
fi; 
	fi; 
fi; 
fi; fi;

#############

if (K[j]-P[i] in P ) then 
if not (IsEmpty(To[p+Position(P, P[i])][Position(P,K[j])][2])=true) then
if (p+Position(P,P[j]-K[i]) > p+Position(P,K[j]-P[i])) then

if (t[Position(P,K[j]-P[i] )] = 1) then
Add(KK1, (+1*To[p+Position(P, P[i])][Position(P,K[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, Position(P, K[j]-P[i]) );
Add(KK2, (-1*To[p+Position(P, P[i])][Position(P,K[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, Position(P, P[j]-K[i]) );
Add(KK3, (+1*To[p+Position(P, P[i])][Position(P,K[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, Position(P, P[j]-K[i]) );
Add(KK4, +1*(+1*To[p+Position(P, P[i])][Position(P,K[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, Position(P, K[j]-P[i]) );
else
Add(KK1, (+1*To[p+Position(P, P[i])][Position(P,K[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, Position(P, K[j]-P[i]) );
Add(KK2, (-1*To[p+Position(P, P[i])][Position(P,K[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, Position(P, P[j]-K[i]) );
Add(KK3, (+1*To[p+Position(P, P[i])][Position(P,K[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, Position(P, P[j]-K[i]) );
Add(KK4, +1*(+1*To[p+Position(P, P[i])][Position(P,K[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, Position(P, K[j]-P[i]) );
fi; 

else
	if (p+Position(P,P[j]-K[i]) = p+Position(P,K[j]-P[i])) then

if (t[Position(P,K[j]-P[i] )] = 1) then
Add(KK1, (+1*To[p+Position(P, P[i])][Position(P,K[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, Position(P, K[j]-P[i]) );
Add(KK4, +1*(+1*To[p+Position(P, P[i])][Position(P,K[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, Position(P, K[j]-P[i]) );
else
Add(KK2, (-1*To[p+Position(P, P[i])][Position(P,K[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, Position(P, P[j]-K[i]) );
Add(KK3, (+1*To[p+Position(P, P[i])][Position(P,K[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, Position(P, P[j]-K[i]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))

if (t[Position(P,K[j]-P[i] )] = 1) then
Add(KK1, (+1*To[p+Position(P, P[i])][Position(P,K[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, Position(P, P[j]-K[i]) );
Add(KK2, -1*(-1*To[p+Position(P, P[i])][Position(P,K[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, Position(P, K[j]-P[i]) );
Add(KK3, -1*(+1*To[p+Position(P, P[i])][Position(P,K[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, Position(P, K[j]-P[i]) );
Add(KK4, +1*(+1*To[p+Position(P, P[i])][Position(P,K[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, Position(P, P[j]-K[i]) );
else
Add(KK1, -1*(+1*To[p+Position(P, P[i])][Position(P,K[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, Position(P, P[j]-K[i]) );
Add(KK2, (-1*To[p+Position(P, P[i])][Position(P,K[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, Position(P, K[j]-P[i]) );
Add(KK3, (+1*To[p+Position(P, P[i])][Position(P,K[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, Position(P, K[j]-P[i]) );
Add(KK4, -1*(+1*To[p+Position(P, P[i])][Position(P,K[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, Position(P, P[j]-K[i]) );
fi; 
	fi; 
fi; 
fi; fi;

fi; fi;

SetEntrySCTable( TT, p+i,  p+j , KK1);
SetEntrySCTable( TT, p+i,  p+Position(P, K[j]), KK2);
SetEntrySCTable( TT, p+Position(P,K[i]), p+j, KK3);
SetEntrySCTable( TT, p+Position(P,K[i]) ,p+Position(P, K[j]), KK4);

else

	if ( p+Position(P, K[j]) = p+Position(P,P[j]) ) then #beta = phi_beta

KK1:=[];
KK2:=[];

if (t[i]=1) then
if (t[j]=1) then

if (P[i]+P[j] in P ) then #qui so che aplha+beta < phi_aplha+phi_beta quindi sono nel caso sicuro
if not (IsEmpty(To[Position(P, P[i])][Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]+K[i]) > p+Position(P,P[j]+P[i])) then

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, P[i]+P[j]) );
Add(KK2, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, K[i]+K[j]) );
else
Add(KK1, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, P[i]+P[j]) );
Add(KK2, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, K[i]+K[j]) );
fi; 
else
	if (p+Position(P,K[j]+K[i]) = p+Position(P,P[j]+P[i])) then

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, P[i]+P[j]) );
else
Add(KK2, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, K[i]+K[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, K[j]+K[i]) );
Add(KK2, -1*(-1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, P[i]+P[j]) );
else
Add(KK1, -1*(-1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, K[j]+K[i]) );
Add(KK2, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, P[i]+P[j]) );
fi; 
	fi; 
fi; 
fi; fi;

###############

if (P[j]-P[i] in P ) then 
if not (IsEmpty(To[p+Position(P, P[i])][Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]-K[i]) > p+Position(P,P[j]-P[i])) then

if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, P[j]-P[i]) );
Add(KK2, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, K[j]-K[i]) );
else
Add(KK1, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, P[j]-P[i]) );
Add(KK2, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, K[j]-K[i]) );
fi; 
else
	if (p+Position(P,K[j]-K[i]) = p+Position(P,P[j]-P[i])) then

if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, P[j]-P[i]) );
else
Add(KK2, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, K[j]-K[i]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))

if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, K[j]-K[i]) );
Add(KK2, -1*(-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, P[j]-P[i]) );
else
Add(KK1, -1*(-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, K[j]-K[i]) );
Add(KK2, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, P[j]-P[i]) );
fi; 
	fi; 
fi; 
fi; fi;

####################
else #t[j]=-1

if (P[i]+P[j] in P ) then #qui so che aplha+beta < phi_aplha+phi_beta quindi sono nel caso sicuro
if not (IsEmpty(To[Position(P, P[i])][Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]+K[i]) > p+Position(P,P[j]+P[i])) then

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, K[j]+K[i]) );
Add(KK2, +1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, P[i]+P[j]) );
else
Add(KK1, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, K[j]+K[i]) );
Add(KK2, +1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, P[i]+P[j]) );
fi; 
else
	if (p+Position(P,K[j]+K[i]) = p+Position(P,P[j]+P[i])) then

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK2, +1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, P[i]+P[j]) );
else
Add(KK1, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, K[j]+K[i]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, -1*(-1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, K[j]+K[i]) );
Add(KK2, +1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, P[i]+P[j]) );
else
Add(KK1, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, P[i]+P[j]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, K[j]+K[i]) );
fi; 
	fi; 
fi; 
fi; fi;

###############

if (P[j]-P[i] in P ) then 
if not (IsEmpty(To[p+Position(P, P[i])][Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]-K[i]) > p+Position(P,P[j]-P[i])) then
if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, K[j]-K[i]) );
Add(KK2, +1*(+1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, P[j]-P[i]) );
else
Add(KK1, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, K[j]-K[i]) );
Add(KK2, +1*(+1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, P[j]-P[i]) );
fi; 

else
	if (p+Position(P,K[j]-K[i]) = p+Position(P,P[j]-P[i])) then

if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK2, +1*(+1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, P[j]-P[i]) );
else
Add(KK1, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, K[j]-K[i]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))

if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, -1*(-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, P[j]-P[i]) );
Add(KK2, +1*(+1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, K[j]-K[i]) );
else
Add(KK1, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, P[j]-P[i]) );
Add(KK2, -1*(+1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, K[j]-K[i]) );
fi; 
	fi; 
fi; 
fi; fi;

fi;

else #t[i]=-1
if (t[j]=1) then

if (P[i]+P[j] in P ) then #qui so che aplha+beta < phi_aplha+phi_beta quindi sono nel caso sicuro
if not (IsEmpty(To[Position(P, P[i])][Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]+K[i]) > p+Position(P,P[j]+P[i])) then
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, P[i]+P[j]) );
Add(KK2, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, K[i]+K[j]) );
else
Add(KK1, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, P[i]+P[j]) );
Add(KK2, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, K[i]+K[j]) );
fi; 

else
	if (p+Position(P,K[j]+K[i]) = p+Position(P,P[j]+P[i])) then

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, P[i]+P[j]) );
else
Add(KK2, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, K[i]+K[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, K[j]+K[i]) );
Add(KK2, -1*(-1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, P[i]+P[j]) );
else
Add(KK1, -1*(-1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, K[j]+K[i]) );
Add(KK2, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, P[i]+P[j]) );
fi; 
	fi; 
fi; 
fi; fi;


#######
############
###############


if (P[j]-P[i] in P ) then 
if not (IsEmpty(To[p+Position(P, P[i])][Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]-K[i]) > p+Position(P,P[j]-P[i])) then

if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, P[j]-P[i]) );
Add(KK2, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, K[j]-K[i]) );
else
Add(KK1, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, P[j]-P[i]) );
Add(KK2, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, K[j]-K[i]) );
fi; 

else
	if (p+Position(P,K[j]-K[i]) = p+Position(P,P[j]-P[i])) then

if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, P[j]-P[i]) );
else

Add(KK2, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, K[j]-K[i]) );

fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))

if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, K[j]-K[i]) );
Add(KK2, -1*(-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, P[j]-P[i]) );
else
Add(KK1, -1*(-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, K[j]-K[i]) );
Add(KK2, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, P[j]-P[i]) );
fi; 
	fi; 

fi; 
fi; fi;


###########
##################
####################


else #t[i]=t[j]=-1


if (P[i]+P[j] in P ) then #qui so che aplha+beta < phi_aplha+phi_beta quindi sono nel caso sicuro
if not (IsEmpty(To[Position(P, P[i])][Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]+K[i]) > p+Position(P,P[j]+P[i])) then
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, K[i]+K[j]) );
Add(KK2, +1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, P[j]+P[i]) );
else
Add(KK1, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, K[i]+K[j]) );
Add(KK2, +1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, P[j]+P[i]) );
fi; 

else
	if (p+Position(P,K[j]+K[i]) = p+Position(P,P[j]+P[i])) then

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK2, +1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, P[j]+P[i]) );
else
Add(KK1, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, K[i]+K[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, -1*(-1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, P[j]+P[i]) );
Add(KK2, +1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, K[i]+K[j]) );
else
Add(KK1, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, P[j]+P[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, K[i]+K[j]) );
fi; 
	fi; 

fi; 
fi; fi;


#######
############
###############

if (P[j]-P[i] in P ) then 
if not (IsEmpty(To[p+Position(P, P[i])][Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]-K[i]) > p+Position(P,P[j]-P[i])) then
if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, K[j]-K[i]) );
Add(KK2, +1*(+1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, P[j]-P[i]) );
else
Add(KK1, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, K[j]-K[i]) );
Add(KK2, +1*(+1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, P[j]-P[i]) );
fi; 

else
	if (p+Position(P,K[j]-K[i]) = p+Position(P,P[j]-P[i])) then

if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK2, +1*(+1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, P[j]-P[i]) );
else
Add(KK1, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, K[j]-K[i]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))

if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, -1*(-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, P[j]-P[i]) );
Add(KK2, +1*(+1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, K[j]-K[i]) );
else
Add(KK1, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, P[j]-P[i]) );
Add(KK2, -1*(+1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, K[j]-K[i]) );
fi; 
	fi; 

fi; 
fi; fi;

fi; fi;

SetEntrySCTable( TT, p+i,  p+j , KK1);
SetEntrySCTable( TT, p+Position(P,K[i]), p+j, KK2);

fi;
fi;

else
	if ( p+Position(P, K[i]) = p+Position(P,P[i]) ) then #alpha< phi_alpha
			if ( p+Position(P, K[j]) > p+Position(P,P[j]) ) then #beta< phi_beta

KK1:=[];
KK2:=[];

if (t[i]=1) then
if (t[j]=1) then

if (P[i]+P[j] in P ) then #qui so che aplha+beta < phi_aplha+phi_beta quindi sono nel caso sicuro
if not (IsEmpty(To[Position(P, P[i])][Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]+K[i]) > p+Position(P,P[j]+P[i])) then

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, P[i]+P[j]) );
Add(KK2, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, K[i]+K[j]) );
else
Add(KK1, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, P[i]+P[j]) );
Add(KK2, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, K[i]+K[j]) );
fi; 

else
	if (p+Position(P,K[j]+K[i]) = p+Position(P,P[j]+P[i])) then

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, P[i]+P[j]) );
else
Add(KK2, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, K[i]+K[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, K[j]+K[i]) );
Add(KK2, -1*(-1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, P[i]+P[j]) );
else
Add(KK1, -1*(-1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, K[j]+K[i]) );
Add(KK2, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, P[i]+P[j]) );
fi; 
	fi; 

fi; 
fi; fi;


#######
############
###############


if (P[j]-P[i] in P ) then 
if not (IsEmpty(To[p+Position(P, P[i])][Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]-K[i]) > p+Position(P,P[j]-P[i])) then
if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, P[j]-P[i]) );
Add(KK2, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, K[j]-K[i]) );
else
Add(KK1, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, P[j]-P[i]) );
Add(KK2, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, K[j]-K[i]) );
fi; 

else
	if (p+Position(P,K[j]-K[i]) = p+Position(P,P[j]-P[i])) then

if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, P[j]-P[i]) );
else
Add(KK2, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, K[j]-K[i]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))

if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, K[j]-K[i]) );
Add(KK2, -1*(-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, P[j]-P[i]) );
else
Add(KK1, -1*(-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, K[j]-K[i]) );
Add(KK2, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, P[j]-P[i]) );
fi; 
	fi; 

fi; 
fi; fi;


###########
##################
####################

else #t[j]=-1

if (P[i]+P[j] in P ) then #qui so che aplha+beta < phi_aplha+phi_beta quindi sono nel caso sicuro
if not (IsEmpty(To[Position(P, P[i])][Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]+K[i]) > p+Position(P,P[j]+P[i])) then
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, P[i]+P[j]) );
Add(KK2, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, K[i]+K[j]) );
else
Add(KK1, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, P[i]+P[j]) );
Add(KK2, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, K[i]+K[j]) );
fi; 

else
	if (p+Position(P,K[j]+K[i]) = p+Position(P,P[j]+P[i])) then
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, P[i]+P[j]) );
else
Add(KK2, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, K[i]+K[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, K[j]+K[i]) );
Add(KK2, -1*(-1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, P[i]+P[j]) );
else
Add(KK1, -1*(-1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, K[j]+K[i]) );
Add(KK2, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, P[i]+P[j]) );
fi; 
	fi; 

fi; 
fi; fi;


#######
############
###############


if (P[j]-P[i] in P ) then 
if not (IsEmpty(To[p+Position(P, P[i])][Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]-K[i]) > p+Position(P,P[j]-P[i])) then
if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, P[j]-P[i]) );
Add(KK2, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, K[j]-K[i]) );
else
Add(KK1, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, P[j]-P[i]) );
Add(KK2, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, K[j]-K[i]) );
fi; 

else
	if (p+Position(P,K[j]-K[i]) = p+Position(P,P[j]-P[i])) then

if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, P[j]-P[i]) );
else
Add(KK2, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, K[j]-K[i]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))

if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, K[j]-K[i]) );
Add(KK2, -1*(-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, P[j]-P[i]) );
else
Add(KK1, -1*(-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, K[j]-K[i]) );
Add(KK2, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, P[j]-P[i]) );
fi; 
	fi; 

fi; 
fi; fi;


###########
##################
####################

fi;

else #t[i]=-1

if (t[j]=1) then

if (P[i]+P[j] in P ) then #qui so che aplha+beta < phi_aplha+phi_beta quindi sono nel caso sicuro
if not (IsEmpty(To[Position(P, P[i])][Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]+K[i]) > p+Position(P,P[j]+P[i])) then
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, K[j]+K[i]) );
Add(KK2, +1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, P[i]+P[j]) );
else
Add(KK1, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, K[j]+K[i]) );
Add(KK2, +1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, P[i]+P[j]) );
fi; 

else
	if (p+Position(P,K[j]+K[i]) = p+Position(P,P[j]+P[i])) then
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK2, +1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, P[i]+P[j]) );
else
Add(KK1, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, K[j]+K[i]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, -1*(-1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, K[j]+K[i]) );
Add(KK2, +1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, P[i]+P[j]) );
else

Add(KK1, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, P[i]+P[j]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, K[j]+K[i]) );
fi; 
	fi; 

fi; 
fi; fi;


#######
############
###############

if (P[j]-P[i] in P ) then 
if not (IsEmpty(To[p+Position(P, P[i])][Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]-K[i]) > p+Position(P,P[j]-P[i])) then

if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, K[j]-K[i]) );
Add(KK2, +1*(+1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, P[j]-P[i]) );
else
Add(KK1, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, K[j]-K[i]) );
Add(KK2, +1*(+1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, P[j]-P[i]) );
fi; 

else
	if (p+Position(P,K[j]-K[i]) = p+Position(P,P[j]-P[i])) then

if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK2, +1*(+1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, P[j]-P[i]) );
else
Add(KK1, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, K[j]-K[i]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))

if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, -1*(-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, P[j]-P[i]) );
Add(KK2, +1*(+1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, K[j]-K[i]) );
else
Add(KK1, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, P[j]-P[i]) );
Add(KK2, -1*(+1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, K[j]-K[i]) );
fi; 
	fi; 

fi; 
fi; fi;


else


if (P[i]+P[j] in P ) then #qui so che aplha+beta < phi_aplha+phi_beta quindi sono nel caso sicuro
if not (IsEmpty(To[Position(P, P[i])][Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]+K[i]) > p+Position(P,P[j]+P[i])) then
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, K[i]+K[j]) );
Add(KK2, +1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, P[j]+P[i]) );
else
Add(KK1, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, K[i]+K[j]) );
Add(KK2, +1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, P[j]+P[i]) );
fi; 

else
	if (p+Position(P,K[j]+K[i]) = p+Position(P,P[j]+P[i])) then

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK2, +1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, P[j]+P[i]) );
else
Add(KK1, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, K[i]+K[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, -1*(-1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, P[j]+P[i]) );
Add(KK2, +1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, K[i]+K[j]) );
else
Add(KK1, (-1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, P[j]+P[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, K[i]+K[j]) );
fi; 
	fi; 

fi; 
fi; fi;


#######
############
###############

if (P[j]-P[i] in P ) then 
if not (IsEmpty(To[p+Position(P, P[i])][Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]-K[i]) > p+Position(P,P[j]-P[i])) then
if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, K[j]-K[i]) );
Add(KK2, +1*(+1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, P[j]-P[i]) );
else
Add(KK1, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, K[j]-K[i]) );
Add(KK2, +1*(+1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, P[j]-P[i]) );
fi; 

else
	if (p+Position(P,K[j]-K[i]) = p+Position(P,P[j]-P[i])) then

if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK2, +1*(+1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, P[j]-P[i]) );
else
Add(KK1, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, K[j]-K[i]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))

if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, -1*(-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, P[j]-P[i]) );
Add(KK2, +1*(+1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, K[j]-K[i]) );
else
Add(KK1, (-1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]+1*To[p+Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, Position(P, P[j]-P[i]) );
Add(KK2, -1*(+1*To[p+Position(P, P[i])][Position(P,P[j])][2][1]-1*To[p+Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, Position(P, K[j]-K[i]) );
fi; 
	fi; 

fi; 
fi; fi;

fi;
fi;

SetEntrySCTable( TT, p+i,  p+j , KK1);
SetEntrySCTable( TT, p+i,  p+Position(P,K[j]), KK2);


##########################
###########################################
#################################################


		else
			if ( p+Position(P, K[j]) = p+Position(P,P[j]) ) then #beta = phi_beta


KK1:=[];

if (t[i]=1) then
if (t[j]=1) then

if (P[i]+P[j] in P ) then #qui so che aplha+beta < phi_aplha + phi_beta quindi sono nel caso sicuro
if not (IsEmpty(To[Position(P, P[i])][Position(P,P[j])][2])=true) then
if (t[Position(P,P[i]+P[j])]=1) then
Add(KK1, -2*To[Position(P, P[i])][Position(P,P[j])][2][1]);
Add(KK1, Position(P, P[i]+P[j]) );
fi; fi; fi;

if (P[j]-P[i] in P ) then
if not (IsEmpty(To[p+Position(P, P[i])][Position(P,P[j])][2])=true) then
if (t[Position(P,P[j]-P[i])]=1) then
if (p+Position(P,K[j]-K[i]) > p+Position(P,P[j]-P[i])) then
Add(KK1, -2*To[p+Position(P, P[i])][Position(P,P[j])][2][1]);
Add(KK1, Position(P, P[j]-P[i]) );
else 
	if (p+Position(P,K[j]-K[i]) = p+Position(P, P[j]-P[i])) then
Add(KK1, -2*To[p+Position(P, P[i])][Position(P,P[j])][2][1]);
Add(KK1, Position(P, P[j]-P[i]) );
	else
Add(KK1, -2*To[p+Position(P, P[i])][Position(P,P[j])][2][1]);
Add(KK1, Position(P, K[j]-K[i]) );
	fi;
fi; fi; fi; fi;

SetEntrySCTable( TT, p+i,  p+j , KK1);

else

if (P[i]+P[j] in P ) then #qui so che aplha+beta < phi_aplha+phi_beta quindi sono nel caso sicuro
if not (IsEmpty(To[Position(P, P[i])][Position(P,P[j])][2])=true) then
if (t[Position(P,P[j]+P[i])]=-1) then
Add(KK1, -2*To[Position(P, P[i])][Position(P,P[j])][2][1]);
Add(KK1, Position(P, P[i]+P[j]) );
fi; fi; fi;

if (P[j]-P[i] in P ) then
if not (IsEmpty(To[p+Position(P, P[i])][Position(P,P[j])][2])=true) then
if (t[Position(P,P[j]-P[i])]=-1) then
if (p+Position(P,K[j]-K[i]) > p+Position(P,P[j]-P[i])) then
Add(KK1, -2*To[p+Position(P, P[i])][Position(P,P[j])][2][1]);
Add(KK1, Position(P, P[j]-P[i]) );
else 
	if (p+Position(P,K[j]-K[i]) = p+Position(P, P[j]-P[i])) then
Add(KK1, -2*To[p+Position(P, P[i])][Position(P,P[j])][2][1]);
Add(KK1, Position(P, P[j]-P[i]) );
	else
Add(KK1, -2*To[p+Position(P, P[i])][Position(P,P[j])][2][1]);
Add(KK1, Position(P, K[j]-K[i]) );
	fi;
fi; fi; fi; fi;
SetEntrySCTable( TT, p+i,  p+j , KK1);
fi;

else
if (t[j]=1) then

if (P[i]+P[j] in P ) then #qui so che aplha+beta < phi_aplha+phi_beta quindi sono nel caso sicuro
if not (IsEmpty(To[Position(P, P[i])][Position(P,P[j])][2])=true) then
if (t[Position(P,P[j]+P[i])]=-1) then
Add(KK1, -2*To[Position(P, P[i])][Position(P,P[j])][2][1]);
Add(KK1, Position(P, P[i]+P[j]) );
fi; fi; fi;

if (P[j]-P[i] in P ) then
if not (IsEmpty(To[p+Position(P, P[i])][Position(P,P[j])][2])=true) then
if (t[Position(P,P[j]-P[i])]=-1) then
if (p+Position(P,K[j]-K[i]) > p+Position(P,P[j]-P[i])) then
Add(KK1, -2*To[p+Position(P, P[i])][Position(P,P[j])][2][1]);
Add(KK1, Position(P, P[j]-P[i]) );
else 
	if (p+Position(P,K[j]-K[i]) = p+Position(P, P[j]-P[i])) then
Add(KK1, -2*To[p+Position(P, P[i])][Position(P,P[j])][2][1]);
Add(KK1, Position(P, P[j]-P[i]) );
	else
Add(KK1, -2*To[p+Position(P, P[i])][Position(P,P[j])][2][1]);
Add(KK1, Position(P, K[j]-K[i]) );
	fi;
fi; fi; fi; fi;

SetEntrySCTable( TT, p+i,  p+j , KK1);


else


if (P[i]+P[j] in P ) then #qui so che aplha+beta < phi_aplha+phi_beta quindi sono nel caso sicuro
if not (IsEmpty(To[Position(P, P[i])][Position(P,P[j])][2])=true) then
if (t[Position(P,P[j]+P[i])]=1) then
Add(KK1, +2*To[Position(P, P[i])][Position(P,P[j])][2][1]);
Add(KK1, Position(P, P[i]+P[j]) );
fi; fi; fi;

if (P[j]-P[i] in P ) then
if not (IsEmpty(To[p+Position(P, P[i])][Position(P,P[j])][2])=true) then
if (t[Position(P,P[j]-P[i])]=1) then
if (p+Position(P,K[j]-K[i]) > p+Position(P,P[j]-P[i])) then
Add(KK1, +2*To[p+Position(P, P[i])][Position(P,P[j])][2][1]);
Add(KK1, Position(P, P[j]-P[i]) );
else 
	if (p+Position(P,K[j]-K[i]) = p+Position(P, P[j]-P[i])) then
Add(KK1, +2*To[p+Position(P, P[i])][Position(P,P[j])][2][1]);
Add(KK1, Position(P, P[j]-P[i]) );
	else
Add(KK1, +2*To[p+Position(P, P[i])][Position(P,P[j])][2][1]);
Add(KK1, Position(P, K[j]-K[i]) );
	fi;
fi; fi; fi; fi;

SetEntrySCTable( TT, p+i,  p+j , KK1);

fi;


fi; fi; fi;

fi; fi;


od; od;

for i in [1..p] do #variazione alpha per gli X
		if ( Position(P, K[i]) > Position(P,P[i]) ) then #alpha< phi_alpha (ragiono con gli X)
KK1:=[];

if (t[i] = 1) then

if (P[i]+K[i] in P ) then
if not  (IsEmpty(To[i][Position(P,K[i])][2] ) = true ) then
if (t[Position(P, P[i]+K[i])]=-1) then
Add(KK1, +1*To[Position(P, P[i])][Position(P,K[i])][2][1]);
Add(KK1, Position(P, K[i]+P[i]) );
fi; fi; fi;

else

if (P[i]+K[i] in P ) then
if not  (IsEmpty(To[i][Position(P,K[i])][2] ) = true ) then
if (t[Position(P, P[i]+K[i])]=-1) then
Add(KK1, -1*To[Position(P, P[i])][Position(P,K[i])][2][1]);
Add(KK1, Position(P, K[i]+P[i]) );
fi; fi; fi;

fi;

SetEntrySCTable( TT, p+i,  p+Position(P,K[i]) , KK1);

		fi;
od;

################################################################################################
########################################################################################################

#############################################################################################################################################
###################################################################################################################

#Prodotti [X,Y] con indici differenti X>Y ovvero alpha > beta !!!!!!!


for i in [2..p] do #variazione alpha per gli X
	for j in [1..i-1] do #variazione beta < alpha per gli Y
		if ( Position(P, K[i]) > Position(P,P[i]) ) then #alpha< phi_alpha (ragiono con gli X)
			if ( p+Position(P, K[j]) > p+Position(P,P[j]) ) then # beta< phi_beta (ragiono con gli Y)

KK1:=[]; # i p+j
KK2:=[]; # i p+jbar
KK3:=[]; # ibar p+j
KK4:=[]; # ibar p+jbar

if (t[i]=1) then
if (t[j]=1) then

if (P[i]+P[j] in P ) then #qui so che aplha+beta < phi_aplha+phi_beta quindi sono nel caso sicuro
if not (IsEmpty(To[Position(P, P[i])][Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]+K[i]) > p+Position(P,P[j]+P[i])) then

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, p+Position(P, P[i]+P[j]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]+K[j]) );
Add(KK3, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]+K[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+P[j]) );
else
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, p+Position(P, P[i]+P[j]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]+K[j]) );
Add(KK3, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]+K[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+P[j]) );
fi; 
else
	if (p+Position(P,K[j]+K[i]) = p+Position(P,P[j]+P[i])) then
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, p+Position(P, P[i]+P[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+P[j]) );
else
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]+K[j]) );
Add(KK3, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]+K[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, p+Position(P, K[j]+K[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, P[i]+P[j]) );
Add(KK3, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, P[i]+P[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P,K[j]+K[i]) );
else
Add(KK1, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, p+Position(P, K[j]+K[i]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, P[i]+P[j]) );
Add(KK3, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, P[i]+P[j]) );
Add(KK4, +1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, K[j]+K[i]) );
fi; 
	fi; 
fi; 
fi; fi;

###############

if (P[i]-P[j] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][p+Position(P,P[j])][2])=true) then
if (p+Position(P,K[i]-K[j]) > p+Position(P,P[i]-P[j])) then

if (t[Position(P,P[i]-P[j] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK1, p+Position(P, P[i]-P[j]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]-K[j]) );
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]-K[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]-P[j]) );
else
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK1, p+Position(P, P[i]-P[j]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]-K[j]) );
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]-K[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]-P[j]) );
fi; 
else
	if (p+Position(P,K[i]-K[j]) = p+Position(P,P[i]-P[j])) then
if (t[Position(P,P[i]-P[j] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK1, p+Position(P, P[i]-P[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]-P[j]) );
else
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]-K[j]) );
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]-K[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[i]-P[j] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK1, p+Position(P, K[i]-K[j]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, P[i]-P[j]) );
Add(KK3, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, P[i]-P[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, K[i]-K[j]) );
else
Add(KK1, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK1, p+Position(P, K[i]-K[j]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, P[i]-P[j]) );
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, P[i]-P[j]) );
Add(KK4, +1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, K[i]-K[j]) );
fi; 
	fi; 
fi; 
fi; fi;

#############

if (P[i]+K[j] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][Position(P,K[j])][2])=true) then
if (p+Position(P,K[i]+P[j]) > p+Position(P,P[i]+K[j])) then

if (t[Position(P,P[i]+K[j] )] = 1) then
Add(KK1, +1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, p+Position(P, P[i]+K[j]) );
Add(KK2, +1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]+P[j]) );
Add(KK3, +1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]+P[j]) );
Add(KK4, -1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+K[j]) );
else
Add(KK1, +1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, p+Position(P, P[i]+K[j]) );
Add(KK2, +1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]+P[j]) );
Add(KK3, +1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]+P[j]) );
Add(KK4, -1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+K[j]) );
fi; 
else
	if (p+Position(P,K[i]+P[j]) = p+Position(P,P[i]+K[j])) then
if (t[Position(P,K[j]+P[i] )] = 1) then
Add(KK1, +1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, p+Position(P, P[i]+K[j]) );
Add(KK4, -1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+K[j]) );
else
Add(KK2, +1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]+P[j]) );
Add(KK3, +1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]+P[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[i]+K[j] )] = 1) then
Add(KK1, +1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, p+Position(P, K[i]+P[j]) );
Add(KK2, -1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, P[i]+K[j]) );
Add(KK3, -1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, P[i]+K[j]) );
Add(KK4, -1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, K[i]+P[j]) );
else
Add(KK1, -1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, p+Position(P, K[i]+P[j]) );
Add(KK2, +1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, P[i]+K[j]) );
Add(KK3, +1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, P[i]+K[j]) );
Add(KK4, +1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, K[i]+P[j]) );
fi; 
	fi; 
fi; 
fi; fi;

#############

if (P[i]-K[j] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][p+Position(P,K[j])][2])=true) then
if (p+Position(P,K[i]-P[j]) > p+Position(P,P[i]-K[j])) then

if (t[Position(P,P[i]-K[j] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK1, p+Position(P, P[i]-K[j]) );
Add(KK2, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]-P[j]) );
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]-P[j]) );
Add(KK4, -1*(-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]-K[j]) );
else
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK1, p+Position(P, P[i]-K[j]) );
Add(KK2, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]-P[j]) );
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]-P[j]) );
Add(KK4, -1*(-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]-K[j]) );
fi; 
else
	if (p+Position(P,K[i]-P[j]) = p+Position(P,P[i]-K[j])) then
if (t[Position(P,P[i]-K[j] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK1, p+Position(P, P[i]-K[j]) );
Add(KK4, -1*(-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]-K[j]) );
else
Add(KK2, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]-P[j]) );
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]-P[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[i]-K[j] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK1, p+Position(P, K[i]-P[j]) );
Add(KK2, -1*(-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, P[i]-K[j]) );
Add(KK3, -1*(+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, P[i]-K[j]) );
Add(KK4, -1*(-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, K[i]-P[j]) );
else
Add(KK1, -1*(+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK1, p+Position(P, K[i]-P[j]) );
Add(KK2, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, P[i]-K[j]) );
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, P[i]-K[j]) );
Add(KK4, -1*(-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, K[i]-P[j]) );
fi; 
	fi; 
fi; 
fi; fi;

#############
else # ovvero t[j]=-1
#########################

if (P[i]+P[j] in P ) then #qui so che aplha+beta < phi_aplha+phi_beta quindi sono nel caso sicuro
if not (IsEmpty(To[Position(P, P[i])][Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]+K[i]) > p+Position(P,P[j]+P[i])) then

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, p+Position(P, P[i]+P[j]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]+K[j]) );
Add(KK3, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]+K[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+P[j]) );
else
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, p+Position(P, P[i]+P[j]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]+K[j]) );
Add(KK3, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]+K[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+P[j]) );
fi; 
else
	if (p+Position(P,K[j]+K[i]) = p+Position(P,P[j]+P[i])) then
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, p+Position(P, P[i]+P[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+P[j]) );
else
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]+K[j]) );
Add(KK3, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]+K[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, p+Position(P, K[j]+K[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, P[i]+P[j]) );
Add(KK3, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, P[i]+P[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P,K[j]+K[i]) );
else
Add(KK1, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, p+Position(P, K[j]+K[i]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, P[i]+P[j]) );
Add(KK3, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, P[i]+P[j]) );
Add(KK4, +1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, K[j]+K[i]) );
fi; 
	fi; 
fi; 
fi; fi;

###############

if (P[i]-P[j] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][p+Position(P,P[j])][2])=true) then
if (p+Position(P,K[i]-K[j]) > p+Position(P,P[i]-P[j])) then

if (t[Position(P,P[i]-P[j] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK1, p+Position(P, P[i]-P[j]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]-K[j]) );
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]-K[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]-P[j]) );
else
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK1, p+Position(P, P[i]-P[j]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]-K[j]) );
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]-K[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]-P[j]) );
fi; 
else
	if (p+Position(P,K[i]-K[j]) = p+Position(P,P[i]-P[j])) then
if (t[Position(P,P[i]-P[j] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK1, p+Position(P, P[i]-P[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]-P[j]) );
else
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]-K[j]) );
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]-K[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[i]-P[j] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK1, p+Position(P, K[i]-K[j]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, P[i]-P[j]) );
Add(KK3, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, P[i]-P[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, K[i]-K[j]) );
else
Add(KK1, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK1, p+Position(P, K[i]-K[j]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, P[i]-P[j]) );
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, P[i]-P[j]) );
Add(KK4, +1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, K[i]-K[j]) );
fi; 
	fi; 
fi; 
fi; fi;

#############

if (P[i]+K[j] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][Position(P,K[j])][2])=true) then
if (p+Position(P,K[i]+P[j]) > p+Position(P,P[i]+K[j])) then

if (t[Position(P,P[i]+K[j] )] = 1) then
Add(KK1, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, p+Position(P, P[i]+K[j]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]+P[j]) );
Add(KK3, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]+P[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+K[j]) );
else
Add(KK1, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, p+Position(P, P[i]+K[j]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]+P[j]) );
Add(KK3, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]+P[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+K[j]) );
fi; 
else
	if (p+Position(P,K[i]+P[j]) = p+Position(P,P[i]+K[j])) then
if (t[Position(P,K[j]+P[i] )] = 1) then
Add(KK1, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, p+Position(P, P[i]+K[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+K[j]) );
else
Add(KK2, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]+P[j]) );
Add(KK3, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]+P[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[i]+K[j] )] = 1) then
Add(KK1, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, p+Position(P, K[i]+P[j]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, P[i]+K[j]) );
Add(KK3, -1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, P[i]+K[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, K[i]+P[j]) );
else
Add(KK1, -1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, p+Position(P, K[i]+P[j]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, P[i]+K[j]) );
Add(KK3, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, P[i]+K[j]) );
Add(KK4, +1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, K[i]+P[j]) );
fi; 
	fi; 
fi; 
fi; fi;

#############

if (P[i]-K[j] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][p+Position(P,K[j])][2])=true) then
if (p+Position(P,K[i]-P[j]) > p+Position(P,P[i]-K[j])) then

if (t[Position(P,P[i]-K[j] )] = 1) then
Add(KK1, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK1, p+Position(P, P[i]-K[j]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]-P[j]) );
Add(KK3, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]-P[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]-K[j]) );
else
Add(KK1, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK1, p+Position(P, P[i]-K[j]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]-P[j]) );
Add(KK3, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]-P[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]-K[j]) );
fi; 
else
	if (p+Position(P,K[i]-P[j]) = p+Position(P,P[i]-K[j])) then
if (t[Position(P,P[i]-K[j] )] = 1) then
Add(KK1, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK1, p+Position(P, P[i]-K[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]-K[j]) );
else
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]-P[j]) );
Add(KK3, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]-P[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[i]-K[j] )] = 1) then
Add(KK1, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK1, p+Position(P, K[i]-P[j]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, P[i]-K[j]) );
Add(KK3, -1*(-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, P[i]-K[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, K[i]-P[j]) );
else
Add(KK1, -1*(-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK1, p+Position(P, K[i]-P[j]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, P[i]-K[j]) );
Add(KK3, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, P[i]-K[j]) );
Add(KK4, +1*(+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, K[i]-P[j]) );
fi; 
	fi; 
fi; 
fi; fi;

fi;

########################à
else #t[i]=-1
if (t[j]=1) then
###############

if (P[i]+P[j] in P ) then #qui so che aplha+beta < phi_aplha+phi_beta quindi sono nel caso sicuro
if not (IsEmpty(To[Position(P, P[i])][Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]+K[i]) > p+Position(P,P[j]+P[i])) then

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, p+Position(P, P[i]+P[j]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]+K[j]) );
Add(KK3, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]+K[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+P[j]) );
else
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, p+Position(P, P[i]+P[j]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]+K[j]) );
Add(KK3, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]+K[j]) );
Add(KK4, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+P[j]) );
fi; 
else
	if (p+Position(P,K[j]+K[i]) = p+Position(P,P[j]+P[i])) then
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, p+Position(P, P[i]+P[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+P[j]) );
else
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]+K[j]) );
Add(KK3, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]+K[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, p+Position(P, K[j]+K[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, P[i]+P[j]) );
Add(KK3, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, P[i]+P[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P,K[j]+K[i]) );
else
Add(KK1, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, p+Position(P, K[j]+K[i]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, P[i]+P[j]) );
Add(KK3, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, P[i]+P[j]) );
Add(KK4, +1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, K[j]+K[i]) );
fi; 
	fi; 
fi; 
fi; fi;
###############

if (P[i]-P[j] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][p+Position(P,P[j])][2])=true) then
if (p+Position(P,K[i]-K[j]) > p+Position(P,P[i]-P[j])) then

if (t[Position(P,P[i]-P[j] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK1, p+Position(P, P[i]-P[j]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]-K[j]) );
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]-K[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]-P[j]) );
else
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK1, p+Position(P, P[i]-P[j]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]-K[j]) );
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]-K[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]-P[j]) );
fi; 
else
	if (p+Position(P,K[i]-K[j]) = p+Position(P,P[i]-P[j])) then
if (t[Position(P,P[i]-P[j] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK1, p+Position(P, P[i]-P[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]-P[j]) );
else
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]-K[j]) );
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]-K[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[i]-P[j] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK1, p+Position(P, K[i]-K[j]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, P[i]-P[j]) );
Add(KK3, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, P[i]-P[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, K[i]-K[j]) );
else
Add(KK1, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK1, p+Position(P, K[i]-K[j]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, P[i]-P[j]) );
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, P[i]-P[j]) );
Add(KK4, +1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, K[i]-K[j]) );
fi; 
	fi; 
fi; 
fi; fi;

#############

if (P[i]+K[j] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][Position(P,K[j])][2])=true) then
if (p+Position(P,K[i]+P[j]) > p+Position(P,P[i]+K[j])) then

if (t[Position(P,P[i]+K[j] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, p+Position(P, P[i]+K[j]) );
Add(KK2, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]+P[j]) );
Add(KK3, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]+P[j]) );
Add(KK4, -1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+K[j]) );
else
Add(KK1, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, p+Position(P, P[i]+K[j]) );
Add(KK2, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]+P[j]) );
Add(KK3, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]+P[j]) );
Add(KK4, -1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+K[j]) );
fi; 
else
	if (p+Position(P,K[i]+P[j]) = p+Position(P,P[i]+K[j])) then
if (t[Position(P,K[j]+P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, p+Position(P, P[i]+K[j]) );
Add(KK4, -1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+K[j]) );
else
Add(KK2, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]+P[j]) );
Add(KK3, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]+P[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[i]+K[j] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, p+Position(P, K[i]+P[j]) );
Add(KK2, -1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, P[i]+K[j]) );
Add(KK3, -1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, P[i]+K[j]) );
Add(KK4, -1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, K[i]+P[j]) );
else
Add(KK1, -1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, p+Position(P, K[i]+P[j]) );
Add(KK2, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, P[i]+K[j]) );
Add(KK3, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, P[i]+K[j]) );
Add(KK4, +1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, K[i]+P[j]) );
fi; 
	fi; 
fi; 
fi; fi;

#############

if (P[i]-K[j] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][p+Position(P,K[j])][2])=true) then
if (p+Position(P,K[i]-P[j]) > p+Position(P,P[i]-K[j])) then

if (t[Position(P,P[i]-K[j] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK1, p+Position(P, P[i]-K[j]) );
Add(KK2, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]-P[j]) );
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]-P[j]) );
Add(KK4, -1*(-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]-K[j]) );
else
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK1, p+Position(P, P[i]-K[j]) );
Add(KK2, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]-P[j]) );
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]-P[j]) );
Add(KK4, -1*(-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]-K[j]) );
fi; 
else
	if (p+Position(P,K[i]-P[j]) = p+Position(P,P[i]-K[j])) then
if (t[Position(P,P[i]-K[j] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK1, p+Position(P, P[i]-K[j]) );
Add(KK4, -1*(-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]-K[j]) );
else
Add(KK2, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]-P[j]) );
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]-P[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[i]-K[j] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK1, p+Position(P, K[i]-P[j]) );
Add(KK2, -1*(-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, P[i]-K[j]) );
Add(KK3, -1*(+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, P[i]-K[j]) );
Add(KK4, -1*(-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, K[i]-P[j]) );
else
Add(KK1, -1*(+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK1, p+Position(P, K[i]-P[j]) );
Add(KK2, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, P[i]-K[j]) );
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, P[i]-K[j]) );
Add(KK4, +1*(-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, K[i]-P[j]) );
fi; 
	fi; 
fi; 
fi; fi;

#####################
else #t[j]=-1
#########################

if (P[i]+P[j] in P ) then #qui so che aplha+beta < phi_aplha+phi_beta quindi sono nel caso sicuro
if not (IsEmpty(To[Position(P, P[i])][Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]+K[i]) > p+Position(P,P[j]+P[i])) then

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, p+Position(P, P[i]+P[j]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]+K[j]) );
Add(KK3, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]+K[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+P[j]) );
else
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, p+Position(P, P[i]+P[j]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]+K[j]) );
Add(KK3, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]+K[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+P[j]) );
fi; 
else
	if (p+Position(P,K[j]+K[i]) = p+Position(P,P[j]+P[i])) then
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, p+Position(P, P[i]+P[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+P[j]) );
else
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]+K[j]) );
Add(KK3, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]+K[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, p+Position(P, K[j]+K[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, P[i]+P[j]) );
Add(KK3, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, P[i]+P[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P,K[j]+K[i]) );
else
Add(KK1, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, p+Position(P, K[j]+K[i]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, P[i]+P[j]) );
Add(KK3, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, P[i]+P[j]) );
Add(KK4, +1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, K[j]+K[i]) );
fi; 
	fi; 
fi; 
fi; fi;

###############

if (P[i]-P[j] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][p+Position(P,P[j])][2])=true) then
if (p+Position(P,K[i]-K[j]) > p+Position(P,P[i]-P[j])) then

if (t[Position(P,P[i]-P[j] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK1, p+Position(P, P[i]-P[j]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]-K[j]) );
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]-K[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]-P[j]) );
else
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK1, p+Position(P, P[i]-P[j]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]-K[j]) );
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]-K[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]-P[j]) );
fi; 
else
	if (p+Position(P,K[i]-K[j]) = p+Position(P,P[i]-P[j])) then
if (t[Position(P,P[i]-P[j] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK1, p+Position(P, P[i]-P[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]-P[j]) );
else
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]-K[j]) );
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]-K[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[i]-P[j] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK1, p+Position(P, K[i]-K[j]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, P[i]-P[j]) );
Add(KK3, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, P[i]-P[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, K[i]-K[j]) );
else
Add(KK1, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK1, p+Position(P, K[i]-K[j]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, P[i]-P[j]) );
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, P[i]-P[j]) );
Add(KK4, +1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, K[i]-K[j]) );
fi; 
	fi; 
fi; 
fi; fi;

#############

if (P[i]+K[j] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][Position(P,K[j])][2])=true) then
if (p+Position(P,K[i]+P[j]) > p+Position(P,P[i]+K[j])) then

if (t[Position(P,P[i]+K[j] )] = 1) then
Add(KK1, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, p+Position(P, P[i]+K[j]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]+P[j]) );
Add(KK3, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]+P[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+K[j]) );
else
Add(KK1, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, p+Position(P, P[i]+K[j]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]+P[j]) );
Add(KK3, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]+P[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+K[j]) );
fi; 
else
	if (p+Position(P,K[i]+P[j]) = p+Position(P,P[i]+K[j])) then
if (t[Position(P,K[j]+P[i] )] = 1) then
Add(KK1, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, p+Position(P, P[i]+K[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+K[j]) );
else
Add(KK2, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]+P[j]) );
Add(KK3, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]+P[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[i]+K[j] )] = 1) then
Add(KK1, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, p+Position(P, K[i]+P[j]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, P[i]+K[j]) );
Add(KK3, -1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, P[i]+K[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, K[i]+P[j]) );
else
Add(KK1, -1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, p+Position(P, K[i]+P[j]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, P[i]+K[j]) );
Add(KK3, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, P[i]+K[j]) );
Add(KK4, +1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, K[i]+P[j]) );
fi; 
	fi; 
fi; 
fi; fi;

#############

if (P[i]-K[j] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][p+Position(P,K[j])][2])=true) then
if (p+Position(P,K[i]-P[j]) > p+Position(P,P[i]-K[j])) then

if (t[Position(P,P[i]-K[j] )] = 1) then
Add(KK1, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK1, p+Position(P, P[i]-K[j]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]-P[j]) );
Add(KK3, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]-P[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]-K[j]) );
else
Add(KK1, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK1, p+Position(P, P[i]-K[j]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]-P[j]) );
Add(KK3, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]-P[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]-K[j]) );
fi; 
else
	if (p+Position(P,K[i]-P[j]) = p+Position(P,P[i]-K[j])) then
if (t[Position(P,P[i]-K[j] )] = 1) then
Add(KK1, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK1, p+Position(P, P[i]-K[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]-K[j]) );
else
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]-P[j]) );
Add(KK3, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]-P[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[i]-K[j] )] = 1) then
Add(KK1, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK1, p+Position(P, K[i]-P[j]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, P[i]-K[j]) );
Add(KK3, -1*(-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, P[i]-K[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, K[i]-P[j]) );
else
Add(KK1, -1*(-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK1, p+Position(P, K[i]-P[j]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, P[i]-K[j]) );
Add(KK3, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, P[i]-K[j]) );
Add(KK4, +1*(+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, K[i]-P[j]) );
fi; 
	fi; 
fi; 
fi; fi;

fi; fi;

if ( Position(P,P[i]) > Position(P,K[j]) ) then
SetEntrySCTable( TT, Position(P,P[i]),  p+Position(P,P[j]) , KK1);
SetEntrySCTable( TT, Position(P,P[i]),  p+Position(P, K[j]), KK2);
SetEntrySCTable( TT, Position(P,K[i]),  p+Position(P,P[j]), KK3);
SetEntrySCTable( TT, Position(P,K[i]),  p+Position(P, K[j]), KK4);
fi;

if ( Position(P,K[j]) > Position(P,K[i]) ) then
SetEntrySCTable( TT, Position(P,P[i]),  p+Position(P,P[j]) , KK1);
SetEntrySCTable( TT, Position(P,K[i]),  p+Position(P,P[j]), KK3);
fi;

if ( Position(P,K[j]) < Position(P,K[i]) ) then
	if ( Position(P,P[i]) < Position(P,K[j]) ) then
SetEntrySCTable( TT, Position(P,P[i]),  p+Position(P,P[j]) , KK1);
SetEntrySCTable( TT, Position(P,K[i]),  p+Position(P,P[j]), KK3);
SetEntrySCTable( TT, Position(P,K[i]),  p+Position(P, K[j]), KK4);
	fi;
fi;

else

	if ( p+Position(P, K[j]) = p+Position(P,P[j]) ) then #beta = phi_beta

KK1:=[];
KK2:=[];

if (t[i]=1) then
if (t[j]=1) then

if (P[i]+P[j] in P ) then #qui so che aplha+beta < phi_aplha+phi_beta quindi sono nel caso sicuro
if not (IsEmpty(To[Position(P, P[i])][Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]+K[i]) > p+Position(P,P[j]+P[i])) then

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, P[i]+P[j]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[i]+K[j]) );
else
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, P[i]+P[j]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[i]+K[j]) );
fi; 
else
	if (p+Position(P,K[j]+K[i]) = p+Position(P,P[j]+P[i])) then
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, P[i]+P[j]) );
else
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[i]+K[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, K[j]+K[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[i]+P[j]) );
else
Add(KK1, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, K[j]+K[i]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[i]+P[j]) );
fi; 
	fi; 
fi; 
fi; fi;

###############

if (P[i]-P[j] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][p+Position(P,P[j])][2])=true) then
if (p+Position(P,K[i]-K[j]) > p+Position(P,P[i]-P[j])) then

if (t[Position(P,P[i]-P[j] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, P[i]-P[j]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[i]-K[j]) );
else
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, P[i]-P[j]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[i]-K[j]) );
fi; 
else
	if (p+Position(P,K[i]-K[j]) = p+Position(P,P[i]-P[j])) then
if (t[Position(P,P[i]-P[j] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, P[i]-P[j]) );
else
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[i]-K[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[i]-P[j] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, K[i]-K[j]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[i]-P[j]) );
else
Add(KK1, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, K[i]-K[j]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[i]-P[j]) );
fi; 
	fi; 
fi; 
fi; fi;

####################
else #t[j]=-1

if (P[i]+P[j] in P ) then #qui so che aplha+beta < phi_aplha+phi_beta quindi sono nel caso sicuro
if not (IsEmpty(To[Position(P, P[i])][Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]+K[i]) > p+Position(P,P[j]+P[i])) then

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, K[j]+K[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[i]+P[j]) );
else
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, K[j]+K[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[i]+P[j]) );
fi; 
else
	if (p+Position(P,K[j]+K[i]) = p+Position(P,P[j]+P[i])) then
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[i]+P[j]) );
else
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, K[j]+K[i]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, K[j]+K[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[i]+P[j]) );
else
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, P[i]+P[j]) );
Add(KK2, +1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[j]+K[i]) );
fi; 
	fi; 
fi; 
fi; fi;

###############

if (P[i]-P[j] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][p+Position(P,P[j])][2])=true) then
if (p+Position(P,K[i]-K[j]) > p+Position(P,P[i]-P[j])) then

if (t[Position(P,P[i]-P[j] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, K[i]-K[j]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[i]-P[j]) );
else
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, K[i]-K[j]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[i]-P[j]) );
fi; 
else
	if (p+Position(P,K[i]-K[j]) = p+Position(P,P[i]-P[j])) then
if (t[Position(P,P[i]-P[j] )] = 1) then
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[i]-P[j]) );
else
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, K[i]-K[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[i]-P[j] )] = 1) then
Add(KK1, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, P[i]-P[j]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[i]-K[j]) );
else
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, P[i]-P[j]) );
Add(KK2, +1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[i]-K[j]) );
fi; 
	fi; 
fi; 
fi; fi;

fi;

else #t[i]=-1
if (t[j]=1) then

if (P[i]+P[j] in P ) then #qui so che aplha+beta < phi_aplha+phi_beta quindi sono nel caso sicuro
if not (IsEmpty(To[Position(P, P[i])][Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]+K[i]) > p+Position(P,P[j]+P[i])) then

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, P[i]+P[j]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[i]+K[j]) );
else
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, P[i]+P[j]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[i]+K[j]) );
fi; 
else
	if (p+Position(P,K[j]+K[i]) = p+Position(P,P[j]+P[i])) then
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, P[i]+P[j]) );
else
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[i]+K[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, K[j]+K[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[i]+P[j]) );
else
Add(KK1, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, K[j]+K[i]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[i]+P[j]) );
fi; 
	fi; 
fi; 
fi; fi;

###############

if (P[i]-P[j] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][p+Position(P,P[j])][2])=true) then
if (p+Position(P,K[i]-K[j]) > p+Position(P,P[i]-P[j])) then

if (t[Position(P,P[i]-P[j] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, P[i]-P[j]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[i]-K[j]) );
else
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, P[i]-P[j]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[i]-K[j]) );
fi; 
else
	if (p+Position(P,K[i]-K[j]) = p+Position(P,P[i]-P[j])) then
if (t[Position(P,P[i]-P[j] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, P[i]-P[j]) );
else
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[i]-K[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[i]-P[j] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, K[i]-K[j]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[i]-P[j]) );
else
Add(KK1, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, K[i]-K[j]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[i]-P[j]) );
fi; 
	fi; 
fi; 
fi; fi;

####################
else #t[i]=t[j]=-1

if (P[i]+P[j] in P ) then #qui so che aplha+beta < phi_aplha+phi_beta quindi sono nel caso sicuro
if not (IsEmpty(To[Position(P, P[i])][Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]+K[i]) > p+Position(P,P[j]+P[i])) then

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, K[i]+K[j]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[j]+P[i]) );
else
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, K[i]+K[j]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[j]+P[i]) );
fi; 
else
	if (p+Position(P,K[j]+K[i]) = p+Position(P,P[j]+P[i])) then
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[j]+P[i]) );
else
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, K[i]+K[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, P[j]+P[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[i]+K[j]) );
else
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, P[j]+P[i]) );
Add(KK2, +1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[i]+K[j]) );
fi; 
	fi; 
fi; 
fi; fi;

###############

if (P[i]-P[j] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][p+Position(P,P[j])][2])=true) then
if (p+Position(P,K[i]-K[j]) > p+Position(P,P[i]-P[j])) then

if (t[Position(P,P[i]-P[j] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, K[i]-K[j]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[i]-P[j]) );
else
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, K[i]-K[j]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[i]-P[j]) );
fi; 
else
	if (p+Position(P,K[i]-K[j]) = p+Position(P,P[i]-P[j])) then
if (t[Position(P,P[i]-P[j] )] = 1) then
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[i]-P[j]) );
else
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, K[i]-K[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[i]-P[j] )] = 1) then
Add(KK1, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, P[i]-P[j]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[i]-K[j]) );
else
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, P[i]-P[j]) );
Add(KK2, +1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[i]-K[j]) );
fi; 
	fi; 
fi; 
fi; fi;

fi; fi;

#if (Position(P,P[j]) > Position(P,K[i]) ) then
SetEntrySCTable( TT, i,  p+j , KK1);
SetEntrySCTable( TT, Position(P,K[i]), p+j, KK2);



fi;
fi;

else
	if ( p+Position(P, K[i]) = p+Position(P,P[i]) ) then #alpha< phi_alpha
			if ( p+Position(P, K[j]) > p+Position(P,P[j]) ) then #beta< phi_beta

KK1:=[];
KK2:=[];

if (t[i]=1) then
if (t[j]=1) then

if (P[i]+P[j] in P ) then #qui so che aplha+beta < phi_aplha+phi_beta quindi sono nel caso sicuro
if not (IsEmpty(To[Position(P, P[i])][Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]+K[i]) > p+Position(P,P[j]+P[i])) then

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, P[i]+P[j]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[i]+K[j]) );
else
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, P[i]+P[j]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[i]+K[j]) );
fi; 
else
	if (p+Position(P,K[j]+K[i]) = p+Position(P,P[j]+P[i])) then
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, P[i]+P[j]) );
else
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[i]+K[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, K[j]+K[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[i]+P[j]) );
else
Add(KK1, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, K[j]+K[i]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[i]+P[j]) );
fi; 
	fi; 
fi; 
fi; fi;

###############

if (P[i]-P[j] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][p+Position(P,P[j])][2])=true) then
if (p+Position(P,K[i]-K[j]) > p+Position(P,P[i]-P[j])) then

if (t[Position(P,P[i]-P[j] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, P[i]-P[j]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[i]-K[j]) );
else
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, P[i]-P[j]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[i]-K[j]) );
fi; 
else
	if (p+Position(P,K[i]-K[j]) = p+Position(P,P[i]-P[j])) then
if (t[Position(P,P[i]-P[j] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, P[i]-P[j]) );
else
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[i]-K[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[i]-P[j] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, K[i]-K[j]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[i]-P[j]) );
else
Add(KK1, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, K[i]-K[j]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[i]-P[j]) );
fi; 
	fi; 
fi; 
fi; fi;

####################
else #t[j]=-1

if (P[i]+P[j] in P ) then #qui so che aplha+beta < phi_aplha+phi_beta quindi sono nel caso sicuro
if not (IsEmpty(To[Position(P, P[i])][Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]+K[i]) > p+Position(P,P[j]+P[i])) then

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, P[i]+P[j]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[i]+K[j]) );
else
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, P[i]+P[j]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[i]+K[j]) );
fi; 
else
	if (p+Position(P,K[j]+K[i]) = p+Position(P,P[j]+P[i])) then
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, P[i]+P[j]) );
else
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[i]+K[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, K[j]+K[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[i]+P[j]) );
else
Add(KK1, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, K[j]+K[i]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[i]+P[j]) );
fi; 
	fi; 
fi; 
fi; fi;

###############

if (P[i]-P[j] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][p+Position(P,P[j])][2])=true) then
if (p+Position(P,K[i]-K[j]) > p+Position(P,P[i]-P[j])) then

if (t[Position(P,P[i]-P[j] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, P[i]-P[j]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[i]-K[j]) );
else
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, P[i]-P[j]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[i]-K[j]) );
fi; 
else
	if (p+Position(P,K[i]-K[j]) = p+Position(P,P[i]-P[j])) then
if (t[Position(P,P[i]-P[j] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, P[i]-P[j]) );
else
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[i]-K[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[i]-P[j] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, K[i]-K[j]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[i]-P[j]) );
else
Add(KK1, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, K[i]-K[j]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[i]-P[j]) );
fi; 
	fi; 
fi; 
fi; fi;

####################
fi;

else #t[i]=-1
if (t[j]=1) then

if (P[i]+P[j] in P ) then #qui so che aplha+beta < phi_aplha+phi_beta quindi sono nel caso sicuro
if not (IsEmpty(To[Position(P, P[i])][Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]+K[i]) > p+Position(P,P[j]+P[i])) then

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, K[j]+K[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[i]+P[j]) );
else
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, K[j]+K[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[i]+P[j]) );
fi; 
else
	if (p+Position(P,K[j]+K[i]) = p+Position(P,P[j]+P[i])) then
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[i]+P[j]) );
else
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, K[j]+K[i]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, K[j]+K[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[i]+P[j]) );
else
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, P[i]+P[j]) );
Add(KK2, +1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[j]+K[i]) );
fi; 
	fi; 
fi; 
fi; fi;

###############

if (P[i]-P[j] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][p+Position(P,P[j])][2])=true) then
if (p+Position(P,K[i]-K[j]) > p+Position(P,P[i]-P[j])) then

if (t[Position(P,P[i]-P[j] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, K[i]-K[j]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[i]-P[j]) );
else
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, K[i]-K[j]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[i]-P[j]) );
fi; 
else
	if (p+Position(P,K[i]-K[j]) = p+Position(P,P[i]-P[j])) then
if (t[Position(P,P[i]-P[j] )] = 1) then
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[i]-P[j]) );
else
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, K[i]-K[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[i]-P[j] )] = 1) then
Add(KK1, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, P[i]-P[j]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[i]-K[j]) );
else
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, P[i]-P[j]) );
Add(KK2, +1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[i]-K[j]) );
fi; 
	fi; 
fi; 
fi; fi;

else

if (P[i]+P[j] in P ) then #qui so che aplha+beta < phi_aplha+phi_beta quindi sono nel caso sicuro
if not (IsEmpty(To[Position(P, P[i])][Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]+K[i]) > p+Position(P,P[j]+P[i])) then

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, K[i]+K[j]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[j]+P[i]) );
else
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, K[i]+K[j]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[j]+P[i]) );
fi; 
else
	if (p+Position(P,K[j]+K[i]) = p+Position(P,P[j]+P[i])) then
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[j]+P[i]) );
else
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, K[i]+K[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, P[j]+P[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[i]+K[j]) );
else
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, P[j]+P[i]) );
Add(KK2, +1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[i]+K[j]) );
fi; 
	fi; 
fi; 
fi; fi;

###############

if (P[i]-P[j] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][p+Position(P,P[j])][2])=true) then
if (p+Position(P,K[i]-K[j]) > p+Position(P,P[i]-P[j])) then

if (t[Position(P,P[i]-P[j] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, K[i]-K[j]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[i]-P[j]) );
else
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, K[i]-K[j]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[i]-P[j]) );
fi; 
else
	if (p+Position(P,K[i]-K[j]) = p+Position(P,P[i]-P[j])) then
if (t[Position(P,P[i]-P[j] )] = 1) then
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[i]-P[j]) );
else
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, K[i]-K[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[i]-P[j] )] = 1) then
Add(KK1, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, P[i]-P[j]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[i]-K[j]) );
else
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, P[i]-P[j]) );
Add(KK2, +1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[i]-K[j]) );
fi; 
	fi; 
fi; 
fi; fi;

fi;
fi;

if ( Position(P, P[i]) > Position(P,K[j]) ) then
SetEntrySCTable( TT, i,  p+j , KK1);
SetEntrySCTable( TT, i,  p+Position(P,K[j]), KK2);
else
	if ( Position(P, P[i]) < Position(P,K[j]) ) then
SetEntrySCTable( TT, i,  p+j , KK1);
	fi;
fi;

#################################################
		else
			if ( p+Position(P, K[j]) = p+Position(P,P[j]) ) then #beta = phi_beta

KK1:=[];

if (t[i]=1) then
if (t[j]=1) then

if (P[i]+P[j] in P ) then #qui so che aplha+beta < phi_aplha + phi_beta quindi sono nel caso sicuro
if not (IsEmpty(To[Position(P, P[i])][Position(P,P[j])][2])=true) then
if (t[Position(P,P[j]+P[i])]=1) then
Add(KK1, +2*To[Position(P, P[i])][Position(P,P[j])][2][1]);
Add(KK1, p+Position(P, P[i]+P[j]) );
fi; fi; fi;

if (P[i]-P[j] in P ) then
if not (IsEmpty(To[Position(P, P[i])][p+Position(P,P[j])][2])=true) then
if (t[Position(P,P[i]-P[j])]=1) then
if (p+Position(P,K[i]-K[j]) > p+Position(P,P[i]-P[j])) then
Add(KK1, +2*To[Position(P, P[i])][p+Position(P,P[j])][2][1]);
Add(KK1, p+Position(P, P[i]-P[j]) );
else 
	if (p+Position(P,K[i]-K[j]) = p+Position(P, P[i]-P[j])) then
Add(KK1, +2*To[Position(P, P[i])][p+Position(P,P[j])][2][1]);
Add(KK1, p+Position(P, P[i]-P[j]) );
	else
Add(KK1, +2*To[Position(P, P[i])][p+Position(P,P[j])][2][1]);
Add(KK1, p+Position(P, K[i]-K[j]) );
	fi;
fi; fi; fi; fi;

SetEntrySCTable( TT, i,  p+j , KK1);

else

if (P[i]+P[j] in P ) then #qui so che aplha+beta < phi_aplha+phi_beta quindi sono nel caso sicuro
if not (IsEmpty(To[Position(P, P[i])][Position(P,P[j])][2])=true) then
if (t[Position(P,P[j]+P[i])]=-1) then
Add(KK1, +2*To[Position(P, P[i])][Position(P,P[j])][2][1]);
Add(KK1, p+Position(P, P[i]+P[j]) );
fi; fi; fi;

if (P[i]-P[j] in P ) then
if not (IsEmpty(To[Position(P, P[i])][p+Position(P,P[j])][2])=true) then
if (t[Position(P,P[i]-P[j])]=-1) then
if (p+Position(P,K[i]-K[j]) > p+Position(P,P[i]-P[j])) then
Add(KK1, +2*To[Position(P, P[i])][p+Position(P,P[j])][2][1]);
Add(KK1, p+Position(P, P[i]-P[j]) );
else 
	if (p+Position(P,K[i]-K[j]) = p+Position(P, P[i]-P[j])) then
Add(KK1, +2*To[Position(P, P[i])][p+Position(P,P[j])][2][1]);
Add(KK1, p+Position(P, P[i]-P[j]) );
	else
Add(KK1, +2*To[Position(P, P[i])][p+Position(P,P[j])][2][1]);
Add(KK1, p+Position(P, K[i]-K[j]) );
	fi;
fi; fi; fi; fi;
SetEntrySCTable( TT, i,  p+j , KK1);
fi;

else
if (t[j]=1) then

if (P[i]+P[j] in P ) then #qui so che aplha+beta < phi_aplha+phi_beta quindi sono nel caso sicuro
if not (IsEmpty(To[Position(P, P[i])][Position(P,P[j])][2])=true) then
if (t[Position(P,P[j]+P[i])]=-1) then
Add(KK1, +2*To[Position(P, P[i])][Position(P,P[j])][2][1]);
Add(KK1, p+Position(P, P[i]+P[j]) );
fi; fi; fi;

if (P[i]-P[j] in P ) then
if not (IsEmpty(To[Position(P, P[i])][p+Position(P,P[j])][2])=true) then
if (t[Position(P,P[i]-P[j])]=-1) then
if (p+Position(P,K[i]-K[j]) > p+Position(P,P[i]-P[j])) then
Add(KK1, +2*To[Position(P, P[i])][p+Position(P,P[j])][2][1]);
Add(KK1, p+Position(P, P[i]-P[j]) );
else 
	if (p+Position(P,K[i]-K[j]) = p+Position(P, P[i]-P[j])) then
Add(KK1, +2*To[Position(P, P[i])][p+Position(P,P[j])][2][1]);
Add(KK1, p+Position(P, P[i]-P[j]) );
	else
Add(KK1, +2*To[Position(P, P[i])][p+Position(P,P[j])][2][1]);
Add(KK1, p+Position(P, K[i]-K[j]) );
	fi;
fi; fi; fi; fi;

SetEntrySCTable( TT, i,  p+j , KK1);


else


if (P[i]+P[j] in P ) then #qui so che aplha+beta < phi_aplha+phi_beta quindi sono nel caso sicuro
if not (IsEmpty(To[Position(P, P[i])][Position(P,P[j])][2])=true) then
if (t[Position(P,P[j]+P[i])]=1) then
Add(KK1, -2*To[Position(P, P[i])][Position(P,P[j])][2][1]);
Add(KK1, p+Position(P, P[i]+P[j]) );
fi; fi; fi;

if (P[i]-P[j] in P ) then
if not (IsEmpty(To[Position(P, P[i])][p+Position(P,P[j])][2])=true) then
if (t[Position(P,P[i]-P[j])]=1) then
if (p+Position(P,K[i]-K[j]) > p+Position(P,P[i]-P[j])) then
Add(KK1, -2*To[Position(P, P[i])][p+Position(P,P[j])][2][1]);
Add(KK1, p+Position(P, P[i]-P[j]) );
else 
	if (p+Position(P,K[i]-K[j]) = p+Position(P, P[i]-P[j])) then
Add(KK1, -2*To[Position(P, P[i])][p+Position(P,P[j])][2][1]);
Add(KK1, p+Position(P, P[i]-P[j]) );
	else
Add(KK1, -2*To[Position(P, P[i])][p+Position(P,P[j])][2][1]);
Add(KK1, p+Position(P, K[i]-K[j]) );
	fi;
fi; fi; fi; fi;

SetEntrySCTable( TT, i,  p+j , KK1);

fi;
fi; fi; fi;
fi; fi;
od; od;

for i in [1..p-1] do #alpha
	for j in [i+1..p] do #beta > alpha
		if ( Position(P,K[i]) > Position(P, P[j]) ) then #quindi alpha<phi_alpha
			if ( Position(P, K[j]) > Position(P, P[j]) ) then
KK3:=[]; # ibar p+j
KK4:=[]; # ibar p+jbar

if (t[i]=1) then
if (t[j]=1) then

if (P[i]+P[j] in P ) then #qui so che aplha+beta < phi_aplha+phi_beta quindi sono nel caso sicuro
if not (IsEmpty(To[Position(P, P[i])][Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]+K[i]) > p+Position(P,P[j]+P[i])) then

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK3, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]+K[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+P[j]) );
else
Add(KK3, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]+K[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+P[j]) );
fi; 
else
	if (p+Position(P,K[j]+K[i]) = p+Position(P,P[j]+P[i])) then
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+P[j]) );
else
Add(KK3, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]+K[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK3, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, P[i]+P[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P,K[j]+K[i]) );
else
Add(KK3, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, P[i]+P[j]) );
Add(KK4, +1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, K[j]+K[i]) );
fi; 
	fi; 
fi; 
fi; fi;

###############

if (P[i]-P[j] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][p+Position(P,P[j])][2])=true) then
if (p+Position(P,K[i]-K[j]) > p+Position(P,P[i]-P[j])) then

if (t[Position(P,P[i]-P[j] )] = 1) then
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]-K[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]-P[j]) );
else
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]-K[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]-P[j]) );
fi; 
else
	if (p+Position(P,K[i]-K[j]) = p+Position(P,P[i]-P[j])) then
if (t[Position(P,P[i]-P[j] )] = 1) then
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]-P[j]) );
else
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]-K[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[i]-P[j] )] = 1) then
Add(KK3, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, P[i]-P[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, K[i]-K[j]) );
else
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, P[i]-P[j]) );
Add(KK4, +1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, K[i]-K[j]) );
fi; 
	fi; 
fi; 
fi; fi;

#############

if (P[i]+K[j] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][Position(P,K[j])][2])=true) then
if (p+Position(P,K[i]+P[j]) > p+Position(P,P[i]+K[j])) then

if (t[Position(P,P[i]+K[j] )] = 1) then
Add(KK3, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]+P[j]) );
Add(KK4, -1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+K[j]) );
else
Add(KK3, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]+P[j]) );
Add(KK4, -1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+K[j]) );
fi; 
else
	if (p+Position(P,K[i]+P[j]) = p+Position(P,P[i]+K[j])) then
if (t[Position(P,K[j]+P[i] )] = 1) then
Add(KK4, -1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+K[j]) );
else
Add(KK3, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]+P[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[i]+K[j] )] = 1) then
Add(KK3, -1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, P[i]+K[j]) );
Add(KK4, -1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, K[i]+P[j]) );
else
Add(KK3, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, P[i]+K[j]) );
Add(KK4, -1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, K[i]+P[j]) );
fi; 
	fi; 
fi; 
fi; fi;

#############

if (P[i]-K[j] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][p+Position(P,K[j])][2])=true) then
if (p+Position(P,K[i]-P[j]) > p+Position(P,P[i]-K[j])) then

if (t[Position(P,P[i]-K[j] )] = 1) then
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]-P[j]) );
Add(KK4, -1*(-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]-K[j]) );
else
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]-P[j]) );
Add(KK4, -1*(-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]-K[j]) );
fi; 
else
	if (p+Position(P,K[i]-P[j]) = p+Position(P,P[i]-K[j])) then
if (t[Position(P,P[i]-K[j] )] = 1) then
Add(KK4, -1*(-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]-K[j]) );
else
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]-P[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[i]-K[j] )] = 1) then
Add(KK3, -1*(+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, P[i]-K[j]) );
Add(KK4, -1*(-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, K[i]-P[j]) );
else
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, P[i]-K[j]) );
Add(KK4, -1*(-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, K[i]-P[j]) );
fi; 
	fi; 
fi; 
fi; fi;

#############
else # ovvero t[j]=-1
#########################

if (P[i]+P[j] in P ) then #qui so che aplha+beta < phi_aplha+phi_beta quindi sono nel caso sicuro
if not (IsEmpty(To[Position(P, P[i])][Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]+K[i]) > p+Position(P,P[j]+P[i])) then

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK3, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]+K[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+P[j]) );
else
Add(KK3, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]+K[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+P[j]) );
fi; 
else
	if (p+Position(P,K[j]+K[i]) = p+Position(P,P[j]+P[i])) then
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+P[j]) );
else
Add(KK3, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]+K[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK3, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, P[i]+P[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P,K[j]+K[i]) );
else
Add(KK3, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, P[i]+P[j]) );
Add(KK4, +1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, K[j]+K[i]) );
fi; 
	fi; 
fi; 
fi; fi;

###############

if (P[i]-P[j] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][p+Position(P,P[j])][2])=true) then
if (p+Position(P,K[i]-K[j]) > p+Position(P,P[i]-P[j])) then

if (t[Position(P,P[i]-P[j] )] = 1) then
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]-K[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]-P[j]) );
else
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]-K[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]-P[j]) );
fi; 
else
	if (p+Position(P,K[i]-K[j]) = p+Position(P,P[i]-P[j])) then
if (t[Position(P,P[i]-P[j] )] = 1) then
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]-P[j]) );
else
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]-K[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[i]-P[j] )] = 1) then
Add(KK3, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, P[i]-P[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, K[i]-K[j]) );
else
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, P[i]-P[j]) );
Add(KK4, +1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, K[i]-K[j]) );
fi; 
	fi; 
fi; 
fi; fi;

#############

if (P[i]+K[j] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][Position(P,K[j])][2])=true) then
if (p+Position(P,K[i]+P[j]) > p+Position(P,P[i]+K[j])) then

if (t[Position(P,P[i]+K[j] )] = 1) then
Add(KK3, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]+P[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+K[j]) );
else
Add(KK3, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]+P[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+K[j]) );
fi; 
else
	if (p+Position(P,K[i]+P[j]) = p+Position(P,P[i]+K[j])) then
if (t[Position(P,K[j]+P[i] )] = 1) then
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+K[j]) );
else
Add(KK3, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]+P[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[i]+K[j] )] = 1) then
Add(KK3, -1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, P[i]+K[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, K[i]+P[j]) );
else
Add(KK3, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, P[i]+K[j]) );
Add(KK4, +1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, K[i]+P[j]) );
fi; 
	fi; 
fi; 
fi; fi;

#############

if (P[i]-K[j] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][p+Position(P,K[j])][2])=true) then
if (p+Position(P,K[i]-P[j]) > p+Position(P,P[i]-K[j])) then

if (t[Position(P,P[i]-K[j] )] = 1) then
Add(KK3, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]-P[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]-K[j]) );
else
Add(KK3, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]-P[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]-K[j]) );
fi; 
else
	if (p+Position(P,K[i]-P[j]) = p+Position(P,P[i]-K[j])) then
if (t[Position(P,P[i]-K[j] )] = 1) then
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]-K[j]) );
else
Add(KK3, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]-P[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[i]-K[j] )] = 1) then
Add(KK3, -1*(-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, P[i]-K[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, K[i]-P[j]) );
else
Add(KK3, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, P[i]-K[j]) );
Add(KK4, +1*(+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, K[i]-P[j]) );
fi; 
	fi; 
fi; 
fi; fi;

fi;

########################à
else #t[i]=-1
if (t[j]=1) then
###############

if (P[i]+P[j] in P ) then #qui so che aplha+beta < phi_aplha+phi_beta quindi sono nel caso sicuro
if not (IsEmpty(To[Position(P, P[i])][Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]+K[i]) > p+Position(P,P[j]+P[i])) then

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK3, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]+K[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+P[j]) );
else
Add(KK3, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]+K[j]) );
Add(KK4, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+P[j]) );
fi; 
else
	if (p+Position(P,K[j]+K[i]) = p+Position(P,P[j]+P[i])) then
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+P[j]) );
else
Add(KK3, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]+K[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK3, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, P[i]+P[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P,K[j]+K[i]) );
else
Add(KK3, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, P[i]+P[j]) );
Add(KK4, +1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, K[j]+K[i]) );
fi; 
	fi; 
fi; 
fi; fi;
###############

if (P[i]-P[j] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][p+Position(P,P[j])][2])=true) then
if (p+Position(P,K[i]-K[j]) > p+Position(P,P[i]-P[j])) then

if (t[Position(P,P[i]-P[j] )] = 1) then
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]-K[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]-P[j]) );
else
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]-K[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]-P[j]) );
fi; 
else
	if (p+Position(P,K[i]-K[j]) = p+Position(P,P[i]-P[j])) then
if (t[Position(P,P[i]-P[j] )] = 1) then
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]-P[j]) );
else
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]-K[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[i]-P[j] )] = 1) then
Add(KK3, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, P[i]-P[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, K[i]-K[j]) );
else
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, P[i]-P[j]) );
Add(KK4, +1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, K[i]-K[j]) );
fi; 
	fi; 
fi; 
fi; fi;

#############

if (P[i]+K[j] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][Position(P,K[j])][2])=true) then
if (p+Position(P,K[i]+P[j]) > p+Position(P,P[i]+K[j])) then

if (t[Position(P,P[i]+K[j] )] = 1) then
Add(KK3, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]+P[j]) );
Add(KK4, -1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+K[j]) );
else
Add(KK3, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]+P[j]) );
Add(KK4, -1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+K[j]) );
fi; 
else
	if (p+Position(P,K[i]+P[j]) = p+Position(P,P[i]+K[j])) then
if (t[Position(P,K[j]+P[i] )] = 1) then
Add(KK4, -1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+K[j]) );
else
Add(KK3, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]+P[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[i]+K[j] )] = 1) then
Add(KK3, -1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, P[i]+K[j]) );
Add(KK4, -1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, K[i]+P[j]) );
else
Add(KK3, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, P[i]+K[j]) );
Add(KK4, +1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, K[i]+P[j]) );
fi; 
	fi; 
fi; 
fi; fi;

#############

if (P[i]-K[j] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][p+Position(P,K[j])][2])=true) then
if (p+Position(P,K[i]-P[j]) > p+Position(P,P[i]-K[j])) then

if (t[Position(P,P[i]-K[j] )] = 1) then
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]-P[j]) );
Add(KK4, -1*(-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]-K[j]) );
else
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]-P[j]) );
Add(KK4, -1*(-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]-K[j]) );
fi; 
else
	if (p+Position(P,K[i]-P[j]) = p+Position(P,P[i]-K[j])) then
if (t[Position(P,P[i]-K[j] )] = 1) then
Add(KK4, -1*(-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]-K[j]) );
else
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]-P[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[i]-K[j] )] = 1) then
Add(KK3, -1*(+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, P[i]-K[j]) );
Add(KK4, -1*(-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, K[i]-P[j]) );
else
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, P[i]-K[j]) );
Add(KK4, +1*(-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, K[i]-P[j]) );
fi; 
	fi; 
fi; 
fi; fi;

#####################
else #t[j]=-1
#########################

if (P[i]+P[j] in P ) then #qui so che aplha+beta < phi_aplha+phi_beta quindi sono nel caso sicuro
if not (IsEmpty(To[Position(P, P[i])][Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]+K[i]) > p+Position(P,P[j]+P[i])) then

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK3, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]+K[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+P[j]) );
else
Add(KK3, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]+K[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+P[j]) );
fi; 
else
	if (p+Position(P,K[j]+K[i]) = p+Position(P,P[j]+P[i])) then
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+P[j]) );
else
Add(KK3, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]+K[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK3, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, P[i]+P[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P,K[j]+K[i]) );
else
Add(KK3, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, P[i]+P[j]) );
Add(KK4, +1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, K[j]+K[i]) );
fi; 
	fi; 
fi; 
fi; fi;

###############

if (P[i]-P[j] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][p+Position(P,P[j])][2])=true) then
if (p+Position(P,K[i]-K[j]) > p+Position(P,P[i]-P[j])) then

if (t[Position(P,P[i]-P[j] )] = 1) then
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]-K[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]-P[j]) );
else
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]-K[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]-P[j]) );
fi; 
else
	if (p+Position(P,K[i]-K[j]) = p+Position(P,P[i]-P[j])) then
if (t[Position(P,P[i]-P[j] )] = 1) then
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]-P[j]) );
else
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]-K[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[i]-P[j] )] = 1) then
Add(KK3, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, P[i]-P[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, K[i]-K[j]) );
else
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, P[i]-P[j]) );
Add(KK4, +1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, K[i]-K[j]) );
fi; 
	fi; 
fi; 
fi; fi;

#############

if (P[i]+K[j] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][Position(P,K[j])][2])=true) then
if (p+Position(P,K[i]+P[j]) > p+Position(P,P[i]+K[j])) then

if (t[Position(P,P[i]+K[j] )] = 1) then
Add(KK3, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]+P[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+K[j]) );
else
Add(KK3, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]+P[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+K[j]) );
fi; 
else
	if (p+Position(P,K[i]+P[j]) = p+Position(P,P[i]+K[j])) then
if (t[Position(P,K[j]+P[i] )] = 1) then
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+K[j]) );
else
Add(KK3, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]+P[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[i]+K[j] )] = 1) then
Add(KK3, -1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, P[i]+K[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, K[i]+P[j]) );
else
Add(KK3, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, P[i]+K[j]) );
Add(KK4, +1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, K[i]+P[j]) );
fi; 
	fi; 
fi; 
fi; fi;

#############

if (P[i]-K[j] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][p+Position(P,K[j])][2])=true) then
if (p+Position(P,K[i]-P[j]) > p+Position(P,P[i]-K[j])) then

if (t[Position(P,P[i]-K[j] )] = 1) then
Add(KK3, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]-P[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]-K[j]) );
else
Add(KK3, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]-P[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]-K[j]) );
fi; 
else
	if (p+Position(P,K[i]-P[j]) = p+Position(P,P[i]-K[j])) then
if (t[Position(P,P[i]-K[j] )] = 1) then
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]-K[j]) );
else
Add(KK3, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]-P[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[i]-K[j] )] = 1) then
Add(KK3, -1*(-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, P[i]-K[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, K[i]-P[j]) );
else
Add(KK3, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, P[i]-K[j]) );
Add(KK4, +1*(+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, K[i]-P[j]) );
fi; 
	fi; 
fi; 
fi; fi;

fi; fi;

if ( Position(P,K[i]) > Position(P,K[j]) ) then
SetEntrySCTable( TT, Position(P,K[i]),  p+Position(P,P[j]), KK3);
SetEntrySCTable( TT, Position(P,K[i]),  p+Position(P, K[j]), KK4);
else
	if ( Position(P,K[i]) < Position(P,K[j]) ) then
SetEntrySCTable( TT, Position(P,K[i]),  p+Position(P,P[j]), KK3);
	fi;
fi;


else

	if ( p+Position(P, K[j]) = p+Position(P,P[j]) ) then #beta = phi_beta

KK1:=[];
KK2:=[];

if (t[i]=1) then
if (t[j]=1) then

if (P[i]+P[j] in P ) then #qui so che aplha+beta < phi_aplha+phi_beta quindi sono nel caso sicuro
if not (IsEmpty(To[Position(P, P[i])][Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]+K[i]) > p+Position(P,P[j]+P[i])) then

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[i]+K[j]) );
else
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[i]+K[j]) );
fi; 
else
	if (p+Position(P,K[j]+K[i]) = p+Position(P,P[j]+P[i])) then
if (t[Position(P,P[j]+P[i] )] = 1) then
else
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[i]+K[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[i]+P[j]) );
else
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[i]+P[j]) );
fi; 
	fi; 
fi; 
fi; fi;

###############

if (P[i]-P[j] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][p+Position(P,P[j])][2])=true) then
if (p+Position(P,K[i]-K[j]) > p+Position(P,P[i]-P[j])) then

if (t[Position(P,P[i]-P[j] )] = 1) then
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[i]-K[j]) );
else
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[i]-K[j]) );
fi; 
else
	if (p+Position(P,K[i]-K[j]) = p+Position(P,P[i]-P[j])) then
if (t[Position(P,P[i]-P[j] )] = 1) then
else
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[i]-K[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[i]-P[j] )] = 1) then
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[i]-P[j]) );
else
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[i]-P[j]) );
fi; 
	fi; 
fi; 
fi; fi;

####################
else #t[j]=-1

if (P[i]+P[j] in P ) then #qui so che aplha+beta < phi_aplha+phi_beta quindi sono nel caso sicuro
if not (IsEmpty(To[Position(P, P[i])][Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]+K[i]) > p+Position(P,P[j]+P[i])) then

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[i]+P[j]) );
else
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[i]+P[j]) );
fi; 
else
	if (p+Position(P,K[j]+K[i]) = p+Position(P,P[j]+P[i])) then
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[i]+P[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[i]+P[j]) );
else
Add(KK2, +1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[j]+K[i]) );
fi; 
	fi; 
fi; 
fi; fi;

###############

if (P[i]-P[j] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][p+Position(P,P[j])][2])=true) then
if (p+Position(P,K[i]-K[j]) > p+Position(P,P[i]-P[j])) then

if (t[Position(P,P[i]-P[j] )] = 1) then
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[i]-P[j]) );
else
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[i]-P[j]) );
fi; 
else
	if (p+Position(P,K[i]-K[j]) = p+Position(P,P[i]-P[j])) then
if (t[Position(P,P[i]-P[j] )] = 1) then
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[i]-P[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[i]-P[j] )] = 1) then
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[i]-K[j]) );
else
Add(KK2, +1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[i]-K[j]) );
fi; 
	fi; 
fi; 
fi; fi;

fi;

else #t[i]=-1
if (t[j]=1) then

if (P[i]+P[j] in P ) then #qui so che aplha+beta < phi_aplha+phi_beta quindi sono nel caso sicuro
if not (IsEmpty(To[Position(P, P[i])][Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]+K[i]) > p+Position(P,P[j]+P[i])) then

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[i]+K[j]) );
else
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[i]+K[j]) );
fi; 
else
	if (p+Position(P,K[j]+K[i]) = p+Position(P,P[j]+P[i])) then
if (t[Position(P,P[j]+P[i] )] = 1) then
else
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[i]+K[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[i]+P[j]) );
else
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[i]+P[j]) );
fi; 
	fi; 
fi; 
fi; fi;

###############

if (P[i]-P[j] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][p+Position(P,P[j])][2])=true) then
if (p+Position(P,K[i]-K[j]) > p+Position(P,P[i]-P[j])) then

if (t[Position(P,P[i]-P[j] )] = 1) then
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[i]-K[j]) );
else
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[i]-K[j]) );
fi; 
else
	if (p+Position(P,K[i]-K[j]) = p+Position(P,P[i]-P[j])) then
if (t[Position(P,P[i]-P[j] )] = 1) then
else
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[i]-K[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[i]-P[j] )] = 1) then
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[i]-P[j]) );
else
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[i]-P[j]) );
fi; 
	fi; 
fi; 
fi; fi;

####################
else #t[i]=t[j]=-1

if (P[i]+P[j] in P ) then #qui so che aplha+beta < phi_aplha+phi_beta quindi sono nel caso sicuro
if not (IsEmpty(To[Position(P, P[i])][Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]+K[i]) > p+Position(P,P[j]+P[i])) then

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[j]+P[i]) );
else
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[j]+P[i]) );
fi; 
else
	if (p+Position(P,K[j]+K[i]) = p+Position(P,P[j]+P[i])) then
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[j]+P[i]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[i]+K[j]) );
else
Add(KK2, +1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[i]+K[j]) );
fi; 
	fi; 
fi; 
fi; fi;

###############

if (P[i]-P[j] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][p+Position(P,P[j])][2])=true) then
if (p+Position(P,K[i]-K[j]) > p+Position(P,P[i]-P[j])) then

if (t[Position(P,P[i]-P[j] )] = 1) then
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[i]-P[j]) );
else
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[i]-P[j]) );
fi; 
else
	if (p+Position(P,K[i]-K[j]) = p+Position(P,P[i]-P[j])) then
if (t[Position(P,P[i]-P[j] )] = 1) then
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[i]-P[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[i]-P[j] )] = 1) then
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[i]-K[j]) );
else
Add(KK2, +1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[i]-K[j]) );
fi; 
	fi; 
fi; 
fi; fi;

fi; fi;

SetEntrySCTable( TT, Position(P,K[i]), p+j, KK2);

fi; fi; fi; 
od; od;

#########################################################################
#####################################################################################

#############################################################################################################################################
###################################################################################################################

#Prodotti [X,Y] con indici differenti X<Y ovvero alpha < beta !!!!!!!


for i in [1..p-1] do #variazione alpha per gli X
	for j in [i+1..p] do #variazione beta >alpha per gli Y
		if ( Position(P, K[i]) > Position(P,P[i]) ) then #alpha< phi_alpha (ragiono con gli X)
			if ( p+Position(P, K[j]) > p+Position(P,P[j]) ) then # beta< phi_beta (ragiono con gli Y)

KK1:=[]; # i p+j
KK2:=[]; # i p+jbar
KK3:=[]; # ibar p+j
KK4:=[]; # ibar p+jbar

if (t[i]=1) then
if (t[j]=1) then

if (P[i]+P[j] in P ) then #qui so che aplha+beta < phi_aplha+phi_beta quindi sono nel caso sicuro
if not (IsEmpty(To[Position(P, P[i])][Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]+K[i]) > p+Position(P,P[j]+P[i])) then

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, p+Position(P, P[i]+P[j]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]+K[j]) );
Add(KK3, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]+K[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+P[j]) );
else
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, p+Position(P, P[i]+P[j]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]+K[j]) );
Add(KK3, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]+K[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+P[j]) );
fi; 
else
	if (p+Position(P,K[j]+K[i]) = p+Position(P,P[j]+P[i])) then

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, p+Position(P, P[i]+P[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+P[j]) );
else
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]+K[j]) );
Add(KK3, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]+K[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, p+Position(P, K[j]+K[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, P[i]+P[j]) );
Add(KK3, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, P[i]+P[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P,K[j]+K[i]) );
else
Add(KK1, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, p+Position(P, K[j]+K[i]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, P[i]+P[j]) );
Add(KK3, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, P[i]+P[j]) );
Add(KK4, +1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, K[j]+K[i]) );
fi; 
	fi; 
fi; 
fi; fi;

###############

if (P[j]-P[i] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][p+Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]-K[i]) > p+Position(P,P[j]-P[i])) then

if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK1, p+Position(P, P[j]-P[i]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, K[j]-K[i]) );
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, K[j]-K[i]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[j]-P[i]) );
else
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK1, p+Position(P, P[j]-P[i]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, K[j]-K[i]) );
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, K[j]-K[i]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[j]-P[i]) );
fi; 

else
	if (p+Position(P,K[j]-K[i]) = p+Position(P,P[j]-P[i])) then

if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK1, p+Position(P, P[j]-P[i]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[j]-P[i]) );
else
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, K[j]-K[i]) );
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, K[j]-K[i]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))

if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK1, p+Position(P, K[j]-K[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, P[j]-P[i]) );
Add(KK3, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, P[j]-P[i]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, K[j]-K[i]) );
else
Add(KK1, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK1, p+Position(P, K[j]-K[i]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, P[j]-P[i]) );
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, P[j]-P[i]) );
Add(KK4, +1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, K[j]-K[i]) );
fi; 
	fi; 
fi; 
fi; fi;

#############

if (P[i]+K[j] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][Position(P,K[j])][2])=true) then
if (p+Position(P,K[i]+P[j]) > p+Position(P,P[i]+K[j])) then

if (t[Position(P,P[i]+K[j] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, p+Position(P, P[i]+K[j]) );
Add(KK2, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]+P[j]) );
Add(KK3, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]+P[j]) );
Add(KK4, -1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+K[j]) );
else
Add(KK1, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, p+Position(P, P[i]+K[j]) );
Add(KK2, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]+P[j]) );
Add(KK3, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]+P[j]) );
Add(KK4, -1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+K[j]) );
fi; 

else
	if (p+Position(P,K[i]+P[j]) = p+Position(P,P[i]+K[j])) then

if (t[Position(P,K[j]+P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, p+Position(P, P[i]+K[j]) );
Add(KK4, -1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+K[j]) );
else
Add(KK2, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]+P[j]) );
Add(KK3, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]+P[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))

if (t[Position(P,P[i]+K[j] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, p+Position(P, K[i]+P[j]) );
Add(KK2, -1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, P[i]+K[j]) );
Add(KK3, -1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, P[i]+K[j]) );
Add(KK4, -1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, K[i]+P[j]) );
else
Add(KK1, -1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, p+Position(P, K[i]+P[j]) );
Add(KK2, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, P[i]+K[j]) );
Add(KK3, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, P[i]+K[j]) );
Add(KK4, +1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, K[i]+P[j]) );
fi; 
	fi; 
fi; 
fi; fi;

##########

if (K[j]-P[i] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][p+Position(P,K[j])][2])=true) then
if (p+Position(P,P[j]-K[i]) > p+Position(P,K[j]-P[i])) then

if (t[Position(P,K[j]-P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK1, p+Position(P, K[j]-P[i]) );
Add(KK2, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, P[j]-K[i]) );
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, P[j]-K[i]) );
Add(KK4, -1*(-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, K[j]-P[i]) );
else
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK1, p+Position(P, K[j]-P[i]) );
Add(KK2, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, P[j]-K[i]) );
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, P[j]-K[i]) );
Add(KK4, -1*(-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, K[j]-P[i]) );
fi; 
else
	if (p+Position(P,P[j]-K[i]) = p+Position(P,K[j]-P[i])) then

if (t[Position(P,P[j]-K[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK1, p+Position(P, K[j]-P[i]) );
Add(KK4, -1*(-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, K[j]-P[i]) );
else
Add(KK2, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, P[j]-K[i]) );
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, P[j]-K[i]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))

if (t[Position(P,K[j]-P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK1, p+Position(P, P[j]-K[i]) );
Add(KK2, -1*(-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, K[j]-P[i]) );
Add(KK3, -1*(+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, K[j]-P[i]) );
Add(KK4, -1*(-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[j]-K[i]) );
else
Add(KK1, -1*(+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK1, p+Position(P, P[j]-K[i]) );
Add(KK2, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, K[j]-P[i]) );
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, K[j]-P[i]) );
Add(KK4, +1*(-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[j]-K[i]) );
fi; 
	fi; 
fi; 
fi; fi;

##################
else # ovvero t[j]=-1
#########################

if (P[i]+P[j] in P ) then #qui so che aplha+beta < phi_aplha+phi_beta quindi sono nel caso sicuro
if not (IsEmpty(To[Position(P, P[i])][Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]+K[i]) > p+Position(P,P[j]+P[i])) then

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, p+Position(P, P[i]+P[j]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]+K[j]) );
Add(KK3, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]+K[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+P[j]) );
else
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, p+Position(P, P[i]+P[j]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]+K[j]) );
Add(KK3, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]+K[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+P[j]) );
fi; 
else
	if (p+Position(P,K[j]+K[i]) = p+Position(P,P[j]+P[i])) then

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, p+Position(P, P[i]+P[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+P[j]) );
else
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]+K[j]) );
Add(KK3, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]+K[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, p+Position(P, K[j]+K[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, P[i]+P[j]) );
Add(KK3, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, P[i]+P[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P,K[j]+K[i]) );
else
Add(KK1, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, p+Position(P, K[j]+K[i]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, P[i]+P[j]) );
Add(KK3, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, P[i]+P[j]) );
Add(KK4, +1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, K[j]+K[i]) );
fi; 
	fi; 
fi; 
fi; fi;

###############

if (P[j]-P[i] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][p+Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]-K[i]) > p+Position(P,P[j]-P[i])) then

if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK1, p+Position(P, P[j]-P[i]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, K[j]-K[i]) );
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, K[j]-K[i]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[j]-P[i]) );
else
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK1, p+Position(P, P[j]-P[i]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, K[j]-K[i]) );
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, K[j]-K[i]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[j]-P[i]) );
fi; 
else
	if (p+Position(P,K[j]-K[i]) = p+Position(P,P[j]-P[i])) then

if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK1, p+Position(P, P[j]-P[i]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[j]-P[i]) );
else
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, K[j]-K[i]) );
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, K[j]-K[i]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))

if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK1, p+Position(P, K[j]-K[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, P[j]-P[i]) );
Add(KK3, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, P[j]-P[i]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, K[j]-K[i]) );
else
Add(KK1, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK1, p+Position(P, K[j]-K[i]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, P[j]-P[i]) );
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, P[j]-P[i]) );
Add(KK4, +1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, K[j]-K[i]) );
fi; 
	fi; 
fi; 
fi; fi;

#############

if (P[i]+K[j] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][Position(P,K[j])][2])=true) then
if (p+Position(P,K[i]+P[j]) > p+Position(P,P[i]+K[j])) then

if (t[Position(P,P[i]+K[j] )] = 1) then
Add(KK1, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, p+Position(P, P[i]+K[j]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]+P[j]) );
Add(KK3, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]+P[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+K[j]) );
else
Add(KK1, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, p+Position(P, P[i]+K[j]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]+P[j]) );
Add(KK3, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]+P[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+K[j]) );
fi; 
else
	if (p+Position(P,K[i]+P[j]) = p+Position(P,P[i]+K[j])) then

if (t[Position(P,K[j]+P[i] )] = 1) then
Add(KK1, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, p+Position(P, P[i]+K[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+K[j]) );
else
Add(KK2, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]+P[j]) );
Add(KK3, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]+P[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))

if (t[Position(P,P[i]+K[j] )] = 1) then
Add(KK1, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, p+Position(P, K[i]+P[j]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, P[i]+K[j]) );
Add(KK3, -1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, P[i]+K[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, K[i]+P[j]) );
else
Add(KK1, -1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, p+Position(P, K[i]+P[j]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, P[i]+K[j]) );
Add(KK3, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, P[i]+K[j]) );
Add(KK4, +1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, K[i]+P[j]) );
fi; 
	fi; 
fi; 
fi; fi;

#############

if (K[j]-P[i] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][p+Position(P,K[j])][2])=true) then
if (p+Position(P,P[j]-K[i]) > p+Position(P,K[j]-P[i])) then

if (t[Position(P,K[j]-P[i] )] = 1) then
Add(KK1, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK1, p+Position(P, K[j]-P[i]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, P[j]-K[i]) );
Add(KK3, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, P[j]-K[i]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, K[j]-P[i]) );
else
Add(KK1, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK1, p+Position(P, K[j]-P[i]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, P[j]-K[i]) );
Add(KK3, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, P[j]-K[i]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, K[j]-P[i]) );
fi; 
else
	if (p+Position(P,P[j]-K[i]) = p+Position(P,K[j]-P[i])) then

if (t[Position(P,K[j]-P[i] )] = 1) then
Add(KK1, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK1, p+Position(P, K[j]-P[i]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, K[j]-P[i]) );
else
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, P[j]-K[i]) );
Add(KK3, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, P[j]-K[i]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))

if (t[Position(P,K[j]-P[i] )] = 1) then
Add(KK1, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK1, p+Position(P, P[j]-K[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, K[j]-P[i]) );
Add(KK3, -1*(-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, K[j]-P[i]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[j]-K[i]) );
else
Add(KK1, -1*(-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK1, p+Position(P, P[j]-K[i]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, K[j]-P[i]) );
Add(KK3, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, K[j]-P[i]) );
Add(KK4, +1*(+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[j]-K[i]) );
fi; 
	fi; 
fi; 
fi; fi;

fi;

########################à
else #t[i]=-1
if (t[j]=1) then
###############

if (P[i]+P[j] in P ) then #qui so che aplha+beta < phi_aplha+phi_beta quindi sono nel caso sicuro
if not (IsEmpty(To[Position(P, P[i])][Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]+K[i]) > p+Position(P,P[j]+P[i])) then
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, p+Position(P, P[i]+P[j]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]+K[j]) );
Add(KK3, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]+K[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+P[j]) );
else
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, p+Position(P, P[i]+P[j]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]+K[j]) );
Add(KK3, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]+K[j]) );
Add(KK4, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+P[j]) );
fi; 
else
	if (p+Position(P,K[j]+K[i]) = p+Position(P,P[j]+P[i])) then
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, p+Position(P, P[i]+P[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+P[j]) );
else
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]+K[j]) );
Add(KK3, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]+K[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, p+Position(P, K[j]+K[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, P[i]+P[j]) );
Add(KK3, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, P[i]+P[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P,K[j]+K[i]) );
else
Add(KK1, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, p+Position(P, K[j]+K[i]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, P[i]+P[j]) );
Add(KK3, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, P[i]+P[j]) );
Add(KK4, +1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, K[j]+K[i]) );
fi; 
	fi; 
fi; 
fi; fi;

###############

if (P[j]-P[i] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][p+Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]-K[i]) > p+Position(P,P[j]-P[i])) then

if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK1, p+Position(P, P[j]-P[i]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, K[j]-K[i]) );
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, K[j]-K[i]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[j]-P[i]) );
else
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK1, p+Position(P, P[j]-P[i]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, K[j]-K[i]) );
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, K[j]-K[i]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[j]-P[i]) );
fi; 
else
	if (p+Position(P,K[j]-K[i]) = p+Position(P,P[j]-P[i])) then
if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK1, p+Position(P, P[j]-P[i]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[j]-P[i]) );
else
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, K[j]-K[i]) );
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, K[j]-K[i]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK1, p+Position(P, K[j]-K[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, P[j]-P[i]) );
Add(KK3, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, P[j]-P[i]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, K[j]-K[i]) );
else
Add(KK1, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK1, p+Position(P, K[j]-K[i]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, P[j]-P[i]) );
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, P[j]-P[i]) );
Add(KK4, +1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, K[j]-K[i]) );
fi; 
	fi; 
fi; 
fi; fi;

#############

if (P[i]+K[j] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][Position(P,K[j])][2])=true) then
if (p+Position(P,K[i]+P[j]) > p+Position(P,P[i]+K[j])) then
if (t[Position(P,P[i]+K[j] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, p+Position(P, P[i]+K[j]) );
Add(KK2, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]+P[j]) );
Add(KK3, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]+P[j]) );
Add(KK4, -1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+K[j]) );
else
Add(KK1, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, p+Position(P, P[i]+K[j]) );
Add(KK2, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]+P[j]) );
Add(KK3, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]+P[j]) );
Add(KK4, -1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+K[j]) );
fi; 
else
	if (p+Position(P,K[i]+P[j]) = p+Position(P,P[i]+K[j])) then
if (t[Position(P,K[j]+P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, p+Position(P, P[i]+K[j]) );
Add(KK4, -1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+K[j]) );
else
Add(KK2, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]+P[j]) );
Add(KK3, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]+P[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[i]+K[j] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, p+Position(P, K[i]+P[j]) );
Add(KK2, -1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, P[i]+K[j]) );
Add(KK3, -1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, P[i]+K[j]) );
Add(KK4, -1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, K[i]+P[j]) );
else
Add(KK1, -1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, p+Position(P, K[i]+P[j]) );
Add(KK2, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, P[i]+K[j]) );
Add(KK3, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, P[i]+K[j]) );
Add(KK4, +1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, K[i]+P[j]) );
fi; 
	fi; 
fi; 
fi; fi;

#############

if (K[j]-P[i] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][p+Position(P,K[j])][2])=true) then
if (p+Position(P,P[j]-K[i]) > p+Position(P,K[j]-P[i])) then

if (t[Position(P,K[j]-P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK1, p+Position(P, K[j]-P[i]) );
Add(KK2, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, P[j]-K[i]) );
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, P[j]-K[i]) );
Add(KK4, -1*(-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, K[j]-P[i]) );
else
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK1, p+Position(P, K[j]-P[i]) );
Add(KK2, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, P[j]-K[i]) );
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, P[j]-K[i]) );
Add(KK4, -1*(-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, K[j]-P[i]) );
fi; 
else
	if (p+Position(P,P[j]-K[i]) = p+Position(P,K[j]-P[i])) then
if (t[Position(P,K[j]-P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK1, p+Position(P, K[j]-P[i]) );
Add(KK4, -1*(-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, K[j]-P[i]) );
else
Add(KK2, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, P[j]-K[i]) );
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, P[j]-K[i]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,K[j]-P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK1, p+Position(P, P[j]-K[i]) );
Add(KK2, -1*(-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, K[j]-P[i]) );
Add(KK3, -1*(+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, K[j]-P[i]) );
Add(KK4, -1*(-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[j]-K[i]) );
else
Add(KK1, -1*(+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK1, p+Position(P, P[j]-K[i]) );
Add(KK2, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, K[j]-P[i]) );
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, K[j]-P[i]) );
Add(KK4, +1*(-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[j]-K[i]) );
fi; 
	fi; 
fi; 
fi; fi;

else #t[j]=-1

###############

if (P[i]+P[j] in P ) then #qui so che aplha+beta < phi_aplha+phi_beta quindi sono nel caso sicuro
if not (IsEmpty(To[Position(P, P[i])][Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]+K[i]) > p+Position(P,P[j]+P[i])) then

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, p+Position(P, P[i]+P[j]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]+K[j]) );
Add(KK3, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]+K[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+P[j]) );
else
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, p+Position(P, P[i]+P[j]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]+K[j]) );
Add(KK3, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]+K[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+P[j]) );
fi; 
else
	if (p+Position(P,K[j]+K[i]) = p+Position(P,P[j]+P[i])) then
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, p+Position(P, P[i]+P[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+P[j]) );
else
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]+K[j]) );
Add(KK3, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]+K[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, p+Position(P, K[j]+K[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, P[i]+P[j]) );
Add(KK3, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, P[i]+P[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P,K[j]+K[i]) );
else
Add(KK1, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK1, p+Position(P, K[j]+K[i]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, P[i]+P[j]) );
Add(KK3, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, P[i]+P[j]) );
Add(KK4, +1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, K[j]+K[i]) );
fi; 
	fi; 
fi; 
fi; fi;

#######

if (P[j]-P[i] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][p+Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]-K[i]) > p+Position(P,P[j]-P[i])) then

if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK1, p+Position(P, P[j]-P[i]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, K[j]-K[i]) );
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, K[j]-K[i]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[j]-P[i]) );
else
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK1, p+Position(P, P[j]-P[i]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, K[j]-K[i]) );
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, K[j]-K[i]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[j]-P[i]) );
fi; 
else
	if (p+Position(P,K[j]-K[i]) = p+Position(P,P[j]-P[i])) then
if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK1, p+Position(P, P[j]-P[i]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[j]-P[i]) );
else
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, K[j]-K[i]) );
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, K[j]-K[i]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK1, p+Position(P, K[j]-K[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, P[j]-P[i]) );
Add(KK3, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, P[j]-P[i]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, K[j]-K[i]) );
else
Add(KK1, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK1, p+Position(P, K[j]-K[i]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, P[j]-P[i]) );
Add(KK3, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK3, p+Position(P, P[j]-P[i]) );
Add(KK4, +1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, K[j]-K[i]) );
fi; 
	fi; 
fi; 
fi; fi;

######

if (P[i]+K[j] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][Position(P,K[j])][2])=true) then
if (p+Position(P,K[i]+P[j]) > p+Position(P,P[i]+K[j])) then

if (t[Position(P,P[i]+K[j] )] = 1) then
Add(KK1, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, p+Position(P, P[i]+K[j]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]+P[j]) );
Add(KK3, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]+P[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+K[j]) );
else
Add(KK1, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, p+Position(P, P[i]+K[j]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]+P[j]) );
Add(KK3, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]+P[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+K[j]) );
fi; 
else
	if (p+Position(P,K[i]+P[j]) = p+Position(P,P[i]+K[j])) then
if (t[Position(P,K[j]+P[i] )] = 1) then
Add(KK1, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, p+Position(P, P[i]+K[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+K[j]) );
else
Add(KK2, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]+P[j]) );
Add(KK3, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, K[i]+P[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[i]+K[j] )] = 1) then
Add(KK1, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, p+Position(P, K[i]+P[j]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, P[i]+K[j]) );
Add(KK3, -1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, P[i]+K[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, K[i]+P[j]) );
else
Add(KK1, -1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK1, p+Position(P, K[i]+P[j]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, P[i]+K[j]) );
Add(KK3, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, P[i]+K[j]) );
Add(KK4, +1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, K[i]+P[j]) );
fi; 
	fi; 
fi; 
fi; fi;

#############

if (K[j]-P[i] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][p+Position(P,K[j])][2])=true) then
if (p+Position(P,P[j]-K[i]) > p+Position(P,K[j]-P[i])) then

if (t[Position(P,K[j]-P[i] )] = 1) then
Add(KK1, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK1, p+Position(P, K[j]-P[i]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, P[j]-K[i]) );
Add(KK3, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, P[j]-K[i]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, K[j]-P[i]) );
else
Add(KK1, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK1, p+Position(P, K[j]-P[i]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, P[j]-K[i]) );
Add(KK3, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, P[j]-K[i]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, K[j]-P[i]) );
fi; 
else
	if (p+Position(P,P[j]-K[i]) = p+Position(P,K[j]-P[i])) then
if (t[Position(P,K[j]-P[i] )] = 1) then
Add(KK1, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK1, p+Position(P, K[j]-P[i]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, K[j]-P[i]) );
else
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, P[j]-K[i]) );
Add(KK3, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, P[j]-K[i]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,K[j]-P[i] )] = 1) then
Add(KK1, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK1, p+Position(P, P[j]-K[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, K[j]-P[i]) );
Add(KK3, -1*(-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, K[j]-P[i]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[j]-K[i]) );
else
Add(KK1, -1*(-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK1, p+Position(P, P[j]-K[i]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, K[j]-P[i]) );
Add(KK3, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK3, p+Position(P, K[j]-P[i]) );
Add(KK4, +1*(+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[j]-K[i]) );
fi; 
	fi; 
fi; 
fi; fi;

fi; fi;

if ( Position(P,K[i]) < Position(P,P[j]) ) then
SetEntrySCTable( TT, Position(P,P[i]),  p+Position(P,P[j]) , KK1);
SetEntrySCTable( TT, Position(P,P[i]),  p+Position(P, K[j]), KK2);
SetEntrySCTable( TT, Position(P,K[i]),  p+Position(P,P[j]), KK3);
SetEntrySCTable( TT, Position(P,K[i]),  p+Position(P, K[j]), KK4);
fi;

if ( Position(P,K[i]) > Position(P,K[j]) ) then
SetEntrySCTable( TT, Position(P,P[i]),  p+Position(P,P[j]) , KK1);
SetEntrySCTable( TT, Position(P,P[i]),  p+Position(P, K[j]), KK2);
fi;

if ( Position(P,P[j]) < Position(P,K[i]) ) then
	if ( Position(P,K[i]) < Position(P,K[j]) ) then
SetEntrySCTable( TT, Position(P,P[i]),  p+Position(P,P[j]) , KK1);
SetEntrySCTable( TT, Position(P,P[i]),  p+Position(P, K[j]), KK2);
#SetEntrySCTable( TT, Position(P,K[i]),  p+Position(P,P[j]), KK3);
SetEntrySCTable( TT, Position(P,K[i]),  p+Position(P, K[j]), KK4);
	fi;
fi;

else
	if ( p+Position(P, K[j]) = p+Position(P,P[j]) ) then #beta = phi_beta

KK1:=[];
KK2:=[];

if (t[i]=1) then
if (t[j]=1) then

if (P[i]+P[j] in P ) then #qui so che aplha+beta < phi_aplha+phi_beta quindi sono nel caso sicuro
if not (IsEmpty(To[Position(P, P[i])][Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]+K[i]) > p+Position(P,P[j]+P[i])) then

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, P[i]+P[j]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[i]+K[j]) );
else
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, P[i]+P[j]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[i]+K[j]) );
fi; 
else
	if (p+Position(P,K[j]+K[i]) = p+Position(P,P[j]+P[i])) then
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, P[i]+P[j]) );
else
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[i]+K[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, K[j]+K[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[i]+P[j]) );
else
Add(KK1, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, K[j]+K[i]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[i]+P[j]) );
fi; 
	fi; 
fi; 
fi; fi;

###############

if (P[j]-P[i] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][p+Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]-K[i]) > p+Position(P,P[j]-P[i])) then

if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, P[j]-P[i]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[j]-K[i]) );
else
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, P[j]-P[i]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[j]-K[i]) );
fi; 
else
	if (p+Position(P,K[j]-K[i]) = p+Position(P,P[j]-P[i])) then

if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, P[j]-P[i]) );
else
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[j]-K[i]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))

if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, K[j]-K[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[j]-P[i]) );
else
Add(KK1, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, K[j]-K[i]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[j]-P[i]) );
fi; 
	fi; 
fi; 
fi; fi;

####################
else #t[j]=-1

if (P[i]+P[j] in P ) then #qui so che aplha+beta < phi_aplha+phi_beta quindi sono nel caso sicuro
if not (IsEmpty(To[Position(P, P[i])][Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]+K[i]) > p+Position(P,P[j]+P[i])) then

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, K[j]+K[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[i]+P[j]) );
else
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, K[j]+K[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[i]+P[j]) );
fi; 
else
	if (p+Position(P,K[j]+K[i]) = p+Position(P,P[j]+P[i])) then
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[i]+P[j]) );
else
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, K[j]+K[i]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, K[j]+K[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[i]+P[j]) );
else
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, P[i]+P[j]) );
Add(KK2, +1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[j]+K[i]) );
fi; 
	fi; 
fi; 
fi; fi;

###############

if (P[j]-P[i] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][p+Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]-K[i]) > p+Position(P,P[j]-P[i])) then

if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, K[j]-K[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[j]-P[i]) );
else
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, K[j]-K[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[j]-P[i]) );
fi; 

else
	if (p+Position(P,K[j]-K[i]) = p+Position(P,P[j]-P[i])) then
if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[j]-P[i]) );
else
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, K[j]-K[i]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))

if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, P[j]-P[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[j]-K[i]) );
else
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, P[j]-P[i]) );
Add(KK2, +1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[j]-K[i]) );
fi; 
	fi; 
fi; 
fi; fi;

fi;

else #t[i]=-1
if (t[j]=1) then

if (P[i]+P[j] in P ) then #qui so che aplha+beta < phi_aplha+phi_beta quindi sono nel caso sicuro
if not (IsEmpty(To[Position(P, P[i])][Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]+K[i]) > p+Position(P,P[j]+P[i])) then

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, P[i]+P[j]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[i]+K[j]) );
else
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, P[i]+P[j]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[i]+K[j]) );
fi; 
else
	if (p+Position(P,K[j]+K[i]) = p+Position(P,P[j]+P[i])) then

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, P[i]+P[j]) );
else
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[i]+K[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, K[j]+K[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[i]+P[j]) );
else
Add(KK1, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, K[j]+K[i]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[i]+P[j]) );
fi; 
	fi; 
fi; 
fi; fi;

###############

if (P[j]-P[i] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][p+Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]-K[i]) > p+Position(P,P[j]-P[i])) then

if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, P[j]-P[i]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[j]-K[i]) );
else
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, P[j]-P[i]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[j]-K[i]) );
fi; 
else
	if (p+Position(P,K[j]-K[i]) = p+Position(P,P[j]-P[i])) then

if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, P[j]-P[i]) );
else
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[j]-K[i]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, K[j]-K[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[j]-P[i]) );
else
Add(KK1, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, K[j]-K[i]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[j]-P[i]) );
fi; 
	fi; 
fi; 
fi; fi;

###########
else #t[i]=t[j]=-1

if (P[i]+P[j] in P ) then #qui so che aplha+beta < phi_aplha+phi_beta quindi sono nel caso sicuro
if not (IsEmpty(To[Position(P, P[i])][Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]+K[i]) > p+Position(P,P[j]+P[i])) then

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, K[i]+K[j]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[j]+P[i]) );
else
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, K[i]+K[j]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[j]+P[i]) );
fi; 
else
	if (p+Position(P,K[j]+K[i]) = p+Position(P,P[j]+P[i])) then

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[j]+P[i]) );
else
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, K[i]+K[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, P[j]+P[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[i]+K[j]) );
else
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, P[j]+P[i]) );
Add(KK2, +1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[i]+K[j]) );
fi; 
	fi; 
fi; 
fi; fi;

###############

if (P[j]-P[i] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][p+Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]-K[i]) > p+Position(P,P[j]-P[i])) then

if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, K[j]-K[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[j]-P[i]) );
else
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, K[j]-K[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[j]-P[i]) );
fi; 
else
	if (p+Position(P,K[j]-K[i]) = p+Position(P,P[j]-P[i])) then

if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[j]-P[i]) );
else
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, K[j]-K[i]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, P[j]-P[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[j]-K[i]) );
else
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, P[j]-P[i]) );
Add(KK2, +1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[j]-K[i]) );
fi; 
	fi; 
fi; 
fi; fi;

fi; fi;

if ( Position(P, K[i]) < Position(P, P[j]) ) then
SetEntrySCTable( TT, i,  p+j , KK1);
SetEntrySCTable( TT, Position(P,K[i]), p+j, KK2);
fi;

if ( Position(P, K[i]) > Position(P, P[j]) ) then
SetEntrySCTable( TT, i,  p+j , KK1);
#SetEntrySCTable( TT, Position(P,K[i]), p+j, KK2);
fi;
fi;

fi;

else
	if ( p+Position(P, K[i]) = p+Position(P,P[i]) ) then #alpha< phi_alpha
			if ( p+Position(P, K[j]) > p+Position(P,P[j]) ) then #beta< phi_beta

KK1:=[];
KK2:=[];

if (t[i]=1) then
if (t[j]=1) then

if (P[i]+P[j] in P ) then #qui so che aplha+beta < phi_aplha+phi_beta quindi sono nel caso sicuro
if not (IsEmpty(To[Position(P, P[i])][Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]+K[i]) > p+Position(P,P[j]+P[i])) then

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, P[i]+P[j]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[i]+K[j]) );
else
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, P[i]+P[j]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[i]+K[j]) );
fi; 
else
	if (p+Position(P,K[j]+K[i]) = p+Position(P,P[j]+P[i])) then

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, P[i]+P[j]) );
else
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[i]+K[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, K[j]+K[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[i]+P[j]) );
else
Add(KK1, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, K[j]+K[i]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[i]+P[j]) );
fi; 
	fi; 
fi; 
fi; fi;

###############

if (P[j]-P[i] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][p+Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]-K[i]) > p+Position(P,P[j]-P[i])) then

if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, P[j]-P[i]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[j]-K[i]) );
else
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, P[j]-P[i]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[j]-K[i]) );
fi; 
else
	if (p+Position(P,K[j]-K[i]) = p+Position(P,P[j]-P[i])) then

if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, P[j]-P[i]) );
else
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[j]-K[i]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))

if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, K[j]-K[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[j]-P[i]) );
else
Add(KK1, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, K[j]-K[i]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[j]-P[i]) );
fi; 
	fi; 
fi; 
fi; fi;

###########

else #t[j]=-1

if (P[i]+P[j] in P ) then #qui so che aplha+beta < phi_aplha+phi_beta quindi sono nel caso sicuro
if not (IsEmpty(To[Position(P, P[i])][Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]+K[i]) > p+Position(P,P[j]+P[i])) then

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, P[i]+P[j]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[i]+K[j]) );
else
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, P[i]+P[j]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[i]+K[j]) );
fi; 
else
	if (p+Position(P,K[j]+K[i]) = p+Position(P,P[j]+P[i])) then

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, P[i]+P[j]) );
else
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[i]+K[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, K[j]+K[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[i]+P[j]) );
else
Add(KK1, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, K[j]+K[i]) );
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[i]+P[j]) );
fi; 
	fi; 
fi; 
fi; fi;

#######

if (P[j]-P[i] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][p+Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]-K[i]) > p+Position(P,P[j]-P[i])) then
if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, P[j]-P[i]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[j]-K[i]) );
else
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, P[j]-P[i]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[j]-K[i]) );
fi; 
else
	if (p+Position(P,K[j]-K[i]) = p+Position(P,P[j]-P[i])) then

if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, P[j]-P[i]) );
else
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[j]-K[i]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, K[j]-K[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[j]-P[i]) );
else
Add(KK1, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, K[j]-K[i]) );
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[j]-P[i]) );
fi; 
	fi; 
fi; 
fi; fi;

####################

fi;

else #t[i]=-1
if (t[j]=1) then

if (P[i]+P[j] in P ) then #qui so che aplha+beta < phi_aplha+phi_beta quindi sono nel caso sicuro
if not (IsEmpty(To[Position(P, P[i])][Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]+K[i]) > p+Position(P,P[j]+P[i])) then

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, K[j]+K[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[i]+P[j]) );
else
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, K[j]+K[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[i]+P[j]) );
fi; 
else
	if (p+Position(P,K[j]+K[i]) = p+Position(P,P[j]+P[i])) then
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[i]+P[j]) );
else
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, K[j]+K[i]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, K[j]+K[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[i]+P[j]) );
else
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, P[i]+P[j]) );
Add(KK2, +1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[j]+K[i]) );
fi; 
	fi; 
fi; 
fi; fi;

############

if (P[j]-P[i] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][p+Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]-K[i]) > p+Position(P,P[j]-P[i])) then

if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, K[j]-K[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[j]-P[i]) );
else
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, K[j]-K[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[j]-P[i]) );
fi; 
else
	if (p+Position(P,K[j]-K[i]) = p+Position(P,P[j]-P[i])) then

if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[j]-P[i]) );
else
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, K[j]-K[i]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, P[j]-P[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[j]-K[i]) );
else
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, P[j]-P[i]) );
Add(KK2, +1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[j]-K[i]) );
fi; 
	fi; 
fi; 
fi; fi;

else

if (P[i]+P[j] in P ) then #qui so che aplha+beta < phi_aplha+phi_beta quindi sono nel caso sicuro
if not (IsEmpty(To[Position(P, P[i])][Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]+K[i]) > p+Position(P,P[j]+P[i])) then

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, K[i]+K[j]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[j]+P[i]) );
else
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, K[i]+K[j]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[j]+P[i]) );
fi; 
else
	if (p+Position(P,K[j]+K[i]) = p+Position(P,P[j]+P[i])) then
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[j]+P[i]) );
else
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, K[i]+K[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK1, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, P[j]+P[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[i]+K[j]) );
else
Add(KK1, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, P[j]+P[i]) );
Add(KK2, +1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[i]+K[j]) );
fi; 
	fi; 
fi; 
fi; fi;

###############

if (P[j]-P[i] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][p+Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]-K[i]) > p+Position(P,P[j]-P[i])) then

if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, K[j]-K[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[j]-P[i]) );
else
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, K[j]-K[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[j]-P[i]) );
fi; 
else
	if (p+Position(P,K[j]-K[i]) = p+Position(P,P[j]-P[i])) then
if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[j]-P[i]) );
else
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, K[j]-K[i]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK1, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, P[j]-P[i]) );
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[j]-K[i]) );
else
Add(KK1, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK1, p+Position(P, P[j]-P[i]) );
Add(KK2, +1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[j]-K[i]) );
fi; 
	fi; 
fi; 
fi; fi;

fi;
fi;

SetEntrySCTable( TT, i,  p+j , KK1);
SetEntrySCTable( TT, i,  p+Position(P,K[j]), KK2);

###########################################

		else
			if ( p+Position(P, K[j]) = p+Position(P,P[j]) ) then #beta = phi_beta

KK1:=[];

if (t[i]=1) then
if (t[j]=1) then

if (P[i]+P[j] in P ) then #qui so che aplha+beta < phi_aplha + phi_beta quindi sono nel caso sicuro
if not (IsEmpty(To[Position(P, P[i])][Position(P,P[j])][2])=true) then
if (t[Position(P,P[j]+P[i])]=1) then
Add(KK1, +2*To[Position(P, P[i])][Position(P,P[j])][2][1]);
Add(KK1, p+Position(P, P[i]+P[j]) );
fi; fi; fi;

if (P[j]-P[i] in P ) then
if not (IsEmpty(To[Position(P, P[i])][p+Position(P,P[j])][2])=true) then
if (t[Position(P,P[j]-P[i])]=1) then
if (p+Position(P,K[j]-K[i]) > p+Position(P,P[j]-P[i])) then
Add(KK1, +2*To[Position(P, P[i])][p+Position(P,P[j])][2][1]);
Add(KK1, p+Position(P, P[j]-P[i]) );
else 
	if (p+Position(P,K[j]-K[i]) = p+Position(P, P[j]-P[i])) then
Add(KK1, +2*To[Position(P, P[i])][p+Position(P,P[j])][2][1]);
Add(KK1, p+Position(P, P[j]-P[i]) );
	else
Add(KK1, +2*To[Position(P, P[i])][p+Position(P,P[j])][2][1]);
Add(KK1, p+Position(P, K[j]-K[i]) );
	fi;
fi; fi; fi; fi;

SetEntrySCTable( TT, i,  p+j , KK1);

else

if (P[i]+P[j] in P ) then #qui so che aplha+beta < phi_aplha+phi_beta quindi sono nel caso sicuro
if not (IsEmpty(To[Position(P, P[i])][Position(P,P[j])][2])=true) then
if (t[Position(P,P[j]+P[i])]=-1) then
Add(KK1, +2*To[Position(P, P[i])][Position(P,P[j])][2][1]);
Add(KK1, p+Position(P, P[i]+P[j]) );
fi; fi; fi;

if (P[j]-P[i] in P ) then
if not (IsEmpty(To[Position(P, P[i])][p+Position(P,P[j])][2])=true) then
if (t[Position(P,P[j]-P[i])]=-1) then
if (p+Position(P,K[j]-K[i]) > p+Position(P,P[j]-P[i])) then
Add(KK1, +2*To[Position(P, P[i])][p+Position(P,P[j])][2][1]);
Add(KK1, p+Position(P, P[j]-P[i]) );
else 
	if (p+Position(P,K[j]-K[i]) = p+Position(P, P[j]-P[i])) then
Add(KK1, +2*To[Position(P, P[i])][p+Position(P,P[j])][2][1]);
Add(KK1, p+Position(P, P[j]-P[i]) );
	else
Add(KK1, +2*To[Position(P, P[i])][p+Position(P,P[j])][2][1]);
Add(KK1, p+Position(P, K[j]-K[i]) );
	fi;
fi; fi; fi; fi;
SetEntrySCTable( TT, i,  p+j , KK1);
fi;

else
if (t[j]=1) then

if (P[i]+P[j] in P ) then #qui so che aplha+beta < phi_aplha+phi_beta quindi sono nel caso sicuro
if not (IsEmpty(To[Position(P, P[i])][Position(P,P[j])][2])=true) then
if (t[Position(P,P[j]+P[i])]=-1) then
Add(KK1, +2*To[Position(P, P[i])][Position(P,P[j])][2][1]);
Add(KK1, p+Position(P, P[i]+P[j]) );
fi; fi; fi;

if (P[j]-P[i] in P ) then
if not (IsEmpty(To[Position(P, P[i])][p+Position(P,P[j])][2])=true) then
if (t[Position(P,P[j]-P[i])]=-1) then
if (p+Position(P,K[j]-K[i]) > p+Position(P,P[j]-P[i])) then
Add(KK1, +2*To[Position(P, P[i])][p+Position(P,P[j])][2][1]);
Add(KK1, p+Position(P, P[j]-P[i]) );
else 
	if (p+Position(P,K[j]-K[i]) = p+Position(P, P[j]-P[i])) then
Add(KK1, +2*To[Position(P, P[i])][p+Position(P,P[j])][2][1]);
Add(KK1, p+Position(P, P[j]-P[i]) );
	else
Add(KK1, +2*To[Position(P, P[i])][p+Position(P,P[j])][2][1]);
Add(KK1, p+Position(P, K[j]-K[i]) );
	fi;
fi; fi; fi; fi;

SetEntrySCTable( TT, i,  p+j , KK1);

else

if (P[i]+P[j] in P ) then #qui so che aplha+beta < phi_aplha+phi_beta quindi sono nel caso sicuro
if not (IsEmpty(To[Position(P, P[i])][Position(P,P[j])][2])=true) then
if (t[Position(P,P[j]+P[i])]=1) then
Add(KK1, -2*To[Position(P, P[i])][Position(P,P[j])][2][1]);
Add(KK1, p+Position(P, P[i]+P[j]) );
fi; fi; fi;

if (P[j]-P[i] in P ) then
if not (IsEmpty(To[Position(P, P[i])][p+Position(P,P[j])][2])=true) then
if (t[Position(P,P[j]-P[i])]=1) then
if (p+Position(P,K[j]-K[i]) > p+Position(P,P[j]-P[i])) then
Add(KK1, -2*To[Position(P, P[i])][p+Position(P,P[j])][2][1]);
Add(KK1, p+Position(P, P[j]-P[i]) );
else 
	if (p+Position(P,K[j]-K[i]) = p+Position(P, P[j]-P[i])) then
Add(KK1, -2*To[Position(P, P[i])][p+Position(P,P[j])][2][1]);
Add(KK1, p+Position(P, P[j]-P[i]) );
	else
Add(KK1, -2*To[Position(P, P[i])][p+Position(P,P[j])][2][1]);
Add(KK1, p+Position(P, K[j]-K[i]) );
	fi;
fi; fi; fi; fi;

SetEntrySCTable( TT, i,  p+j , KK1);

fi;
fi; fi; fi;
fi; fi;
od; od;

###############
###########à
###########
for i in [2..p] do 
	for j in [1..i-1] do
	if ( Position(P, K[j]) > Position(P, P[i]) ) then 		
		if ( Position(P, K[i]) > Position(P,P[i]) ) then #alpha< phi_alpha (ragiono con gli X)

KK2:=[]; # i p+jbar
KK4:=[]; # ibar p+jbar

if (t[i]=1) then
if (t[j]=1) then

if (P[i]+P[j] in P ) then #qui so che aplha+beta < phi_aplha+phi_beta quindi sono nel caso sicuro
if not (IsEmpty(To[Position(P, P[i])][Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]+K[i]) > p+Position(P,P[j]+P[i])) then

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]+K[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+P[j]) );
else
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]+K[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+P[j]) );
fi; 
else
	if (p+Position(P,K[j]+K[i]) = p+Position(P,P[j]+P[i])) then

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+P[j]) );
else
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]+K[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, P[i]+P[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P,K[j]+K[i]) );
else
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, P[i]+P[j]) );
Add(KK4, +1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, K[j]+K[i]) );
fi; 
	fi; 
fi; 
fi; fi;

###############

if (P[j]-P[i] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][p+Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]-K[i]) > p+Position(P,P[j]-P[i])) then

if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, K[j]-K[i]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[j]-P[i]) );
else
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, K[j]-K[i]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[j]-P[i]) );
fi; 

else
	if (p+Position(P,K[j]-K[i]) = p+Position(P,P[j]-P[i])) then

if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[j]-P[i]) );
else
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, K[j]-K[i]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))

if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, P[j]-P[i]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, K[j]-K[i]) );
else
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, P[j]-P[i]) );
Add(KK4, +1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, K[j]-K[i]) );
fi; 
	fi; 
fi; 
fi; fi;

#############

if (P[i]+K[j] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][Position(P,K[j])][2])=true) then
if (p+Position(P,K[i]+P[j]) > p+Position(P,P[i]+K[j])) then

if (t[Position(P,P[i]+K[j] )] = 1) then
Add(KK2, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]+P[j]) );
Add(KK4, -1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+K[j]) );
else
Add(KK2, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]+P[j]) );
Add(KK4, -1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+K[j]) );
fi; 

else
	if (p+Position(P,K[i]+P[j]) = p+Position(P,P[i]+K[j])) then

if (t[Position(P,K[j]+P[i] )] = 1) then
Add(KK4, -1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+K[j]) );
else
Add(KK2, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]+P[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))

if (t[Position(P,P[i]+K[j] )] = 1) then
Add(KK2, -1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, P[i]+K[j]) );
Add(KK4, -1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, K[i]+P[j]) );
else
Add(KK2, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, P[i]+K[j]) );
Add(KK4, +1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, K[i]+P[j]) );
fi; 
	fi; 
fi; 
fi; fi;

##########

if (K[j]-P[i] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][p+Position(P,K[j])][2])=true) then
if (p+Position(P,P[j]-K[i]) > p+Position(P,K[j]-P[i])) then

if (t[Position(P,K[j]-P[i] )] = 1) then
Add(KK2, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, P[j]-K[i]) );
Add(KK4, -1*(-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, K[j]-P[i]) );
else
Add(KK2, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, P[j]-K[i]) );
Add(KK4, -1*(-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, K[j]-P[i]) );
fi; 
else
	if (p+Position(P,P[j]-K[i]) = p+Position(P,K[j]-P[i])) then

if (t[Position(P,P[j]-K[i] )] = 1) then
Add(KK4, -1*(-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, K[j]-P[i]) );
else
Add(KK2, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, P[j]-K[i]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))

if (t[Position(P,K[j]-P[i] )] = 1) then
Add(KK2, -1*(-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, K[j]-P[i]) );
Add(KK4, -1*(-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[j]-K[i]) );
else
Add(KK2, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, K[j]-P[i]) );
Add(KK4, +1*(-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[j]-K[i]) );
fi; 
	fi; 
fi; 
fi; fi;

##################
else # ovvero t[j]=-1
#########################

if (P[i]+P[j] in P ) then #qui so che aplha+beta < phi_aplha+phi_beta quindi sono nel caso sicuro
if not (IsEmpty(To[Position(P, P[i])][Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]+K[i]) > p+Position(P,P[j]+P[i])) then

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]+K[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+P[j]) );
else
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]+K[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+P[j]) );
fi; 
else
	if (p+Position(P,K[j]+K[i]) = p+Position(P,P[j]+P[i])) then
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+P[j]) );
else
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]+K[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, P[i]+P[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P,K[j]+K[i]) );
else
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, P[i]+P[j]) );
Add(KK4, +1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, K[j]+K[i]) );
fi; 
	fi; 
fi; 
fi; fi;

###############

if (P[j]-P[i] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][p+Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]-K[i]) > p+Position(P,P[j]-P[i])) then

if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, K[j]-K[i]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[j]-P[i]) );
else
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, K[j]-K[i]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[j]-P[i]) );
fi; 
else
	if (p+Position(P,K[j]-K[i]) = p+Position(P,P[j]-P[i])) then
if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[j]-P[i]) );
else
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, K[j]-K[i]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))

if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, P[j]-P[i]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, K[j]-K[i]) );
else
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, P[j]-P[i]) );
Add(KK4, +1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, K[j]-K[i]) );
fi; 
	fi; 
fi; 
fi; fi;

#############

if (P[i]+K[j] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][Position(P,K[j])][2])=true) then
if (p+Position(P,K[i]+P[j]) > p+Position(P,P[i]+K[j])) then

if (t[Position(P,P[i]+K[j] )] = 1) then
Add(KK2, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]+P[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+K[j]) );
else
Add(KK2, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]+P[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+K[j]) );
fi; 
else
	if (p+Position(P,K[i]+P[j]) = p+Position(P,P[i]+K[j])) then
if (t[Position(P,K[j]+P[i] )] = 1) then
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+K[j]) );
else
Add(KK2, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]+P[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))

if (t[Position(P,P[i]+K[j] )] = 1) then
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, P[i]+K[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, K[i]+P[j]) );
else
Add(KK2, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, P[i]+K[j]) );
Add(KK4, +1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, K[i]+P[j]) );
fi; 
	fi; 
fi; 
fi; fi;

#############

if (K[j]-P[i] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][p+Position(P,K[j])][2])=true) then
if (p+Position(P,P[j]-K[i]) > p+Position(P,K[j]-P[i])) then

if (t[Position(P,K[j]-P[i] )] = 1) then
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, P[j]-K[i]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, K[j]-P[i]) );
else
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, P[j]-K[i]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, K[j]-P[i]) );
fi; 
else
	if (p+Position(P,P[j]-K[i]) = p+Position(P,K[j]-P[i])) then
if (t[Position(P,K[j]-P[i] )] = 1) then
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, K[j]-P[i]) );
else
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, P[j]-K[i]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,K[j]-P[i] )] = 1) then
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, K[j]-P[i]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[j]-K[i]) );
else
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, K[j]-P[i]) );
Add(KK4, +1*(+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[j]-K[i]) );
fi; 
	fi; 
fi; 
fi; fi;

fi;

########################à
else #t[i]=-1
if (t[j]=1) then
###############

if (P[i]+P[j] in P ) then #qui so che aplha+beta < phi_aplha+phi_beta quindi sono nel caso sicuro
if not (IsEmpty(To[Position(P, P[i])][Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]+K[i]) > p+Position(P,P[j]+P[i])) then
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]+K[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+P[j]) );
else
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]+K[j]) );
Add(KK4, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+P[j]) );
fi; 
else
	if (p+Position(P,K[j]+K[i]) = p+Position(P,P[j]+P[i])) then
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+P[j]) );
else
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]+K[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, P[i]+P[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P,K[j]+K[i]) );
else
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, P[i]+P[j]) );
Add(KK4, +1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, K[j]+K[i]) );
fi; 
	fi; 
fi; 
fi; fi;

###############

if (P[j]-P[i] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][p+Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]-K[i]) > p+Position(P,P[j]-P[i])) then

if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, K[j]-K[i]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[j]-P[i]) );
else
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, K[j]-K[i]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[j]-P[i]) );
fi; 
else
	if (p+Position(P,K[j]-K[i]) = p+Position(P,P[j]-P[i])) then
if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[j]-P[i]) );
else
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, K[j]-K[i]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, P[j]-P[i]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, K[j]-K[i]) );
else
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, P[j]-P[i]) );
Add(KK4, +1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, K[j]-K[i]) );
fi; 
	fi; 
fi; 
fi; fi;

#############

if (P[i]+K[j] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][Position(P,K[j])][2])=true) then
if (p+Position(P,K[i]+P[j]) > p+Position(P,P[i]+K[j])) then
if (t[Position(P,P[i]+K[j] )] = 1) then
Add(KK2, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]+P[j]) );
Add(KK4, -1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+K[j]) );
else
Add(KK2, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]+P[j]) );
Add(KK4, -1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+K[j]) );
fi; 
else
	if (p+Position(P,K[i]+P[j]) = p+Position(P,P[i]+K[j])) then
if (t[Position(P,K[j]+P[i] )] = 1) then
Add(KK4, -1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+K[j]) );
else
Add(KK2, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]+P[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[i]+K[j] )] = 1) then
Add(KK2, -1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, P[i]+K[j]) );
Add(KK4, -1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, K[i]+P[j]) );
else
Add(KK2, (-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, P[i]+K[j]) );
Add(KK4, +1*(-1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, K[i]+P[j]) );
fi; 
	fi; 
fi; 
fi; fi;

#############

if (K[j]-P[i] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][p+Position(P,K[j])][2])=true) then
if (p+Position(P,P[j]-K[i]) > p+Position(P,K[j]-P[i])) then

if (t[Position(P,K[j]-P[i] )] = 1) then
Add(KK2, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, P[j]-K[i]) );
Add(KK4, -1*(-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, K[j]-P[i]) );
else
Add(KK2, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, P[j]-K[i]) );
Add(KK4, -1*(-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, K[j]-P[i]) );
fi; 
else
	if (p+Position(P,P[j]-K[i]) = p+Position(P,K[j]-P[i])) then
if (t[Position(P,K[j]-P[i] )] = 1) then
Add(KK4, -1*(-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, K[j]-P[i]) );
else
Add(KK2, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, P[j]-K[i]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,K[j]-P[i] )] = 1) then
Add(KK2, -1*(-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, K[j]-P[i]) );
Add(KK4, -1*(-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[j]-K[i]) );
else
Add(KK2, (-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, K[j]-P[i]) );
Add(KK4, +1*(-1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[j]-K[i]) );
fi; 
	fi; 
fi; 
fi; fi;

else #t[j]=-1

###############

if (P[i]+P[j] in P ) then #qui so che aplha+beta < phi_aplha+phi_beta quindi sono nel caso sicuro
if not (IsEmpty(To[Position(P, P[i])][Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]+K[i]) > p+Position(P,P[j]+P[i])) then

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]+K[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+P[j]) );
else
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]+K[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+P[j]) );
fi; 
else
	if (p+Position(P,K[j]+K[i]) = p+Position(P,P[j]+P[i])) then
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+P[j]) );
else
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]+K[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, P[i]+P[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P,K[j]+K[i]) );
else
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, P[i]+P[j]) );
Add(KK4, +1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, K[j]+K[i]) );
fi; 
	fi; 
fi; 
fi; fi;

#######

if (P[j]-P[i] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][p+Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]-K[i]) > p+Position(P,P[j]-P[i])) then

if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, K[j]-K[i]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[j]-P[i]) );
else
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, K[j]-K[i]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[j]-P[i]) );
fi; 
else
	if (p+Position(P,K[j]-K[i]) = p+Position(P,P[j]-P[i])) then
if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, P[j]-P[i]) );
else
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, K[j]-K[i]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, P[j]-P[i]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, K[j]-K[i]) );
else
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK2, p+Position(P, P[j]-P[i]) );
Add(KK4, +1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1])/2);
Add(KK4, p+Position(P, K[j]-K[i]) );
fi; 
	fi; 
fi; 
fi; fi;

######

if (P[i]+K[j] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][Position(P,K[j])][2])=true) then
if (p+Position(P,K[i]+P[j]) > p+Position(P,P[i]+K[j])) then

if (t[Position(P,P[i]+K[j] )] = 1) then
Add(KK2, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]+P[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+K[j]) );
else
Add(KK2, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]+P[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+K[j]) );
fi; 
else
	if (p+Position(P,K[i]+P[j]) = p+Position(P,P[i]+K[j])) then
if (t[Position(P,K[j]+P[i] )] = 1) then
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[i]+K[j]) );
else
Add(KK2, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, K[i]+P[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[i]+K[j] )] = 1) then
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, P[i]+K[j]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]+1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, K[i]+P[j]) );
else
Add(KK2, (+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, P[i]+K[j]) );
Add(KK4, +1*(+1*To[Position(P, P[i])][Position(P,K[j])][2][1]-1*To[Position(P, K[i])][Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, K[i]+P[j]) );
fi; 
	fi; 
fi; 
fi; fi;

#############

if (K[j]-P[i] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][p+Position(P,K[j])][2])=true) then
if (p+Position(P,P[j]-K[i]) > p+Position(P,K[j]-P[i])) then

if (t[Position(P,K[j]-P[i] )] = 1) then
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, P[j]-K[i]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, K[j]-P[i]) );
else
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, P[j]-K[i]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, K[j]-P[i]) );
fi; 
else
	if (p+Position(P,P[j]-K[i]) = p+Position(P,K[j]-P[i])) then
if (t[Position(P,K[j]-P[i] )] = 1) then
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, K[j]-P[i]) );
else
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, P[j]-K[i]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,K[j]-P[i] )] = 1) then
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, K[j]-P[i]) );
Add(KK4, -1*(+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[j]-K[i]) );
else
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK2, p+Position(P, K[j]-P[i]) );
Add(KK4, +1*(+1*To[Position(P, P[i])][p+Position(P,K[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,P[j])][2][1])/2);
Add(KK4, p+Position(P, P[j]-K[i]) );
fi; 
	fi; 
fi; 
fi; fi;

fi; fi;

if ( Position(P,K[i]) < Position(P,K[j]) ) then
SetEntrySCTable( TT, Position(P,P[i]),  p+Position(P, K[j]), KK2);
SetEntrySCTable( TT, Position(P,K[i]),  p+Position(P, K[j]), KK4);
else
 	if ( Position(P,K[j]) < Position(P,K[i]) ) then
SetEntrySCTable( TT, Position(P,P[i]),  p+Position(P, K[j]), KK2);
	fi;
fi;

else
	if ( Position(P, K[i]) = Position(P,P[i]) ) then #alpha< phi_alpha (ragiono con gli X)

KK2:=[];

if (t[i]=1) then
if (t[j]=1) then

if (P[i]+P[j] in P ) then #qui so che aplha+beta < phi_aplha+phi_beta quindi sono nel caso sicuro
if not (IsEmpty(To[Position(P, P[i])][Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]+K[i]) > p+Position(P,P[j]+P[i])) then

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[i]+K[j]) );
else
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[i]+K[j]) );
fi; 
else
	if (p+Position(P,K[j]+K[i]) = p+Position(P,P[j]+P[i])) then

if (t[Position(P,P[j]+P[i] )] = 1) then
else
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[i]+K[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[i]+P[j]) );
else
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[i]+P[j]) );
fi; 
	fi; 
fi; 
fi; fi;

###############

if (P[j]-P[i] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][p+Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]-K[i]) > p+Position(P,P[j]-P[i])) then

if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[j]-K[i]) );
else
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[j]-K[i]) );
fi; 
else
	if (p+Position(P,K[j]-K[i]) = p+Position(P,P[j]-P[i])) then

if (t[Position(P,P[j]-P[i] )] = 1) then
else
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[j]-K[i]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))

if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[j]-P[i]) );
else
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[j]-P[i]) );
fi; 
	fi; 
fi; 
fi; fi;

###########

else #t[j]=-1

if (P[i]+P[j] in P ) then #qui so che aplha+beta < phi_aplha+phi_beta quindi sono nel caso sicuro
if not (IsEmpty(To[Position(P, P[i])][Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]+K[i]) > p+Position(P,P[j]+P[i])) then

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[i]+K[j]) );
else
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[i]+K[j]) );
fi; 
else
	if (p+Position(P,K[j]+K[i]) = p+Position(P,P[j]+P[i])) then

if (t[Position(P,P[j]+P[i] )] = 1) then
else
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[i]+K[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[i]+P[j]) );
else
Add(KK2, (+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[i]+P[j]) );
fi; 
	fi; 
fi; 
fi; fi;

#######

if (P[j]-P[i] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][p+Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]-K[i]) > p+Position(P,P[j]-P[i])) then
if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[j]-K[i]) );
else
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[j]-K[i]) );
fi; 
else
	if (p+Position(P,K[j]-K[i]) = p+Position(P,P[j]-P[i])) then

if (t[Position(P,P[j]-P[i] )] = 1) then
else
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[j]-K[i]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[j]-P[i]) );
else
Add(KK2, (+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[j]-P[i]) );
fi; 
	fi; 
fi; 
fi; fi;

####################

fi;

else #t[i]=-1
if (t[j]=1) then

if (P[i]+P[j] in P ) then #qui so che aplha+beta < phi_aplha+phi_beta quindi sono nel caso sicuro
if not (IsEmpty(To[Position(P, P[i])][Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]+K[i]) > p+Position(P,P[j]+P[i])) then

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[i]+P[j]) );
else
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[i]+P[j]) );
fi; 
else
	if (p+Position(P,K[j]+K[i]) = p+Position(P,P[j]+P[i])) then
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[i]+P[j]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[i]+P[j]) );
else
Add(KK2, +1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[j]+K[i]) );
fi; 
	fi; 
fi; 
fi; fi;

############

if (P[j]-P[i] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][p+Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]-K[i]) > p+Position(P,P[j]-P[i])) then

if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[j]-P[i]) );
else
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[j]-P[i]) );
fi; 
else
	if (p+Position(P,K[j]-K[i]) = p+Position(P,P[j]-P[i])) then
if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[j]-P[i]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[j]-K[i]) );
else
Add(KK2, +1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[j]-K[i]) );
fi; 
	fi; 
fi; 
fi; fi;

else

if (P[i]+P[j] in P ) then #qui so che aplha+beta < phi_aplha+phi_beta quindi sono nel caso sicuro
if not (IsEmpty(To[Position(P, P[i])][Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]+K[i]) > p+Position(P,P[j]+P[i])) then

if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[j]+P[i]) );
else
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[j]+P[i]) );
fi; 
else
	if (p+Position(P,K[j]+K[i]) = p+Position(P,P[j]+P[i])) then
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[j]+P[i]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[j]+P[i] )] = 1) then
Add(KK2, -1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]+1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[i]+K[j]) );
else
Add(KK2, +1*(+1*To[Position(P, P[i])][Position(P,P[j])][2][1]-1*To[Position(P, K[i])][Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[i]+K[j]) );
fi; 
	fi; 
fi; 
fi; fi;

###############

if (P[j]-P[i] in P ) then 
if not (IsEmpty(To[Position(P, P[i])][p+Position(P,P[j])][2])=true) then
if (p+Position(P,K[j]-K[i]) > p+Position(P,P[j]-P[i])) then

if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[j]-P[i]) );
else
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[j]-P[i]) );
fi; 
else
	if (p+Position(P,K[j]-K[i]) = p+Position(P,P[j]-P[i])) then
if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, P[j]-P[i]) );
fi; 
	else  #(p+Position(P,K[j]+K[i]) < p+Position(P,P[j]+P[i]))
if (t[Position(P,P[j]-P[i] )] = 1) then
Add(KK2, -1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]+1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[j]-K[i]) );
else
Add(KK2, +1*(+1*To[Position(P, P[i])][p+Position(P,P[j])][2][1]-1*To[Position(P, K[i])][p+Position(P,K[j])][2][1]));
Add(KK2, p+Position(P, K[j]-K[i]) );
fi; 
	fi; 
fi; 
fi; fi;

fi;
fi;

SetEntrySCTable( TT, i,  p+Position(P,K[j]), KK2);

	fi;
fi; fi;

od; od;
##############################################################################
#########################################################################################


L:= LieAlgebraByStructureConstants( F0, TT );



K0:= [ ]; P0:= [ ];
H0:= [ ];
if h = () then

  for i in [1..Length(t)] do
      if t[i] = 1 then
         Add( K0, Basis(L)[i] );
      else
         Add( P0, Basis(L)[i] );
      fi;
  od;

  for i in [1..Length(t)] do
      if t[i] = 1 then
         Add( K0, Basis(L)[Length(t)+i] );
      else
         Add( P0, Basis(L)[Length(t)+i] );
      fi;
  od;

  Append( K0, Basis(L){[2*Length(t)+1..Length(Basis(L))]} );
  H0:= Basis(L){[2*Length(t)+1..Length(Basis(L))]};

else

  for i in [1..Length(t)] do
      pos:= Position( P, K[i] );
      if pos = i then
         if t[i] > 0 then
            Add( K0, Basis(L)[i] );
            Add( K0, Basis(L)[i+p] );
         else
            Add( P0, Basis(L)[i] );
            Add( P0, Basis(L)[i+p] );
         fi;
      elif pos > i then
         Add( K0, Basis(L)[i] );
         Add( P0, Basis(L)[pos] );        
         Add( K0, Basis(L)[i+p] );
         Add( P0, Basis(L)[pos+p] );
      fi;
  od;

   for i in [1..s] do
       if i^h = i then
          Add( K0, Basis(L)[ 2*Length(t)+i ] );
          Add( H0, Basis(L)[ 2*Length(t)+i ] ); 
       elif i^h > i then
          Add( K0, Basis(L)[ 2*Length(t)+i ] );
          Add( H0, Basis(L)[ 2*Length(t)+i ] ); 
          Add( P0, Basis(L)[ 2*Length(t)+i^h ] );
       fi;
   od;
fi;




  makeCartInv := function(L,K,P)
  local bas;
    bas := BasisNC(L,Concatenation(Basis(K),Basis(P)));
  return function(v)
     local k, p, cf, i;
        k   := Length(Basis(K));
        p   := Length(Basis(P));
        cf  := List(Coefficients(bas,v),x->x);
        for i in [k+1..k+p] do cf[i] := -cf[i]; od;
           return cf*bas;
        end;
   end; 

   K0 := SubalgebraNC(L,K0,"basis");
   P0 := SubspaceNC( L, P0 );
   
SetCartanDecomposition( L, rec( K:=K0 ,P:=P0, CartanInv := makeCartInv(L,K0,P0)  ) );
SetIsCompactForm( L, false );

SetCartanSubalgebra(L,SubalgebraNC(L,Basis(L){[2*p+1..2*p+s]},"basis") );
SetMaximallyCompactCartanSubalgebra( L, CartanSubalgebra(L) );
SetCartanSubalgebra( CartanDecomposition(L).K, SubalgebraNC(L,H0) );

# fare un elenco con tre elenchi: [x_alpha], [x_{-alpha}], [h_1,...,h_l], alpha > 0, tutti
# i vettori espressi in termini della base di L.
# Se b:= Basis(L); allora b[1] è il primo elemento della base, ecc. 

bb:=Basis(L);
BB:=[[],[],[]];
KK1:=0;
KK2:=0;
KK3:=0;
KK4:=0;

for i in [1..p] do
Add(BB[1], 0);
Add(BB[2], 0);
od;

for i in [1..s] do
Add(BB[3], 0);
od;

for j in [1..s] do

	if (OnPoints(j,h) > j) then
		BB[3][j]:= 1/2*One(F0)*(-1*E(4)*One(F0)*bb[2*p+j]+1*bb[2*p+OnPoints(j,h)]);
		BB[3][OnPoints(j,h)]:= -1/2*One(F0)*(+1*E(4)*One(F0)*bb[2*p+j]+1*bb[2*p+OnPoints(j,h)]);
	elif (OnPoints(j,h) = j) then
		BB[3][j]:= 1/2*One(F0)*(-1*E(4)*One(F0)*bb[2*p+j]);
	fi;
od;



for i in [1..p] do

	if (t[i] = 1) then

		if (Position(P, K[i]) > Position(P,P[i])) then
			KK1:=1/2*One(F0)*(bb[i]-1*E(4)*One(F0)*bb[Position(P,K[i])]); #X_alpha
			KK2:=1/2*One(F0)*(bb[i]+1*E(4)*One(F0)*bb[Position(P,K[i])]); #X_phi(alpha)
			KK3:=1/2*One(F0)*(bb[p+i]-1*E(4)*One(F0)*bb[p+Position(P,K[i])]); #Y_alpha
			KK4:=1/2*One(F0)*(bb[p+i]+1*E(4)*One(F0)*bb[p+Position(P,K[i])]); #Y_phi(alpha)
			BB[1][i]:= 1/2*One(F0)*(KK1-1*E(4)*One(F0)*KK3);
			BB[1][Position(P,K[i])]:= 1/2*One(F0)*(KK2-1*E(4)*One(F0)*KK4);
			BB[2][i]:=1/2*One(F0)*(-1*KK1-1*E(4)*One(F0)*KK3);
			BB[2][Position(P,K[i])]:=1/2*One(F0)*(-1*One(F0)*KK2-1*E(4)*One(F0)*KK4);
		elif (Position(P, K[i]) = Position(P,P[i])) then
			KK1:=1/2*One(F0)*(bb[i]); #X_alpha
			KK2:=1/2*One(F0)*(bb[p+i]); #Y_alpha
			BB[1][i]:= 1/2*One(F0)*(KK1-1*E(4)*One(F0)*KK2);
			BB[2][i]:= 1/2*One(F0)*(-1*KK1-1*E(4)*One(F0)*KK2);
		fi;
	else
		if (Position(P, K[i]) > Position(P,P[i]) ) then

                          
			KK1:=1/2*One(F0)*(bb[i]-1*E(4)*One(F0)*bb[Position(P,K[i])]); #X_alpha
			KK2:=1/2*One(F0)*(-bb[i]-1*E(4)*One(F0)*bb[Position(P,K[i])]); #X_phi(alpha)
			KK3:=1/2*One(F0)*(bb[p+i]-1*E(4)*One(F0)*bb[p+Position(P,K[i])]); #Y_alpha
			KK4:=1/2*One(F0)*(-bb[p+i]-1*E(4)*One(F0)*bb[p+Position(P,K[i])]); #Y_phi(alpha)
			BB[1][i]:= 1/2*One(F0)*(KK1-1*E(4)*One(F0)*KK3);
			BB[1][Position(P,K[i])]:=1/2*One(F0)*(KK2-1*E(4)*One(F0)*KK4);
			BB[2][i]:=1/2*One(F0)*(-1*KK1-1*One(F0)*E(4)*KK3);
                        
			BB[2][Position(P,K[i])]:=1/2*One(F0)*(-1*KK2-1*One(F0)*E(4)*KK4);
		elif (Position(P, K[i]) = Position(P,P[i]) ) then
			KK1:=-1*E(4)*One(F0)*1/2*(bb[i]); #X_alpha
			KK2:=-1*E(4)*One(F0)*1/2*(bb[p+i]); #Y_alpha
			BB[1][i]:= 1/2*One(F0)*(KK1-1*One(F0)*E(4)*KK2);
			BB[2][i]:= 1/2*One(F0)*(-1*KK1-1*E(4)*One(F0)*KK2);
		fi;
	fi;

od;


    rts:=[ ]; 
    for v in BB[1] do 
        sp:= Basis(SubspaceNC(L,[v],"basis"),[v]); 
        Add( rts, List( BB[3], t -> Coefficients(sp,t*v)[1] ) ); 
    od;

    

    R:= Objectify( NewType( NewFamily( "RootSystemFam", IsObject ),
                IsAttributeStoringRep and IsRootSystemFromLieAlgebra ), 
                rec() );
    SetCanonicalGenerators( R, [ BB[1]{[1..n]}, BB[2]{[1..n]}, BB[3] ] );
    SetUnderlyingLieAlgebra( R, L );
    SetPositiveRootVectors( R, BB[1] );
    SetNegativeRootVectors( R, BB[2] );

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

    allrts:= Concatenation( rts, -rts );
    fundr:= rts{[1..n]};
    SetCartanMatrix( R, List( fundr, x -> List( fundr, y -> CartInt( allrts, x, y ) ) ) );
    
    
  #roots are rationals
    if IsSqrtField(F0) then
      rts := List(rts, x-> List(x, SqrtFieldEltToCyclotomic));
    fi;

   if IsQQBarField(F0) then
      rts := List(rts, x-> List(x, corelg.QQBarRatToGap));
   fi; 
 

    SetPositiveRoots( R, rts );

   #Print("pos roots ",rts[1],rts[1][1] in SqrtField,"\n");


    SetNegativeRoots( R, -rts ); 
    SetSimpleSystem( R, rts{[1..n]} );

    SetRootSystem(L,R);
    SetRootSystem(MaximallyCompactCartanSubalgebra(L),R);

    SetChevalleyBasis( L, BB );

   #Print("still in sqrtfield? ",PositiveRoots(RootSystem(L))[1][1] in SqrtField,"\n");

    return L;

end;
