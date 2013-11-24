corelg.isom:= function( L, H1, H2, grad )

    # MUST have f that respects grading!!!
    # For Z-gradings, take in both cases pos roots that have pos degree...

    local b1, b2, c1, c2, R1, R2, t, tp, en, i;

    R1:= RootsystemOfCartanSubalgebra(L,H1);
    R1:= corelg.rtsys_withgrad( L, Concatenation(PositiveRootVectors(R1),NegativeRootVectors(R1)),
                        H1, grad );
    R2:= RootsystemOfCartanSubalgebra(L,H2);
    R2:= corelg.rtsys_withgrad( L, Concatenation(PositiveRootVectors(R2),NegativeRootVectors(R2)),
                        H2, grad );
    

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

    b1:= SLAfcts.canbas( L, c1 );
    b2:= SLAfcts.canbas( L, c2 );

    return AlgebraHomomorphismByImagesNC( L, L, Flat(b1), Flat(b2) );

end;

Zgrading:= function( L, d )

    # L is semisimple, split, over SqrtField.

    local R, C, h, gr, K0, D, H, cc, cd;

    R:= RootSystem(L);
    C:= CartanMatrix(R);
    h:= (C^-1*d)*ChevalleyBasis(L)[3];
    gr:= SL2Grading( L, h );  # INEFFICIENT!

    K0:= Subalgebra( L, gr[3] );
  
    cd:= CartanDecomposition(L);
    # check...
    if ForAny( Basis(K0), x -> not cd.CartanInv(x) in K0 ) then
       Print("ERROR!!!! K0 not stable under Cartan inv.\n");
    fi; 

    D:= LieDerivedSubalgebra( K0 );
    SetCartanDecomposition( D, rec( CartanInv:= cd.CartanInv,
                                    K:= Intersection( cd.K, D ),
                                    P:= Intersection( cd.P, D ) ) );
    H:= CartanSubalgebras(D);
    cc:= BasisVectors( Basis( LieCentre(K0) ) );
    SetCartanSubalgebras( K0, List( H, U -> 
                          Subalgebra(K0,Concatenation(BasisVectors(Basis(U)),cc))));
    return rec( grading:= rec( g0:= gr[3], gp:= gr[1], gn:= gr[2] ), K0:= K0 );

end;

SubWeylGens:= function( L, H, K0 )

    # H a CSA contained in K0, get: root system of K0, i.e., roots of L wrt H
    # that are contained in K0, get simple roots, and corr permutations of the
    # roots of L wrt H.
    # Assume that K0 is theta stable!

    local R, p, rv, pr0, i, j, sums, bas0, perms, B, rts, norm, lst, rt, ch, h, theta,
          pth; 

    R:= RootsystemOfCartanSubalgebra( L, H );
    p:= PositiveRootsNF(R);
    rv:= PositiveRootVectors(R);
    pr0:= [ ];
    for i in [1..Length(p)] do 
        if rv[i] in K0 then
           Add( pr0, p[i] );
        fi;
    od;

    sums:= [ ];
    for i in [1..Length(pr0)] do
        for j in [1..Length(pr0)] do
            Add( sums, pr0[i]+pr0[j] );
        od;
    od;

    bas0:= Filtered( pr0, x -> not x in sums );

    perms:= [ ];
    B:= BilinearFormMatNF(R);
    rts:= Concatenation( p, -p );
    for i in [1..Length(bas0)] do
        norm:= bas0[i]*B*bas0[i];
        lst:= [ ];
        for j in [1..Length(rts)] do
            rt:= rts[j] - (2*(rts[j]*B*bas0[i])/norm)*bas0[i];
            lst[j]:= Position( rts, rt );
        od;
        Add( perms, PermList( lst ) );
    od;


    ch:= ChevalleyBasis(R);
    h:= List( [1..Length(ch[1])], i -> ch[1][i]*ch[2][i] );
    Append( h, -h ); # for the negative roots...
    theta:= CartanDecomposition(L).CartanInv;
    pth:= PermList( List(h,x -> Position(h,theta(x)) ) );

    return rec( W0:= perms, th:= pth );   

end;


imagesCarr:= function( L, H2, f, car, grp, C )

    # f is an autom  L --> L mappring the root vecs rel to H1 to the root vecs
    # rel to H2. car is a carrier algebra wrt H1 (so a record with comps g0,
    # negdeg, posdeg. grp are permutations in the W0 rel to H2. 
    # we compute all conjugates of car that are defined over R, rel to the elts in grp.

    # here we assume that grp is a set of coset reps wrt subgroup C
    # (so we only use the next function!). 

    # output is a record with fields car, the carrier algebras, and conjs: for each
    # elt in car a list of conjugates under C.

    local R, rv, bas, posi, negi, i, j, k, p, n, cf, ind, g0i, algs, sig, g, 
          p0, g0, n0, U, y, conjs, conj0, h, c;

    R:= RootsystemOfCartanSubalgebra( L, H2 );
    rv:= Concatenation( PositiveRootVectors(R), NegativeRootVectors(R) );
    bas:= Basis( Subspace( L, rv ), rv );

    posi:= [ ];
    negi:= [ ];
    
    for i in [1..Length(car.posdeg)] do
        p:= [ ]; n:= [ ];
        for j in [1..Length(car.posdeg[i])] do
            cf:= Coefficients( bas, Image(f, car.posdeg[i][j] ) );
            # check...
            ind:= 0;
            for k in [1..Length(cf)] do
                if not IsZero(cf[k]) then
                   if ind > 0 then
                      Print("ERROR!!!! more than one nonzero.\n");
                   else
                      ind:= k;
                   fi;
                fi;
            od;
            Add( p, ind );

            cf:= Coefficients( bas, Image(f, car.negdeg[i][j] ) );
            # check...
            ind:= 0;
            for k in [1..Length(cf)] do
                if not IsZero(cf[k]) then
                   if ind > 0 then
                      Print("ERROR!!!! more than one nonzero.\n");
                   else
                      ind:= k;
                   fi;
                fi;
            od;
            Add( n, ind );
        od;
        Add( posi, p ); Add( negi, n );
    od;

    g0i:= [ ];
    for j in [1..Length(car.g0)] do
        y:= Image( f, car.g0[j] );
        if not y in H2 then

           cf:= Coefficients( bas, y );
           # check...
           ind:= 0;
           for k in [1..Length(cf)] do
               if not IsZero(cf[k]) then
                   if ind > 0 then
                      Print("ERROR!!!! more than one nonzero.\n");
                   else
                      ind:= k;
                   fi;
               fi;
           od;
           Add( g0i, ind );
        fi;
    od;

    algs:= [ ];
    sig:= function(u)
       return List( Coefficients( Basis(L), u ), ComplexConjugate )*Basis(L); 
    end;

    conjs:= [ ];
    for g in grp do
        p0:= List( posi, u -> List( u, k -> rv[k^g] ) );
        n0:= List( negi, u -> List( u, k -> rv[k^g] ) );
        g0:= List( g0i, k -> rv[k^g] );
        U:= Subalgebra( L, Concatenation( Flat(g0), Flat(p0), Flat(n0) ) );
        if ForAll( Basis(U), x -> sig(x) in U ) then
           # found one...
           Add( algs, rec( U:= U, posdeg:= p0, negdeg:= n0, g0:= g0 ) );
           conj0:= [ ];
           for c in C do
               h:= g*c;  
               p0:= List( posi, u -> List( u, k -> rv[k^h] ) );
               n0:= List( negi, u -> List( u, k -> rv[k^h] ) );
               g0:= List( g0i, k -> rv[k^h] );
               Add( conj0, rec( posdeg:= p0, negdeg:= n0, g0:= g0 ) );
           od;
           Add( conjs, conj0 );
        fi;
    od;
           
    return rec( cars:= algs, conjs:= conjs );

end;

carrAlgs:= function( L, H1, H2, rc, car )

    # car: a carrier alg wrt H1, find all carr algs rel to H2 that are images of car,
    # up to equivalence.
    # rc: output of Zgrading
    # For the moment we ASSUME that the real Weyl group of H2 is the centralizer of theta...

    local f,sb,W0,S,cos,elms,crs,crs0,i,N,NP,NK,HN,mcpt,s0,bag,equal,CN,DN,HD;

    equal:= function( c1, c2 )  # for carrier algs

         local i, V1, V2; 

         for i in [1..Length(c1.posdeg)] do
             V1:= Subspace( L, c1.posdeg[i] );
             V2:= Subspace( L, c2.posdeg[i] );
             if V1 <> V2 then return false; fi;
         od;

         for i in [1..Length(c1.negdeg)] do
             V1:= Subspace( L, c1.negdeg[i] );
             V2:= Subspace( L, c2.negdeg[i] );
             if V1 <> V2 then return false; fi;
         od;
         V1:= Subspace( L, c1.g0 );
         V2:= Subspace( L, c2.g0 );
         if V1 <> V2 then return false; fi;
         return true;

    end;
         
 

    f:= corelg.isom(L, H1, H2, rc.grading );
    sb:= SubWeylGens( L, H2, rc.K0 );
    W0:= Group(sb.W0);
    S:= Centralizer(W0,sb.th);
    cos:= CosetDecomposition( W0, S );
    cos:= List( cos, a -> a[1]^-1*a*a[1] );
    elms:= List( cos, x -> x[1] );
    s0:= imagesCarr( L, H2, f, car, elms, S );
    crs:= s0.cars;
Print("carrier candidates... ",Length(crs),"\n");
    crs0:= [ ];
    bag:= [ ]; # big bag of carrier algs    
    for i in [1..Length(crs)] do
        if ForAll( bag, c1 -> not equal(c1,crs[i]) ) then 
           N:= LieNormalizer( L, crs[i].U );
           NP:= Intersection( N, CartanDecomposition(L).P );
           NK:= Intersection( N, CartanDecomposition(L).K );
           if Dimension(NP) + Dimension(NK) <> Dimension(N) then
              Print("Error normaliser not stable under theta!!!!\n");
           fi;
           SetCartanDecomposition( N, rec( CartanInv:= CartanDecomposition(L).CartanInv,
                         P:= NP, K:= NK ) );
           CN:= LieCentre(N);
           if Dimension(CN) > 0 then
              DN:= LieDerivedSubalgebra(N);
              SetCartanDecomposition( DN, rec( CartanInv:= CartanDecomposition(L).CartanInv,
                         P:= Intersection( DN, CartanDecomposition(L).P ), 
                         K:= Intersection( DN, CartanDecomposition(L).K ) ) );
              HD:= CartanSubalgebras(DN);
              SetCartanSubalgebras( N, List( HD, U -> Subalgebra( N, Concatenation( Basis(U),
                         Basis(CN) ) ) ) );
           fi;
           HN:= CartanSubalgebras(N);
           mcpt:= Maximum( List( HN, u -> Dimension(Intersection(u,NP)) ) );
Print( List( HN, u -> Dimension(Intersection(u,NP))), " ",Dimension(Intersection(H2,NP)),"\n" ); 
           if mcpt = Dimension(Intersection(H2,NP)) then
              Add( crs0, crs[i] );
           fi;
           Append( bag, s0.conjs[i] );
        fi; 
    od;

    return crs0;

end;

ZgradOrbs:= function( L, grading )

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
         Add( cars, rec( g0:= g0, posdeg:= gp, negdeg:= gn ) );

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


realOrbits:= function( type, rank, grading ) # At the moment only for Z-gradings!

   local L0, r, sl2, car, L, toSqrtField, i, j, c, zg, H, cc;

   L0:= SimpleLieAlgebra( type, rank, Rationals );
   r:=  ZgradOrbs( L0, grading );

   sl2:= [ ];
   car:= [ ];
   L:= SimpleLieAlgebra( type, rank, SqrtField );

   toSqrtField:= function( lst )
      return List( lst, x -> (Coefficients(Basis(L0),x)*One(SqrtField))*Basis(L) );
   end;

   for i in [1..Length(r.sl2)] do
       Add( sl2, toSqrtField(r.sl2[i]) );
       
       c:= rec( g0:= toSqrtField(r.carr[i].g0),
                posdeg:= List( r.carr[i].posdeg, toSqrtField ),
                negdeg:= List( r.carr[i].negdeg, toSqrtField ) );
       Add( car, c );
   od;         

   Print("Computed the orbits...\n");

   zg:= Zgrading( L, grading );
   Print("Computed the Z-grading...\n");

   H:= CartanSubalgebras( zg.K0 );

   if not H[1] = CartanSubalgebra(L) then 
      Print("Something WRONG with the Cartan subalgs of K0...\n");
   fi;

   Print("Start computing the other carrier algebras...\n");

   cc:= [ ];
   for i in [1..Length(car)] do 

       c:= [ [ car[i] ] ];
       for j in [2..Length(H)] do
           Add( c, carrAlgs( L, H[1], H[j], zg, car[i] ) );
       od;
       Add( cc, c );
   od;

   return rec( L:= L, sl2:= sl2, cars:= cc );


end;
