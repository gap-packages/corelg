##  Main functions:
##
##
#F   RealWeylGroup( <L> )
#F   RealWeylGroup( <L>, <H> )
##
##   <L> is a semisimple lie algebra over Gaussian rationals or SqrtField;
##   <H> is a Cartan subalgebra of <L>; returns the real Weyl group wrt <H>
##   if <H> is not given as input, then CartanSubalgebra(<L>) is used
##



################################################
# returns the coroots of the root sys of L wrt H
#
corelg.coroots:= function(L,H)
local R, ch, h;

   R:= RootsystemOfCartanSubalgebra( L, H );
   ch:= ChevalleyBasis(R);
   h:= List( [1..Length(ch[1])], i -> ch[1][i]*ch[2][i] );
   Append( h, -h ); # for the negative roots...
   return h;

end;

######################
# i : index of a pos root, return the perm of the corresp reflection
#
corelg.refl:= function( R, i ) 
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
# test whether the i-th pos root of R is real, L a real Lie algebra, R root sys wrt CSA H
corelg.isreal:= function( L, HK, R, i ) 
    local ch, x;
    ch := ChevalleyBasis(R);
    x  := ch[1][i];
    return ForAll( HK, h -> IsZero(h*x) );
end;

######################
# test whether the i-th pos root of R is imaginary, L a real Lie algebra, R root sys wrt CSA H

corelg.isimag:= function( L, HP, R, i ) 
    local cd, ch, x;
    ch:= ChevalleyBasis(R);
    x:= ch[1][i];
    return ForAll( HP, h -> IsZero(h*x) );
end;

######################
# inds : indices of pos rts of a subsys, find basis
#
corelg.base:= function(R,inds) 
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


############################################################################
##
#F   RealWeylGroup( <L> )
#F   RealWeylGroup( <L>, <H> )
##
##   <L> is a semisimple lie algebra over Gaussian rationals or SqrtField;
##   <H> is a Cartan subalgebra of <L>; returns the real Weyl group wrt <H>
##   if <H> is not given as input, then CartanSubalgebra(<L>) is used. Here
##   the real Weyl group is N_G(H) / Z_G(H) where G is the connected component
##   of the group of real points of the complex adjoint group of L.
##
InstallGlobalFunction(  RealWeylGroup, function( arg )

local R, re, br, gr, im, bi, gi, ch, hr, hi, ex, bex, gex, h, theta, p, 
    gex0, cpt, gcpt, bcpt, HP, HK, rts, C, rank, imgs, imfund, A, sol, B,
    ips, ip, i, j, val, gens, d, r, P, Q, S, pr, forA, rA, rho, res, L, H,
    tmp, phiImRho, phiImRhoBase, imphi;

    L := arg[1];
    if Length(arg) = 2 then H := arg[2]; else H := CartanSubalgebra(L); fi;

    if HascorelgRealWG(H) then return corelgRealWG(H); fi;


    HK := BasisVectors( Basis( Intersection( H, CartanDecomposition( L ).K ) ) );
    HP := BasisVectors( Basis( Intersection( H, CartanDecomposition( L ).P ) ) );
    R  := RootsystemOfCartanSubalgebra( L, H );

    ## do everything wrt simple roots
    B := BilinearFormMatNF(R);
    pr:= PositiveRootsNF(R);

    ## real roots, basis, and reflections
    re := Filtered([1..Length(PositiveRoots(R))], i-> corelg.isreal( L, HK, R, i ));
    br := corelg.base( R, re );
    gr := List( br, x-> corelg.refl( R, x ));

    ## imaginary roots, basis, and reflections
    im := Filtered([1..Length(PositiveRoots(R))], i-> corelg.isimag( L, HP, R, i ));
    bi := corelg.base( R, im );
    gi := List( bi, x-> corelg.refl( R, x ));

    ## compact roots (those with x_\alpha \in K
    ch   := ChevalleyBasis( R );
    cpt  := Filtered( im, i->ch[1][i] in CartanDecomposition( L ).K);
    bcpt := corelg.base( R, cpt );
    gcpt := List( bcpt, x-> corelg.refl( R, x ));


 ## now look forA (preliminary stuff)
 ## here forA will be a basis of Phi^{im,phi}... it's superorthogonal
 ## we should have W(forA) = W(B) where B is the superorthogonal set of simple rts

   if cpt = [] then
      forA := im;
   else 
      rho  := Sum( pr{cpt} );
     ## note that imphi consists of pos roots in Phi^{im,rho}
     ## and then we simply compute a base by using corelg.base
     ## the theory says this base is superorthogonal, so modulo signs 
     ## this is unique... 
      imphi := Filtered(im, i -> IsZero( pr[i]*B*rho ));
      forA  := corelg.base(R,imphi);
   fi;

 ## now get (W^c)^\theta

    ## now ex consists of positive roots in Phi^c
    hr := Sum( re, i -> pr[i] );
    hi := Sum( im, i -> pr[i] );
    ex := Filtered( [ 1 .. Length( pr ) ], i -> IsZero(pr[i]*B*hr) and IsZero(pr[i]*B*hi) );

    ## check which of those reflections commute with theta; 
    ## now gex will be generators for (W^c)^\theta
    ## note h_\alpha = [x_\alpha,x_{-\alpha}] mapped to h_{theta(\alpha)}

    h     := corelg.coroots( L, H );   
    theta := CartanDecomposition( L ).CartanInv;
    if Length( ex ) = 0 then
        gex:= [ ]; 
    else
        bex := corelg.base( R, ex );
        gex := List( bex, x->corelg.refl( R, x ));
        p   := PermList( List( h, x-> Position( h, theta( x ) )));
        gex := GeneratorsOfGroup( Centralizer( Group( gex ), p ) );
    fi;


  ## if forA=[] then A=1 and we are done
  ## otherwise we need to compute A

    if Length(forA) = 0 then
       res:= [ gex, [], gcpt, gr ];
    else

      ## first get matrix describing how theta permutes the simple roots
      ## this mat is linear map given wrt simple roots
       rts := Concatenation( PositiveRootsNF(R), -PositiveRootsNF(R) );
       C   := CartanMatrix(R);
       rank:= Length( C );
       imgs := List( [1..rank], i -> rts[ Position( h, theta( h[i] ) ) ] );


      ## now get a integer basis of P^\theta 
      ## weight space is generated by fundamental roots \lambda_i
      ## Humphrey p68: Cartan C is base change from alpha_i to lambda_i
      ## so imfund is action of theta in terms of fund weights
     
      
     ### LET'S DO NullspaceMatIntMat instead of NullspaceMat... 
     ### since everything is given wrt fundamendal roots, we really only
     ### want integer lincombs of those in solution space
       imfund := C^-1*imgs*C;
       imfund := imfund;
       A      := imfund-imfund^0;
       sol    := NullspaceIntMat( A );  

   
       ## now set up equation system for the sum to be equal 0 mod 2

       ips := [ ];
       for i in [1..Length(sol)] do               
           ip  := [ ];
           for j in [1..Length(forA)] do
             ### first, take solution vector and write it wrt simple roots     
               tmp := sol[i]*C^-1;  
             ### now compute <tmp, forA>
               val := (2*tmp*B*rts[forA[j]]) / (rts[forA[j]]*B*rts[forA[j]]); 
               Add( ip, val );
           od;
           Add( ips, ip );
       od;

     ## now get basis of solution space and then the corresponding elements
     ## need to transpose - and look at everything mod 2:
       sol := NullspaceMat( TransposedMat(ips)*One(GF(2)));
       sol := List( sol, u -> List( u, Int ) );
       rA  := List( forA, i -> corelg.refl( R, i ) );
       gens:= List( sol, u -> Product( [1..Length(u)], i -> rA[i]^u[i] ) );
       res:= [ gex, gens, gcpt, gr ];
    fi;

    tmp := Group( Concatenation( res ), () );
    tmp!.ex := res[1];
    tmp!.a  := res[2];
    tmp!.r  := res[4];
    tmp!.ci := res[3];
    SetcorelgRealWG(H,tmp);
    return tmp;

end);



StructureOfRealWeylGroup := function(W)
   if not IsBound(W!.ex) then
       Error("The input group has not been constructed via RealWeylGroup");
   fi;
   Print("The real Weyl group is  Wc  |x  ( (A  |x  Wci) x Wr )   where \n");
   Print("Wc is of size ",Size( Group(W!.ex,()) ),"\n" );
   Print("A is of size ",Size( Group(W!.a,()) ),"\n" );
   Print("Wci is of size ",Size( Group(W!.ci, () )),"\n" );
   Print("Wr is of size ",Size( Group(W!.r,()) ),"\n" );
end;




