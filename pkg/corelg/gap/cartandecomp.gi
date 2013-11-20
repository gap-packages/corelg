## This file contains the functions to construct Cartan decompositions and
## maximally (non-)compact Cartan subalgebras
##
## functions contained in this file:
##   MaximallyCompactCartanSubalgebra
##   MaximallyNonCompactCartanSubalgebra
##   CartanDecomposition
##   corelg.mncptCSA
##   corelg.specialrtsys
##   RealStructure


#########################################################################
#
#
#
InstallMethod( RealStructure,
    "for a Lie algebra",
    true,
    [ IsLieAlgebra ], 0, 
 function(L)

    local bas, sigma;
    bas := ValueOption( "basis" );
    if not IsBasis(bas) then bas := Basis(L); fi;
    sigma := function(v) 
       return List(Coefficients(bas,v),ComplexConjugate)*bas; 
    end;
    return sigma;
 end );

#################################################################################
#
# a special function for computing a root system used in the construction of
# maximally (non-)compact CSA.
#
corelg.specialrtsys:= function( L, H, spaces, hh, h0 )

     # hh \cup {h0} is a Cartan subalgebra
     # spaces is a list of (bases of) subspaces of L, invariant under hh and h0,
     # decomposing them under h0 yields a root space dec.
     
    local F,          # coefficients domain of `L'
          BL,         # basis of `L'
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
          facs0, num, fam, f, b, c, r, F0, Mold, one, t1, t2, t3; 

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

    F   := LeftActingDomain( L );
    one := One(F);

    # First we compute the common eigenvectors of the adjoint action of a
    # Cartan subalgebra H. Here B will be a list of bases of subspaces
    # of L such that H maps each element of B into itself.
    # Furthermore, B has maximal length w.r.t. this property.

    BL   := Basis( L );
    B    := spaces;

    newB := [ ];
    for j in B do

        if Length(j) = 1 then
           Add( newB, j ); 
        else

           V    := Basis( VectorSpace( F, j, "basis" ), j );
           Mold := List( j, x -> Coefficients( V, h0*x ) );
           if IsSqrtField(F) then
              M    := SqrtFieldMakeRational(Mold);
              if M = false then 
                 #Error("matrix we want to compute char pol of cannot be made rationals");
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
              V := NullspaceMat( Value( l, Mold ) );
              Add( newB, List( V, x -> LinearCombination( j, x ) ) );
           od;
        fi;

  od;
  B:= newB;

  # Now we throw away the subspace H.
  #B:= Filtered( B, x -> ( not corelg.eltInSubspace(L,BasisVectors(Basis(H)),x[1])));
   B:= Filtered( B, x -> not x[1] in H );

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

   basH:= Basis(H);
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
end;

#########################################################################################


#########################################################################################
#
# constructs a maximally compact CSA
#
InstallMethod( MaximallyCompactCartanSubalgebra,
   "for a Lie algebra",
   true,
   [ IsLieAlgebra ], 0, 

function(L)
local F, sigma, H, R, pr, cb, cg, cbH, testIt, decomposeCSA, realRoot,
      decH, alpha, i, x,  y, j, K, prv, nrv, ev, rts, rt, rrr, rrv, spc, pos, cf, sx;

  
   Info(InfoCorelg,1,"start MaximallyCompactCartanSubalgebra");
   F  := LeftActingDomain(L);
   if not E(4)*One(F) in F then Error("need E(4) in field"); fi;
 
   sigma := RealStructure(L); 
   H     := CartanSubalgebra(L);

   if not ForAll(Basis(H),x->x=sigma(x)) then 
      Error("need basis of CSA which is fixed by sigma"); 
   fi;
   
   R   := RootsystemOfCartanSubalgebra(L,H);
   SetRootSystem(L,R);
   pr  := PositiveRoots(R);
   prv := PositiveRootVectors(R);
   nrv := NegativeRootVectors(R);
   #cb  := ChevalleyBasis(R);
   cg  := CanonicalGenerators(R);
   cbH := Basis(H,cg[3]);

   #### small test if all is compatible
   testIt := function()
   local C, basH, h, i, r, cf, cb;

      cb:= ChevalleyBasis(R);
      C   := LieCentraliser(L,H);
      if not IsAbelian(C) or not Dimension(C)=Dimension(H) then
         Error("not a CSA");
      fi;
      if not cg[3] = cb[3]{[1..Length(cg[3])]} then Error("error 1"); fi;
      basH := Basis(H,cg[3]);
      for h in cg[3] do
         for i in [1..Length(pr)] do      
            r  := pr[i]*Coefficients(basH,h);
            cf := Coefficients(BasisNC(SubspaceNC(L,[cb[1][i]],"basis"),[cb[1][i]]),h*cb[1][i])[1];
            if not r=cf then Error("error cf"); fi;
         od;
      od;
      Print("test ok\n");
   end;
   #####
  #testIt();


  #decompose H=(H\cap K) + (H\cap P), determined by the action of sigma
   decomposeCSA := function(H,cbH)
   local theta, ev, esp, HK, HP, h, ad;
      theta := List(cbH,x->Coefficients(cbH,-sigma(x)));
      ev     := [];
      ev[1] := List(NullspaceMat(theta-theta^0),x->x*cbH);
      ev[2] := List(NullspaceMat(theta+theta^0),x->x*cbH);  
      if Length(ev[2])=0 then
         Info(InfoCorelg,3,"  now have CSA with compact dimension ",Dimension(H));
         return rec(basHK := cbH, basHP := Basis(SubspaceNC(H,[],"basis")));
      elif Length(ev[1])=0 then
         Info(InfoCorelg,3,"  now have CSA with compact dimension 0");
         return rec(basHP := cbH, basHK := Basis(SubspaceNC(H,[],"basis")));
      fi;
      Info(InfoCorelg,3,"  now have CSA with compact dimension ",Length(ev[1]));
      return rec(basHK := ev[1], basHP := ev[2]);
   end;
   

  #return a real root alpha; this is a root which vanishes on HK
   realRoot := function(pr,basH,basHK)
   local alpha, real, h, cf;
      for alpha in pr do     
         cf  := List(basHK,h-> Coefficients(basH,h)*alpha);
         if ForAll(cf,x->x=Zero(F)) then return alpha; fi;
      od;
     return false;
   end;

   decH  := decomposeCSA(H,cbH);
   alpha := realRoot(pr,cbH,decH.basHK);

  #after the while-loop H is a maximally compact Cartan subalgebra
   while not alpha=false do
      i  := Position(pr,alpha);
      x  := prv[i];
      sx:= sigma(x);
      if sx <> x then
         if sx <> -x then
            x:= x+sx;
         else
            x:= E(4)*One(F)*x;
         fi;
      fi; 

      y  := nrv[i];
      cf := Coefficients( Basis( SubspaceNC( L, [x],"basis" ), [x] ), (x*y)*x )[1];
      y  := (2/cf)*y;

      #if sigma(x)=x then
      #   if not sigma(y)=y then Error("ups"); fi;
      #elif sigma(x)=-x then
      #   if not sigma(y)=-y then Error("ups"); fi;
      #   x :=  E(4)*One(F)*x;
      #   y := -E(4)*One(F)*y;
      #fi;
      if not (x*y)*x = 2*x or not (x*y)*y = -2*y then Error("not triple"); fi;
   
      K  := List(NullspaceMat(TransposedMat([alpha])),x->(One(F)*x)*cbH); 
      rts:= [ ]; spc:= [ ];
      rrv:= Concatenation( prv, nrv );
      rrr:= Concatenation( pr, -pr );
      for i in [1..Length(rrv)] do
         rt:= List( K, h-> Coefficients(cbH,h)*rrr[i] );
         pos:= Position( rts, rt );
         if pos = fail then
            Add( rts, rt );
            Add( spc, [ rrv[i] ] );
         else
            Add( spc[pos], rrv[i] );
         fi;
      od;

      pos:= Position( rts, List( K, h -> 0 ) );
      Append( spc[pos], Basis(H) );
     
      H   := SubspaceNC(L,Concatenation( K, [ x-y ] ),"basis");
    
     #use special method for computing RS; know already large part of CSA + RS!
     #R   := RootsystemOfCartanSubalgebra(L,H);
      R   := corelg.specialrtsys( L, H, spc, K, x-y ); 
      pr  := PositiveRoots(R);
      prv := PositiveRootVectors(R);
      nrv := NegativeRootVectors(R);   
      cg  := CanonicalGenerators(R);
      cbH := Basis(H,cg[3]);
      decH  := decomposeCSA(H,cbH);
      alpha := realRoot(pr,cbH,decH.basHK);
   od;
   Info(InfoCorelg,1,"end MaximallyCompactCartanSubalgebra; found CSA with compact dim ",Length(decH.basHK));

   SetcorelgCompactDimOfCSA(H,Length(decH.basHK));

   if not HasCartanSubalgebra(L) then 
      SetCartanSubalgebra(L,H); 
   fi;
   
   return H;
end );


##########################################################################
#
# function for computing a maximally non-compact CSA
# used by MaximallyNonCompactCartanSubalgebra
#
#
corelg.mncptCSA:= function(L)
local F, sigma, H, R, pr, cb, cg, cbH, testIt, decomposeCSA, imagRoot, nr, tau, h,yy,
      decH, alpha, i, x,  y, j, K, prv, nrv, ev, rts, rt, rrr, rrv, spc, pos, cf, srt;

  
   Info(InfoCorelg,1,"start MaximallyNonCompactCartanSubalgebra");
   F  := LeftActingDomain(L);
   if not E(4)*One(F) in F then Error("need E(4) in field"); fi;
   
   sigma := RealStructure(L); 
   H     := CartanSubalgebra(L);

   if not ForAll(Basis(H),x->x=sigma(x)) then 
      Error("need basis of CSA which is fixed by sigma"); 
   fi;
   
   R   := RootsystemOfCartanSubalgebra(L,H);
   SetRootSystem(L,R);
   pr  := PositiveRoots(R);
   prv := PositiveRootVectors(R);
   nrv := NegativeRootVectors(R);
   #cb  := ChevalleyBasis(R);
   cg  := CanonicalGenerators(R);
   cbH := Basis(H,cg[3]);

  #decompose H=(H\cap K) + (H\cap P), determined by the action of sigma
   decomposeCSA := function(H,cbH)
   local theta, ev, esp, HK, HP, h, ad;
      theta := List(cbH,x->Coefficients(cbH,-sigma(x))); #theta on h as -sigma
     #Print(theta,"\n");
      ev     := [];
      ev[1] := List(NullspaceMat(theta-theta^0),x->x*cbH);
      ev[2] := List(NullspaceMat(theta+theta^0),x->x*cbH);  
      if Length(ev[2])=0 then
         Info(InfoCorelg,3,"  now have CSA with compact dimension ",Dimension(H));
         return rec(basHK := cbH, basHP := Basis(SubspaceNC(H,[],"basis")));
      elif Length(ev[1])=0 then
         Info(InfoCorelg,3,"  now have CSA with compact dimension 0");
         return rec(basHP := cbH, basHK := Basis(SubspaceNC(H,[],"basis")));
      fi;
      Info(InfoCorelg,3,"  now have CSA with compact dimension ",Length(ev[1]));
      return rec(basHK := ev[1], basHP := ev[2]);
   end;

  #return a imag root alpha; this is a root which vanishes on HP
  #here NEED ALSO  non-compact, that is, g_\alpha in \fp 
   imagRoot := function(pr,basH,basHP)
   local alpha, real, h, cf,i,x,y,yy,tau;
      for i in [1..Length(pr)] do
         alpha := pr[i];     
         cf    := List(basHP,h-> Coefficients(basH,h)*alpha);
         if ForAll(cf,x->x=Zero(F)) then 
          ##check if we can norm everything (a.k.a. noncompact root)
             x   := prv[i];
             y   := sigma(prv[i]);
             yy  := nrv[i];
             cf  := Coefficients( Basis( SubspaceNC( L, [nrv[i]],"basis" ), 
                                                        [nrv[i]] ), (x*yy)*yy )[1];
             yy  := (-2/cf)*yy;
           ##now we have triple x,(x*yy),yy:
             if not (x*yy)*x=2*x or not (x*yy)*yy=-2*yy then Error("not first triple"); fi;
             tau := Coefficients( Basis( SubspaceNC( L, [yy],"basis" ), [yy] ), y )[1]; 
             if tau=ComplexConjugate(tau) and tau>0 then return [alpha,tau,i]; fi;           
         fi;
      od;
     return false;
   end;

   decH  := decomposeCSA(H,cbH);
   alpha := imagRoot(pr,cbH,decH.basHP);

   nr:=1;
  #after the while-loop H is a maximally compact Cartan subalgebra
   while not alpha=false do
      
      i   := alpha[3];
      x   := prv[i];
      y   := sigma(prv[i]);
      tau := alpha[2];
    
      srt := Sqrt(tau^-1);
      if not srt in F then Error("cannot do this over ",F); fi;
      x   := srt*x;
      y   := sigma(x);      
      if not (x*y)*x = 2*x or not (x*y)*y = -2*y then Error("not triple"); fi;
      
      K := List(NullspaceMat(TransposedMat([alpha[1]])),i->(One(F)*i)*cbH);
 
      rts:= [ ]; spc:= [ ];
      rrv:= Concatenation( prv, nrv );
      rrr:= Concatenation( pr, -pr );
      for i in [1..Length(rrv)] do
         rt:= List( K, h-> Coefficients(cbH,h)*rrr[i] );
         pos:= Position( rts, rt );
         if pos = fail then
            Add( rts, rt );
            Add( spc, [ rrv[i] ] );
         else
            Add( spc[pos], rrv[i] );
         fi;
      od;

      pos:= Position( rts, List( K, h -> 0 ) );
      Append( spc[pos], Basis(H) );
  
      H   := SubspaceNC(L,Concatenation( K, [ x+y ] ),"basis");
    
     #use special method for computing RS; know already large part of CSA + RS!
     #R   := RootsystemOfCartanSubalgebra(L,H);
      R   := corelg.specialrtsys( L, H, spc, K, x+y ); 
      pr  := PositiveRoots(R);
      prv := PositiveRootVectors(R);
      nrv := NegativeRootVectors(R);
      cg  := CanonicalGenerators(R);
      cbH := Basis(H,cg[3]);
  
      decH  := decomposeCSA(H,cbH);
      alpha := imagRoot(pr,cbH,decH.basHP);
      
   od;
   Info(InfoCorelg,1,"end MaximallyNonCompactCartanSubalgebra; found CSA with compact dim ",Length(decH.basHK));
  
   SetcorelgCompactDimOfCSA(H,Length(decH.basHK));

   if not HasCartanSubalgebra(L) then SetCartanSubalgebra(L,H); fi;
  
   return H;
end;


#########################################################################
#
#
#
InstallMethod( MaximallyNonCompactCartanSubalgebra,
    "for a Lie algebra",
    true,
    [ IsLieAlgebra ], 0, 
 function(L)

    local m, cs,csa; 

    m:= ValueOption( "method" );
    if m = "CayleyTransforms" then
       return corelg.mncptCSA(L);
    fi;
    cs  := CartanSubspace(L);
    csa := CartanSubalgebra( LieCentralizer( L, cs ) );
    SetcorelgCompactDimOfCSA(csa,Dimension(csa)-Dimension(cs));
    if not HasCartanSubalgebra(L) then 
       SetCartanSubalgebra(L,csa); 
    fi;
    return csa;

 end );

############################################################################
##
#F   CompactDimensionOfCartanSubalgebra( <L> )
#F   CompactDimensionOfCartanSubalgebra( <L>, <H> )
##
##   <L> is a semisimple lie algebra over Gaussian rationals or SqrtField;
##   this function returns the compact dimension of <H>, and
##   of CartanSubalgebra(<L>) if <H> is not provided
##
InstallGlobalFunction( CompactDimensionOfCartanSubalgebra, function( arg ) 
local L, H, sigma, tmp, cbH, cg;

   L := arg[1];
   if Length(arg)=2 then H := arg[2]; else H := CartanSubalgebra(L); fi;
   if HascorelgCompactDimOfCSA(H) then return corelgCompactDimOfCSA(H); fi;
   sigma := RealStructure(L); 
   cg    := CanonicalGenerators(RootsystemOfCartanSubalgebra(L,H));
   cbH   := BasisNC(H,cg[3]);
   tmp   := List(cbH,x->Coefficients(cbH,-sigma(x)));
   tmp   := Length(NullspaceMat(tmp-tmp^0));
   SetcorelgCompactDimOfCSA(H,tmp);
   return tmp;
end);


############################################################################################
InstallMethod( CartanDecomposition,
   "for a Lie algebra",
   true,
   [ IsLieAlgebra ], 0, 
function(L)

local csa, h, R, cb, sigma, hs, es, found, h0, vals, pr, posr, i, sums, base, B,
      C, ct, en, newcg, r, s, fs, new, theta, im, pos, F, cf, esp, mat,iso,cd, tmp,
      decomposeCSA, cg, rH, bas;

   Info(InfoCorelg,1,"start CartanDecomposition");
   F     := LeftActingDomain(L);
   h     := MaximallyCompactCartanSubalgebra(L);
   R     := RootsystemOfCartanSubalgebra( L, h );
   cb    := corelg.myChevalleyBasis(L,R);
   sigma := RealStructure(L);
   cg    := CanonicalGenerators(R);

   decomposeCSA := function(H,cbH)
   local theta, ev, esp, HK, HP, h, ad;
      theta := List(cbH,x->Coefficients(cbH,-sigma(x)));
      ev    := Eigenvalues(F,theta);
      if ev = [1]*One(F) then
         return rec(basHK := cbH, basHP := Basis(SubspaceNC(H,[],"basis")));
      elif ev = [-1]*One(F) then
         return rec(basHP := cbH, basHK := Basis(SubspaceNC(H,[],"basis")));
      fi;
      esp   := Eigenspaces(F,theta);
      HK    := List(Basis(esp[Position(ev,1*One(F))]),x->x*cbH);
      HP    := List(Basis(esp[Position(ev,-1*One(F))]),x->x*cbH);
      return rec(basHK := HK, basHP := HP);
   end;

   rH:= decomposeCSA( h, Basis(h,cg[3]) );
   
  #find h0 to define new root ordering
   hs    := rH.basHK;
   if not ForAll(hs,x->x in h) then Error("ups..CSA"); fi;
   found := false;
   es    := PositiveRootVectors(R);
   while not found do 
      h0 := Sum( hs, h -> Random([-100..100])*h );
      if ForAll( es, x -> not IsZero( h0*x ) ) then found:= true; fi;
   od;

   #find new basis of simple roots (def by root ordering induced by h0)
    vals := List( es, x -> Coefficients( Basis( SubspaceNC( L, [x],"basis" ), [x] ), h0*x )[1] );
    pr   := PositiveRootsNF(R);
    posr := [ ];
    for i in [1..Length(pr)] do
       if vals[i] > vals[i]*0 then Add( posr, pr[i] ); else Add( posr, -pr[i] ); fi; ###^0
    od;
    sums := [];
    for r in posr do for s in posr do Add( sums, r+s ); od; od;
    base := Filtered( posr, x -> not x in sums );
    B    := BilinearFormMatNF(R);
    C    := List( base, x -> List( base, y -> 2*(x*B*y)/(y*B*y) ) );
    ct   := CartanType(C);
    en   := Concatenation( CartanType(C).enumeration );
    base := base{en};
     
   #now construct corresponding canonical generators
   newcg := corelg.makeCanGenByBase(pr,cb,base);
   es    := newcg[1];
   fs    := newcg[2];
   hs    := List([1..Length(es)],x->es[x]*fs[x]);

   pos := Set(List([1..Length(hs)],i->Set([i,Position(hs,-sigma(hs[i]))])));
   new := [[],[],[]];
   cf  := [];
   for i in pos do
      if Length(i) = 1 then 
         new[1][i[1]] := newcg[1][i[1]];
         new[2][i[1]] := newcg[2][i[1]];
         new[3][i[1]] := newcg[3][i[1]];
         cf[i[1]] := Coefficients(Basis(SubspaceNC(L,[newcg[2][i[1]]],"basis"),
                                 [newcg[2][i[1]]]),sigma(newcg[1][i[1]]))[1];
      else
         cf[i[1]]     := One(F);
         cf[i[2]]     := One(F);
         new[1][i[1]] := sigma(newcg[2][i[2]]);
         new[2][i[1]] := sigma(newcg[1][i[2]]);
         new[3][i[1]] := newcg[3][i[1]];
         new[1][i[2]] := newcg[1][i[2]];
         new[2][i[2]] := newcg[2][i[2]];
         new[3][i[2]] := newcg[3][i[2]];
      fi;
   od;
   if F=SqrtField then
      cf := List(cf,x-> -One(F)* SignInt(SqrtFieldMakeRational(x)));
   else
      cf := List(cf,x-> -SignInt(x));
   fi;
   im := [[],[]];
   for i in pos do
      if Length(i)=1 then
         im[1][i[1]] := cf[i[1]]*new[1][i[1]];
         im[2][i[1]] := cf[i[1]]*new[2][i[1]];
      else
         im[1][i[1]] := cf[i[1]]*new[1][i[2]];
         im[1][i[2]] := cf[i[2]]*new[1][i[1]];
         im[2][i[1]] := cf[i[1]]*new[2][i[2]];
         im[2][i[2]] := cf[i[2]]*new[2][i[1]];
      fi;
   od;
   im[3] := List([1..Length(new[1])],x->im[1][x]*im[2][x]);
   theta := LieAlgebraIsomorphismByCanonicalGenerators(L,new,L,im);  
   mat   := List(Basis(L),x->Coefficients(Basis(L),Image(theta,x)));    
   
   esp := [];
   esp[1] := List(NullspaceMat(mat-mat^0),x->x*Basis(L));
   esp[2] := List(NullspaceMat(mat+mat^0),x->x*Basis(L));   

   bas := Basis(L,Concatenation(esp[1],esp[2]));
   theta := function(v)
   local k, p, cf, i;
      k   := Length(esp[1]);
      p   := Length(esp[2]);
      cf  := List(Coefficients(bas,v),x->x);
      for i in [k+1..k+p] do cf[i] := -cf[i]; od;
      return cf*bas;
   end;

   tmp  := SubalgebraNC(L,esp[1],"basis");
   SetCartanSubalgebra(tmp,SubalgebraNC(tmp,rH.basHK,"basis")); 
   cd   := rec(CartanInv:=theta, K:=tmp, P := SubspaceNC(L,esp[2],"basis"));
   SetCartanDecomposition(L,cd);

   Info(InfoCorelg,1,"end CartanDecomposition");
   return cd;

end);


