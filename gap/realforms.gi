# This file contains functions for constructing real simple LAs, Vogan diagrams, Satake diagrams, etc
#
# These are the functions contained here:
#
#  CartanSubalgebrasOfRealForm
#  CartanSubspace
#  VoganDiagram
#  SatakeDiagram
#  IdRealForm
#  RealFormById
#  NumberRealForms
#  AllRealForms
#  RealFormsInformation
#  IsomorphismOfRealSemisimpleLieAlgebras#
#
#  corelg.getDirectSumOfPureLA 
#  corelg.getPureLA
#  corelg.realification 
#  PositiveRootsNF
#  BilinearFormMatNF
#  PositiveRootsAsWeights
#  SignatureTable
#  corelg.ConjugationFct
#  corelg.SOSets  
#  corelg.conj_func
#  corelg.so_sets
#  corelg.signs
#  corelg.signsandperm
#  corelg.Sub3
#  corelg.RealFormsOfSimpleLieAlgebra
#  corelg.MakeSqrtFieldCopyOfLieAlgebra
#  corelg.NonCompactRealFormsOfSimpleLieAlgebra
#  corelg.ParametersOfNonCompactRealForm
#  corelg.RealFormByInnerInvolutiveAutomorphism
#  corelg.makeCanGenByBase
#  corelg.enumOfBase
#  corelg.VoganDiagramOfRealForm
#  corelg.VoganDiagramRealification
#  corelg.getRootsystem
#  corelg.SingleVoganDiagram
#  corelg.makeBlockDiagMat
#  corelg.computeIdRealForm
#  corelg.splitRealFormOfSL
#  corelg.prntdg


############################################################################
###########################################################################
#
# first a few functions from QuaGroup, SLA and the library:
# 


# From QuaGroup
InstallMethod( PositiveRootsNF,
        "for a root system",
        true, [ IsRootSystem ], 0,
        function( R )

    local b, st;

    st:= SimpleSystem(R);
    b:= Basis( VectorSpace( DefaultFieldOfMatrix(st), st ), st );
    return List( PositiveRoots(R), x -> Coefficients( b, x ) );
end );


# From QuaGroup
InstallMethod( BilinearFormMatNF,
        "for a root system",
        true, [ IsRootSystem ], 0,
        function( R )

    local m;

    m:= Minimum( List([1..Length(CartanMatrix(R))], i -> 
            BilinearFormMat(R)[i][i] ) );
    return BilinearFormMat(R)*(2/m);
end );


# from GAP library:
InstallMethod( PositiveRootsAsWeights,
    "for a root system",
    true, [ IsRootSystem ], 0,
    function( R )

      local posR,V,lcombs;

      posR:= PositiveRoots( R );
      V:= VectorSpace( DefaultFieldOfMatrix(SimpleSystem(R) ), SimpleSystem( R ) );
      lcombs:= List( posR, r ->
                       Coefficients( Basis( V, SimpleSystem(R) ), r ) );
      return List( lcombs, c -> LinearCombination( CartanMatrix(R), c ) );

end );


# From SLA:
InstallMethod( SignatureTable,
"for Lie algebra", true, [IsLieAlgebra], 0,
function( L )

    local o, R, p, tab, x, w, max, dims, r, u, wt, dc, char, it, i, 
          Ci, h, ev, pos, tp, res, en;
    
    o:= NilpotentOrbits(L);
    R:= RootSystem(L);

    tp:= CartanType( CartanMatrix(R) ).types[1];
    if tp[1] in [ "A", "B", "C", "E", "F", "G" ] then

       p:= PositiveRootsNF(R);
       tab:= [ ];
       for x in o do
           w:= WeightedDynkinDiagram(x);
           max:= p[Length(p)]*w;
           if not IsInt( max ) then # hack to make it work with SqrtField...
              max:= max![1][1][1];
           fi;
           dims:= List([1..max+1], u -> 0 );
           for r in p do
               u:= r*w+1;
               if not IsInt( u ) then
                  u:= u![1][1][1];
               fi;
               dims[u]:= dims[u]+1;
           od;
           dims[1]:= 2*dims[1]+Length(CartanMatrix(R));
           Add( tab, [ dims, w ] );
       od;

       return rec( tipo:= "notD", tab:= tab );

    else

       en:= CartanType( CartanMatrix(R) ).enumeration[1];
       wt:= List( CartanMatrix( R ), x -> 0 );
       wt[en[1]]:= 1;
       dc:= DominantCharacter( L, wt );
       char:= [[],[]];
       for i in [1..Length(dc[1])] do
           it:= WeylOrbitIterator( WeylGroup(R), dc[1][i] );
           while not IsDoneIterator( it ) do
              Add( char[1], NextIterator( it ) );
              Add( char[2], dc[2][i] );
           od;
       od;

       Ci:= FamilyObj(o[1])!.invCM;
       tab:= [ ];
       for x in o do

           h:= Ci*WeightedDynkinDiagram(x);

           dims:= [ ];
           for i in [1..Length(char[1])] do
               ev:= h*char[1][i];

               pos:= PositionProperty( dims, y -> y[1]=ev);
               if pos = fail then
                  Add( dims, [ev, char[2][i]] );
               else
                  dims[pos][2]:= dims[pos][2]+char[2][i];
               fi;
           od;
           Sort( dims, function(a,b) return a[1] < b[1]; end );

           Add( tab, [dims, WeightedDynkinDiagram(x)] );
       od;

       res:= rec( tipo:= "D", char1:= char, tab1:= tab, V1:= 
                HighestWeightModule( L, wt ) );

       wt:= List( CartanMatrix( R ), x -> 0 );
       wt[en[Length(wt)]]:= 1;
       dc:= DominantCharacter( L, wt );
       char:= [[],[]];
       for i in [1..Length(dc[1])] do
           it:= WeylOrbitIterator( WeylGroup(R), dc[1][i] );
           while not IsDoneIterator( it ) do
              Add( char[1], NextIterator( it ) );
              Add( char[2], dc[2][i] );
           od;
       od;

       Ci:= FamilyObj(o[1])!.invCM;
       tab:= [ ];
       for x in o do

           h:= Ci*WeightedDynkinDiagram(x);

           dims:= [ ];
           for i in [1..Length(char[1])] do
               ev:= h*char[1][i];

               pos:= PositionProperty( dims, y -> y[1]=ev);
               if pos = fail then
                  Add( dims, [ev, char[2][i]] );
               else
                  dims[pos][2]:= dims[pos][2]+char[2][i];
               fi;
           od;
           Sort( dims, function(a,b) return a[1] < b[1]; end );

           Add( tab, [dims, WeightedDynkinDiagram(x)] );
       od;

       res.char2:= char; res.tab2:= tab;
       res.V2:= HighestWeightModule( L, wt );
       return res;
    fi;   

end );



#############################
#
# computes the realification of a simple complex LA over F
#
corelg.realification := function(arg)
local sc, scn, i, j, k, F, type, rank, L, cg, cb, rs, bas, dim, rts, cbn, n, prs, nrs,posp,posn,
      l1,l2,l3,en,Ln, basn, l, csa, K, P, bascd, theta,cd, cgn,tmp,v, sp, R, CartInt,allrts, fundr;

   type := arg[1];
   rank := arg[2];
   F    := SqrtField;
   if Length(arg)=3 then F:=arg[3]; fi;
  
 ##this is complex simple LA and its data
   L    := SimpleLieAlgebra(type,rank,GaussianRationals);
   rs   := RootSystem(L);;
   cg   := CanonicalGenerators(rs);;
   cb   := ChevalleyBasis(L);; ##changed this!
   bas  := Basis(L);
   sc   := StructureConstantsTable(bas);;
   dim  := Dimension(L);
 
 ##now create structure constants of realification of L
 ##take as basis the elements of bas and \imath*bas
 ##
   scn  := EmptySCTable( 2*dim,  Zero(F),  "antisymmetric" );;
   for i in [1..dim-1] do
      SetEntrySCTable( scn, i, i+dim, []);
      for j in [i+1..dim] do
         en := sc[i][j];
         l1 := [];
         l2 := [];
         l3 := [];
         for k in [1..Length(en[1])] do
            Add(l1,en[2][k]*One(F));  Add(l1,en[1][k]);     
            Add(l2,en[2][k]*One(F));  Add(l2,en[1][k]+dim);
            Add(l3,-en[2][k]*One(F)); Add(l3,en[1][k]);
         od;
         SetEntrySCTable( scn, i, j, l1 );        ## prod of two old basis vecs
         SetEntrySCTable( scn, i, j+dim, l2);     ## prod of old + new basis 
         SetEntrySCTable( scn, i+dim,j , l2);     ## prod of old + new basis
         SetEntrySCTable( scn, i+dim, j+dim, l3); ## prod of two new basis
      od;
   od;
   SetEntrySCTable( scn, dim, 2*dim, []);
   
 ##now construct realification and set CSA
 ##
   Ln   := LieAlgebraByStructureConstants(F,scn);
   basn := Basis(Ln);
   csa  := basn{Concatenation([dim-rank+1..dim],[2*dim-rank+1..2*dim])};
   csa  := SubalgebraNC(Ln,csa);
   SetCartanSubalgebra(Ln,csa);

   SetMaximallyCompactCartanSubalgebra(Ln,csa);
 
 ##Cartan decomposition: K is compact real form of L (onishik p 26)
 ##that is, spanned by \imath*h_i, x_a - x_{-a}, \imath*(x_a+x_{-a})
 ##consequently, P is spanned by \imath times these elts
 ##
   K    := basn{[2*dim-rank+1..2*dim]};  ## ih_1,..,ih_rank
   l    := Length(cb[1]);
   for i in [1..l] do
      Add(K, basn[i]-basn[l+i]);         ## (x_a-x_{-a})
      Add(K, basn[dim+i]+basn[dim+l+i]); ## i(x_a+x_{-a})
   od;
   K    := SubalgebraNC(Ln,K,"basis");
   P    := basn{[dim-rank+1..dim]};      ## h_1,..,h_rank
   for i in [1..l] do
      Add(P, basn[dim+i]-basn[dim+l+i]); ## i(x_a-x_{-a})
      Add(P, basn[i]+basn[l+i]);         ## (x_a+x_{-a})
   od;
   P     := Subspace(Ln,P,"basis");
   SetCartanSubalgebra(K,SubalgebraNC(K,basn{[2*dim-rank+1..2*dim]}));

 ##set Cartan decomposition; 
 ##create corresponding Cartan involution
   bascd := BasisNC(Ln,Concatenation(Basis(K),Basis(P)));
   theta := function(v)
      local k, p, cf, i;   
      k   := Length(Basis(K));
      p   := Length(Basis(P));
      cf  := List(Coefficients(bascd,v),x->x);
      for i in [k+1..k+p] do cf[i] := -cf[i]; od;
      return cf*bascd;
   end;
   
   SetCartanDecomposition(Ln,rec( K:= K, P:= P, CartanInv :=theta));

  ##new chevalley basis
   cbn := [[],[],[]];
   l   := Length(cb[1]);
   for i in [1..l] do
      Add(cbn[1], 1/2*One(F)*( basn[i]+E(4)*One(F)*basn[i+dim]   ) );
      Add(cbn[1], 1/2*One(F)*( basn[i]-E(4)*One(F)*basn[i+dim]   ) );
      Add(cbn[2],  1/2*One(F)*( basn[i+l]+E(4)*One(F)*basn[i+dim+l]   ) );
      Add(cbn[2],  1/2*One(F)*( basn[i+l]-E(4)*One(F)*basn[i+dim+l]   ) );
      if i <= rank then
         Add(cbn[3], (1/2)*One(F)*(basn[2*l+i]+E(4)*One(F)*basn[2*l+i+dim]));
         Add(cbn[3],  (1/2)*One(F)*(basn[2*l+i]-E(4)*One(F)*basn[2*l+i+dim]));
      fi;
   od;

    n   := 2*rank;
    rts := [ ]; 
    for v in cbn[1] do 
        sp:= BasisNC(SubspaceNC(Ln,[v],"basis"),[v]); 
        Add( rts, List( cbn[3], t -> Coefficients(sp,t*v)[1] ) ); 
    od;

    R:= Objectify( NewType( NewFamily( "RootSystemFam", IsObject ),
                IsAttributeStoringRep and IsRootSystemFromLieAlgebra ), 
                rec() );
    SetCanonicalGenerators( R, [ cbn[1]{[1..n]}, cbn[2]{[1..n]}, cbn[3] ] );
    SetUnderlyingLieAlgebra( R, Ln );
    SetPositiveRootVectors( R, cbn[1] );
    SetNegativeRootVectors( R, cbn[2] );

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
    if IsSqrtField(F) then
       rts := List(rts, x-> List(x, SqrtFieldEltToCyclotomic));
    fi;

    SetPositiveRoots( R, rts );
    SetNegativeRoots( R, -rts ); 
    SetSimpleSystem( R, rts{[1..n]} );

    SetRootSystem(L,R);
    SetChevalleyBasis(R,cbn);
    SetRootSystem(MaximallyCompactCartanSubalgebra(Ln),R);
    SetRootSystem(CartanSubalgebra(Ln),R);
    SetChevalleyBasis( Ln, cbn );
 
    return Ln;
end;




########################################################################

corelg.signs:= function( type, n )

     local sgn, i, m, s;

     sgn:= [ ];
     if type ="A" then
        if IsEvenInt(n) then
           m:= n/2;
        else
           m:= (n+1)/2;
        fi;
        for i in [1..m] do
            s:= List( [1..n], x -> 1 );
            s[i]:= -1;
            Add( sgn, s );
        od;
     elif type = "B" then
        for i in [1..n] do
            s:= List( [1..n], x -> 1 );
            s[i]:= -1;
            Add( sgn, s );
        od;
     elif type = "C" then
        if IsEvenInt(n-1) then
           m:= (n-1)/2;
        else
           m:= n/2;
        fi;
        for i in [1..m] do
            s:= List( [1..n], x -> 1 );
            s[i]:= -1;
            Add( sgn, s );
        od;
        s:= List( [1..n], x -> 1 ); s[n]:= -1;
        Add( sgn, s );
     elif type = "D" then
        if n = 4 then
           return [ [ 1, 1, -1, 1 ], [ 1, -1, 1, 1 ] ];
        fi;
        if IsEvenInt(n-3) then
           m:= 1+(n-3)/2;
        else
           m:= 1+(n-2)/2;
        fi;
        for i in [1..m] do
            s:= List( [1..n], x -> 1 );
            s[i]:= -1;
            Add( sgn, s );
        od;
        s:= List( [1..n], x -> 1 ); s[n-1]:= -1;
        Add( sgn, s );
     elif type = "E" then
        if n = 6 then
           sgn:= [ [ -1, 1, 1, 1, 1, 1 ], [ 1, -1, 1, 1, 1, 1 ] ];
        elif n = 7 then 
           sgn:= [ [ -1, 1, 1, 1, 1, 1, 1 ], [ 1, -1, 1, 1, 1, 1, 1 ], [ 1, 1, 1, 1, 1, 1, -1 ] ];
        elif n = 8 then
           sgn:= [ [ -1, 1, 1, 1, 1, 1, 1, 1 ], [ 1, 1, 1, 1, 1, 1, 1, -1 ] ];
        fi;
     elif type = "F" then
        sgn:= [ [ 1, 1, 1, -1 ], [ 1, 1, -1, 1 ] ];
     else
        sgn:= [ [ 1, -1 ] ];
     fi;

     return sgn;

end;




########################################################################
corelg.signsandperm:= function( type, n )

   local p, sgn, s, i, m;

    if type = "A" then
       p:= PermList( [n,n-1..1] );
       if IsEvenInt(n) then
          # there is only one...
          sgn:= [ List( [1..n], x -> 1 ) ];
       elif n > 1 then
          sgn:= [ List( [1..n], x -> 1 ), List( [1..n], x -> 1 ) ];
          sgn[2][ (n+1)/2 ]:= -1;
       fi;
    elif type = "D" then
       p:= (n-1,n);
       sgn:= [ List( [1..n], x -> 1 ) ];
       if IsEvenInt(n) then
          m:= n/2-1;
       else
          m:= (n-1)/2;
       fi;
       for i in [1..m] do
           s:= List( [1..n], x -> 1 );
           s[i]:= -1;
           Add( sgn, s );
       od;
    elif type ="E" and n=6 then
       p:= (1,6)*(3,5);
       sgn:= [ [1,1,1,1,1,1], [1,1,1,-1,1,1] ];
    else
       Error("no outer auts");
    fi;

    return rec( sg:= sgn, perm:= p );

end;





########################################################################
corelg.Sub3:=function( arg )

local L, R, P, S, T, s, p, TT, i, j, w, g, F, makeCartInv, a, n, F0, bb, BB, KK1, KK2, rts, v, sp, 
      CartInt, allrts, fundr; 

a:= arg[1];
n:= arg[2];
if Length(arg)=3 then
   F0:= arg[3];
else
   F0:= GaussianRationals;
fi;

#R:=[]; P:=[]; S:=[]; C:=[]; T:=[]; TT:=[]; D:=[]; U:=[]; c:=[]; w:=[];  
L:= SimpleLieAlgebra( a,  n,  Rationals);
R:= RootSystem(L);
P:= PositiveRoots(R);
S:= SimpleSystem(R);
#C:= ChevalleyBasis(L);
#V:= VectorSpace(Rationals,  S);
#B:= Basis(V, S);
T:= StructureConstantsTable(Basis(L));;
s:= Length(S);;
p:= Length(P);;
#D:=[];;
#U:=[];;
TT:=EmptySCTable( 2*p+s,  Zero(F0),  "antisymmetric" );;


# Cerchiamo ora di assegnare i valori dei bracket nella tabella moltiplicativa
#Quelli fra i generatori H sono nulli,  quindi non devo fare niente,  è automatico in TT


#Cerco di sistemare i bracket fra i generatori [H, X] e [H.Y]
 
for i in [1..s] do
	for j in [1..p] do
		if not IsEmpty(T[2*p+i][j][2]) then
SetEntrySCTable( TT,  2*p+i,    j, Flat( [ T[2*p+ i][j][2][1]  ,  p+j    ] )  );
		fi;
		if not IsEmpty(T[2*p+i][p+j][2]) then
SetEntrySCTable( TT,  2*p+i,  p+j, Flat( [ T[2*p+ i][p+j][2][1],    j    ] )  );
		fi;
	od;
od;



#setto i prodotti [X,X] [Y,Y] con indici diversi 

for i in [1..p] do
	for j in [1..i-1] do
		if P[i]+P[j] in P then
			if P[i]-P[j] in P then  #ricordo che i>j in questo caso
SetEntrySCTable( TT, i, j, Flat([ T[i][j][2], Position(P, P[i]+P[j]),  T[p+i][j][2],  Position(P, P[i]-P[j] )   ]));
SetEntrySCTable( TT, p+i, p+j,Flat([ T[p+i][p+j][2], Position(P, P[i]+P[j]), T[p+i][j][2], Position(P, P[i]-P[j] )  ]));
			else #i-j non fa radice
SetEntrySCTable( TT, i, j, Flat([ T[i][j][2], Position(P, P[i]+P[j] )  ]));
SetEntrySCTable(TT, p+i, p+j, Flat([ T[p+i][p+j][2], Position(P, P[i]+P[j] )   ]));
			fi;
		else #i+j non radice
			if P[i]-P[j] in P then
SetEntrySCTable( TT, i, j, Flat([ T[p+i][j][2],  Position(P, P[i]-P[j] )   ]));
SetEntrySCTable( TT, p+i, p+j, Flat([ T[p+i][j][2], Position(P, P[i]-P[j] )   ]));
			fi;
		fi;
	od;

	for j in [i+1..p] do
		if P[i]+P[j] in P then
			if P[j]-P[i] in P then
SetEntrySCTable( TT, i, j, Flat([ T[i][j][2], Position(P, P[i]+P[j]),  T[i][p+j][2],  Position(P, P[j]-P[i] )   ]));
SetEntrySCTable(TT, p+i, p+j, Flat([ T[p+i][p+j][2], Position(P, P[i]+P[j]), T[i][p+j][2], Position(P, P[j]-P[i] )  ]));
 			else #j-i non fa radice
SetEntrySCTable( TT, i, j, Flat([ T[i][j][2], Position(P, P[j]+P[i] )  ]));
SetEntrySCTable(TT, p+i, p+j, Flat([ T[p+i][p+j][2], Position(P, P[j]+P[i] )   ]));
			fi;
		else
			if P[j]-P[i] in P then
SetEntrySCTable( TT, i, j, Flat([ T[i][p+j][2],  Position(P, P[j]-P[i] )   ]));
SetEntrySCTable(TT, p+i, p+j, Flat([ T[i][p+j][2], Position(P, P[j]-P[i] )   ]));
			fi;
		fi;
	od;

od;

#metto a posto i generatori [X, Y] con lo stesso indice

for i in [1..p] do
	g:=T[i][p+i];
	w:=[];
	for j in [1..Length(g[2])] do
		Add( w, 2*g[2][j]); 
		Add( w, g[1][j] );
	od;
SetEntrySCTable( TT,  i, p+i,  w);
od;


#cerco di sistemare i prodotti [X,Y] con indici differenti


for i in [1..p] do
	for j in [1..i-1] do
		if  P[i]+P[j] in P then
			if P[i]-P[j] in P  then
SetEntrySCTable( TT, i, p+j, Flat([ T[i][j][2], p+Position(P, P[i]+P[j] ), T[i][p+j][2], p+Position(P, P[i]-P[j] )  ]));
			else
SetEntrySCTable( TT, i, p+j, Flat([ T[i][j][2], p+Position(P, P[i]+P[j] )   ]));
			fi;
		else
			if P[i]-P[j] in P then
SetEntrySCTable( TT, i, p+j, Flat([  T[i][p+j][2],  p+Position(P, P[i]-P[j] )   ]));
			fi;
		fi;
	od;

	for j in [i+1..p] do
		if  P[i]+P[j] in P then
			if P[j]-P[i] in P  then
SetEntrySCTable( TT, i, p+j, Flat([ T[i][j][2], p+Position(P, P[i]+P[j] ), T[i][p+j][2], p+Position(P, P[j]-P[i] )  ]));
			else
SetEntrySCTable( TT, i, p+j, Flat([ T[i][j][2],  p+Position(P, P[i]+P[j] )   ]));
			fi;
		else
			if P[j]-P[i] in P then
SetEntrySCTable( TT, i, p+j, Flat([ T[i][p+j][2],  p+Position(P, P[j]-P[i] )   ]));
			fi;
		fi;
	od;
od;


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
 

  L:=LieAlgebraByStructureConstants(F0, TT);
  SetCartanDecomposition( L, rec( K:= L, P:= SubspaceNC( L, [ ],"basis" ),
                          CartanInv := makeCartInv(L,L,SubspaceNC(L,[],"basis"))));
  SetIsCompactForm( L, true );

  # fare un elenco con tre elenchi: [x_alpha], [x_{-alpha}], [h_1,...,h_l], alpha > 0, tutti
  # i vettori espressi in termini della base di L.
  # Se b:= Basis(L); allora b[1] è il primo elemento della base, ecc. 

  bb:=Basis(L);
  BB:=[[],[],[]];
  KK1:=0;
  KK2:=0;

  for j in [1..s] do
      BB[3][j]:=(-1*E(4)*One(F0)*bb[2*p+j]);
  od;

  for i in [1..p] do

      KK1:=(bb[i]); #X_alpha
      KK2:=(bb[p+i]); #Y_alpha
      BB[1][i]:= 1/2*One(F0)*(KK1-1*E(4)*One(F0)*KK2);
      BB[2][i]:= 1/2*One(F0)*(-1*One(F0)*KK1-1*E(4)*One(F0)*KK2);

  od;

  SetCartanSubalgebra(L,Subalgebra(L,BB[3]) );
  SetMaximallyCompactCartanSubalgebra( L, CartanSubalgebra(L) );

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
 

    SetPositiveRoots( R, rts );
    SetNegativeRoots( R, -rts ); 
    SetSimpleSystem( R, rts{[1..n]} );

    SetRootSystem(L,R);
    SetRootSystem(MaximallyCompactCartanSubalgebra(L),R);
    SetRootSystem(CartanSubalgebra(L),R);  ###!!! added this recently

    SetChevalleyBasis( L, BB );


  return L;

end;



##############################################################################
##
##  returns all real forms of simple Lie algebras of type <type> and rank <n>
##  up to isomorphism
##
corelg.RealFormsOfSimpleLieAlgebra := function( arg )

  local forms, s, i, tmp, type, n, F;

  type := arg[1];   
  n    := arg[2];
  if Length(arg)=3 then F:=arg[3]; else F:=GaussianRationals; fi;

  forms:= [ corelg.Sub3( type, n, F ) ]; # so the compact form...
  SetIsRealFormOfInnerType(forms[1],true);
  SetRealFormParameters(forms[1],[type,n,ListWithIdenticalEntries(n,1),()]);
  
  s:= corelg.signs( type, n );
  for i in [1..Length(s)] do
      tmp := corelg.SuperLie( type, n, s[i], (), F );
      SetRealFormParameters(tmp, [type,n,s[i],()]);
      SetIsRealFormOfInnerType(tmp,true);
      Add( forms, tmp );
  od;

  if type in ["A","D","E"] and (type <> "E" or n = 6) and not (type = "A" and n = 1) then

     s:= corelg.signsandperm( type, n ); 
     for i in [1..Length(s.sg)] do
         tmp := corelg.SuperLie( type, n, s.sg[i], s.perm, F );
         SetRealFormParameters(tmp, [type,n,s.sg[i],s.perm]);
         SetIsRealFormOfInnerType(tmp,false);
         Add( forms, tmp );
     od;

  fi;
  
  return forms;
end;


########################################################################

InstallOtherMethod( CartanSubspace,
   "for a Lie algebra with Cartan decomposition",
   true, [ IsLieAlgebra ], 0,
   function( L )

   # L = K + P, note that P does not have nilpotent elements, as a nilpotent
   # e would lie in a hom sl_2 triple, with h\in K, not possible. So a subspace C
   # is a Cartan subspace iff its centralizer in P is equal to C.

   local P, found, b, V, C, k;

   P:= CartanDecomposition(L).P;
   # first we determine the rank by computing any Cartan subspace...

   found:= false;
   b:= ShallowCopy( BasisVectors( Basis( Intersection( P, CartanSubalgebra(L) ) ) ) );
   # first try with just basis elements...
   V:= SubspaceNC( P, b );
   C:= Filtered( Basis(P), x -> ForAll( b, y -> IsZero(x*y) ) and not x in V ); 
   while Length(C) > 0 do
      Add( b, C[1] );
      V:= SubspaceNC( P, b );
      C:= Filtered( C, x -> ForAll( b, y -> IsZero(x*y) ) and not x in V ); 
   od;

   if Dimension( Intersection( LieCentralizer( L, V ), P ) ) = Length(b) then
      return V;
   fi;

   b:= ShallowCopy( BasisVectors( Basis( Intersection( P, CartanSubalgebra(L) ) ) ) );
   V:= SubspaceNC( P, b );
   C:= Intersection( LieCentralizer( L, V ), P );
   while not found do
      k:= 1;
      while k <= Dimension(C) do
         if not Basis(C)[k] in V then
            Add( b, Basis(C)[k] );
            break;
         else
            k:= k+1;
         fi;
      od;
      C:= Intersection( LieCentralizer( L, SubalgebraNC( L, b ) ), P );
      if Dimension(C) = Length(b) then
         found:= true;
      else
         V:= SubspaceNC( P, b );
      fi;
   od;

   return C;

end );


########################################################################

corelg.MakeSqrtFieldCopyOfLieAlgebra := function(L)
local MSF, RSF, writeToSF, T, R, rank, csa, ct, cb, ci, K, P, k, p, v, vnew, tmp, mkWhere,TT, i,j, bas;

   ct := Runtime();
   T  := StructureConstantsTable(Basis(L)); 
   if not ForAll(Flat(T),IsRat) then
      Error("SCTable not rational");
   fi; 
   
   TT := ShallowCopy(T);
   for i in [1..Length(TT)] do 
      if IsList(TT[i]) then
         TT[i] := ShallowCopy(TT[i]); 
         for j in [1..Length(TT[i])] do
            TT[i][j] := ShallowCopy(TT[i][j]);
            TT[i][j][2] := TT[i][j][2]*One(SqrtField);
         od; 
      fi;
   od;
   TT[Length(TT)] := Zero(SqrtField);
   T := TT;

   if not HasRootSystem(L) then
      Error("Liealg has no rootsystem attached");
   fi;
   R    := RootSystem(L);
   MSF  := LieAlgebraByStructureConstants( SqrtField, T);
   writeToSF := function(v)
   local er;
     #er := ExtRepOfObj(v)*Sqroot(1);
      er := List(ExtRepOfObj(v),SqrtFieldEltByCyclotomic);
      return ObjByExtRep(FamilyObj(Zero(MSF)),er);
   end;            
   csa := BasisVectors(Basis(CartanSubalgebra(L)));
   csa := List(csa,writeToSF);
   csa := SubalgebraNC(MSF, csa,"basis");
   SetCartanSubalgebra(MSF, csa);

   RSF := Objectify( NewType( NewFamily( "RootSystemFam", IsObject ),
               IsAttributeStoringRep and IsRootSystemFromLieAlgebra ), 
               rec() );
   SetCanonicalGenerators( RSF, List(CanonicalGenerators(R),x->List(x,writeToSF)));
   SetUnderlyingLieAlgebra( RSF, MSF );
   SetPositiveRootVectors( RSF, List(PositiveRootVectors(R),writeToSF));
   SetNegativeRootVectors( RSF, List(NegativeRootVectors(R),writeToSF));
   SetCartanMatrix( RSF,  CartanMatrix(R) );
   SetPositiveRoots( RSF, PositiveRoots(R));
   SetNegativeRoots( RSF, NegativeRoots(R));
   SetSimpleSystem( RSF, SimpleSystem(R));
   SetRootSystem(MSF,RSF);
   SetChevalleyBasis(MSF,List(ChevalleyBasis(L),x->List(x,writeToSF)));
   K  := CartanDecomposition(L).K;
   P  := CartanDecomposition(L).P;
   ci := CartanDecomposition(L).CartanInv;
   K  := SubalgebraNC(MSF,List(Basis(K),writeToSF), "basis");
   SetCartanSubalgebra(K,SubalgebraNC(K,
             List(Basis(CartanSubalgebra(CartanDecomposition(L).K)),writeToSF)));
   P  := SubspaceNC(MSF, List(Basis(P),writeToSF), "basis");
   if HasRealFormParameters(L) then SetRealFormParameters(MSF,RealFormParameters(L)); fi;
   bas := BasisNC(MSF,Concatenation(Basis(K),Basis(P)));
   ci  := function(v)
         local k, p, cf, i;
             k   := Length(Basis(K));
             p   := Length(Basis(P));
             cf  := List(Coefficients(bas,v),x->x);
             for i in [k+1..k+p] do cf[i] := -cf[i]; od;
             return cf*bas;
          end;
   SetCartanDecomposition(MSF, rec(K:=K, P:=P, CartanInv:=ci));


   mkWhere := function(signs,mv)
   local i, new;
      new :=[];
      for i in [1..Length(signs)] do
         if signs[i]=-1 then Add(new,"P");
         elif i in Flat(mv) then Add(new,"?");
         else Add(new,"K");
         fi;
     od;
     return new;
   end;

   if HasVoganDiagram(L) then
      v   := VoganDiagram(L);
      tmp := corelg.VoganDiagramOfRealForm(MSF,
                  rec(cg      := List(CanonicalGenerators(v),x->List(x,writeToSF)), 
                      base    := ShallowCopy(BasisOfSimpleRoots(v)), 
                      mv      := ShallowCopy(MovedPoints(v)),
                      signs   := mkWhere(Signs(v),MovedPoints(v)),
                      cfsigma := ShallowCopy(CoefficientsOfSigmaAndTheta(v).cfsigma), 
                      cftheta := ShallowCopy(CoefficientsOfSigmaAndTheta(v).cftheta)));
     #SetPermInvolution(tmp,PermInvolution(v));
      SetVoganDiagram(MSF,tmp);
   fi;

   return rec(liealg := MSF, writeToSF := writeToSF);
end;
######################################################################




##############################################################################
##
##  returns all lists [<type>,<n>, signs, perm] parametrising the real forms
##  of simple Lie algebras of type <type> and rank <n> up to isomorphism
## 
corelg.ParametersOfNonCompactRealForm := function(type,n)
local params, s, i;
  params := [];
  s      := corelg.signs( type, n );
  for i in [1..Length(s)] do
      Add(params, [type,n,s[i],()]);
  od;
  if type in ["A","D","E"] and (type <> "E" or n = 6) then
     s := corelg.signsandperm( type, n ); 
     for i in [1..Length(s.sg)] do
         Add(params,[type,n,s.sg[i],s.perm]);
     od;
  fi;
  return params;
end;



##############################################################################
##  ONLY USED FOR NILPOTENT ORBITS (RECONSTRUCTION OF DATABASE)
##  returns all noncompact real forms of simple Lie algebras of type <type> 
##  and rank <n> up to isomorphism. If <params> is given, then it has to be
##  a sublist of corelgParametersOfNonCompactRealForm( <type>, <n> ); in this case
##  only the real forms parametrised by these entries are constructed.
##  The output is a list with the following entries:
##            liealg    : the real form defined over Gaussian Rationals, 
##            liealgSF  : the real form defined over SqrtField,
##            writeToSF : function from liealg to liealgSF,
##            rank      : <n>,
##            type      : <type>,
##  all Lie algebras have a rootsystem, CartanSubalgebra and CartanDecompositon
##  attached.
##
corelg.NonCompactRealFormsOfSimpleLieAlgebra := function(arg)
local type, n, rforms, L, LSF, forms, tmp, F, withField, i, sigma;

   withField := false;
   if IsField(arg[Length(arg)]) then 
      F := arg[Length(arg)]; 
      withField := true;
      arg := arg{[1..Length(arg)-1]};
   else 
      F := GaussianRationals; 
   fi;
   if Length(arg) = 2 then   
      type := arg[1];
      n    := arg[2];
      rforms := corelg.RealFormsOfSimpleLieAlgebra( type, n, F);
      rforms := rforms{[2..Length(rforms)]};
   elif Length(arg)=1 then
      arg  := arg[1];
      type := arg[1];
      n    := arg[2];
   fi;
   if Length(arg) = 4 then
      tmp := corelg.SuperLie( type, n, arg[3], arg[4],F );
      SetRealFormParameters(tmp, [type,n,arg[3],arg[4]]);
      rforms := [tmp];
   fi;

   if withField then 
      if E(4) in F or IsSqrtField(F) then      
         for i in rforms do sigma := RealStructure(i); od;
      fi;
      if Length(rforms)=1 then return rforms[1]; else return rforms; fi; 
   fi;

   forms := [];
   for L in rforms do
     #Print("now make copies...\n");
      LSF := corelg.MakeSqrtFieldCopyOfLieAlgebra(L);
      SetIsCompactForm(L,false);
      SetIsCompactForm(LSF.liealg,false);
      if RealFormParameters(LSF.liealg)[4]=() then
         SetIsRealFormOfInnerType(LSF.liealg,true);
         SetIsRealFormOfInnerType(L,true);
      else
         SetIsRealFormOfInnerType(LSF.liealg,false);
         SetIsRealFormOfInnerType(L,false);
      fi;
      sigma := RealStructure(LSF);
      Add(forms, rec( liealg    := L, 
                      liealgSF  := LSF.liealg, 
                      writeToSF := LSF.writeToSF,
                      rank      := n,
                      type      := type));
   od;

   if Length(arg)=4 then return forms[1]; else return forms; fi;
end;


##################################################################


##################################################################
# input:  output of FiniteOrderInnerAutomorphism(type,rank,2)
#         (assumes that theta is in std form, that is, it maps
#          (h_i,x_i,y_i) to (h_i,\mu_i x_i, \mu_i^{-1} y_i)
# output: real form with attached CSA, RS and CartanDecomposition;
#         defined over GaussianRationals wrt theta^tau
#
corelg.RealFormByInnerInvolutiveAutomorphism := function(theta)
local makeCartInv, L, ch, i, k0, p0, k, K, P, bas, T, M, R, im, cg,F;

   if IsList(theta) then
      if Length(theta)=4 then F:=theta[4]; else F:=GaussianRationals; fi;
      L  := SimpleLieAlgebra(theta[1],theta[2],F);
      cg := CanonicalGenerators(RootSystem(L));
      im := [List([1..theta[2]],x-> theta[3][x]*cg[1][x]),
             List([1..theta[2]],x-> theta[3][x]*cg[2][x]),
             List([1..theta[2]],x->cg[3][x])];
      theta := LieAlgebraIsomorphismByCanonicalGenerators(L,cg,L,im);
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
 
   L     := Source(theta); 
   F     := LeftActingDomain(L);
   ch    := ChevalleyBasis(L);
   i     := E(4)*One(F);
   k0    := List( ch[3], x -> i*x );
   p0    := [ ];
   for k in [1..Length(ch[1])] do
      if Image( theta, ch[1][k] ) = ch[1][k] then
         Append( k0, [ ch[1][k]-ch[2][k], i*(ch[1][k]+ch[2][k]) ] );
      else
         Append( p0, [ i*(ch[1][k]-ch[2][k]), ch[1][k]+ch[2][k] ] );
      fi;
   od;
   bas := Concatenation( k0, p0 );
   T   := StructureConstantsTable( Basis(L,bas) );
   M   := LieAlgebraByStructureConstants( F , T );
   SetCartanSubalgebra( M, SubalgebraNC( M, Basis(M){[1..Length(ch[3])]}) );
   R   := RootsystemOfCartanSubalgebra(M);
   SetRootSystem(M,R);
   K   := SubalgebraNC(M,Basis(M){[1..Length(k0)]});
   P   := SubspaceNC(M,Basis(M){[Length(k0)+1..Length(bas)]});
   SetCartanDecomposition(M, rec(K:=K, P:=P, CartanInv := makeCartInv(M,K,P)));
   return M;
end;



##################################################################
# input:  positive roots "pr" with corresponding Chev. basis "cb",
#         and a (new) base of simple roots "bas" contained in pr cat -pr
# output: canonical generators wrt base contained in "cb" 
#
corelg.makeCanGenByBase := function(pr,cb,bas)
local tmp, j, pos;
   tmp := [[],[],[]];
   for j in [1..Length(bas)] do
      pos := Position(pr,bas[j]);
      if not pos = fail then
         Add(tmp[1],cb[1][pos]);
         Add(tmp[2],cb[2][pos]);
         Add(tmp[3],cb[1][pos]*cb[2][pos]);
      else
         pos := Position(pr,-bas[j]);
         Add(tmp[1],cb[2][pos]);
         Add(tmp[2],cb[1][pos]);
         Add(tmp[3],cb[2][pos]*cb[1][pos]);
      fi;
   od; 
   return tmp;
end;



##################################################################
# input:  R a root system, base a base:
# output: enumeration wrt can ordering
#
corelg.enumOfBase := function(R,base)
local tmp, bbas, C, en, rank, B;
   rank := Length(SimpleSystem(R));    
   tmp  := BasisNC(VectorSpace(Rationals,IdentityMat(rank)),SimpleSystemNF(R));
   B    := BilinearFormMatNF(R);
   bbas := List(base,x->Coefficients(tmp,x));
   C    := List( bbas, x -> List( bbas, y -> 2*(x*B*y)/(y*B*y) ) );
   en   := Concatenation( CartanType(C).enumeration );
   return en;
end;



##################################################################
# 
# the stuff for Vogan diagrams
#
corelg.VoganDiagramOfRealForm := function(L, list)
      local o, fam, H,R,en,base,C,tmp,signs,i;

      if not IsBound( L!.voganDiagType ) then
         fam:= NewFamily( "vogandiagfam", IsVoganDiagramOfRealForm );
         L!.voganDiagType:= NewType( fam, IsVoganDiagramOfRealForm and IsAttributeStoringRep );
      fi;
     #these just for getting CartanType!
      C :=  corelg.CartanMatrixOfCanonicalGeneratingSet(L,list.cg);
      o := Objectify( L!.voganDiagType, rec(param:=CartanType(C).types) ); #!!!
      SetCanonicalGenerators(o,List(list.cg,x->List(x,y->y)));
      SetBasisOfSimpleRoots(o,list.base);
      SetMovedPoints(o,list.mv);
      tmp := [1..Length(C)];
      for i in list.mv do tmp[i[1]] := i[2]; tmp[i[2]] := i[1]; od;
      SetPermInvolution(o,PermList(tmp));
      signs := [];
      for i in list.signs do if i="P" then Add(signs,-1); else Add(signs,1); fi; od;
      SetSigns(o,signs);    
      
      SetCartanMatrix(o,C);
      
      SetCoefficientsOfSigmaAndTheta(o,rec(cfsigma:=list.cfsigma, cftheta:=list.cftheta));
      return o;
end ;



################################################################################

#this is display:
InstallMethod( PrintObj,
   "for Vogan diagram",
   true,
   [ IsVoganDiagramOfRealForm ], 0,
   function( o )
   local r,t,m,i,minus,signs, tmp;
      r := Sum(List(o!.param,x->x[2]));
      m := MovedPoints(o);
      signs := Signs(o);
      minus := Filtered([1..r],x->Signs(o)[x]=-1);
      corelg.prntdg(CartanMatrix(o),minus);
      Print("\nInvolution: ",PermInvolution(o));
      
      if IsBound(o!.sstypes) then 
         Print("\nTypes of direct summands:\n");
         Print(o!.sstypes); 
      fi;
end );


InstallMethod( ViewObj,
   "for Vogan diagram",
   true,
   [ IsVoganDiagramOfRealForm ], 0,
   function( o )
   local i,tmp,new;
      tmp := List(o!.param,x->Concatenation(x[1],String(x[2])));
      if Length(tmp)>1 then
         new :="";
         for i in [1..Length(tmp)-1] do new := Concatenation([new,tmp[i],"+"]); od;
         tmp := Concatenation(new,tmp[Length(tmp)]);
      else
         tmp := tmp[1];
      fi;
      Print(Concatenation(["<Vogan diagram in Lie algebra of type ",tmp,">"]));
end );



##############################################################################
# input:  realification,
#         defined over Gaussian rationals, with attached Cartan decompositions
#         (record with entries K, P and a function CartanInv)
# output: vogan diagram of Lie algebra
#
corelg.VoganDiagramRealification := function(L)
local res, rank, cd, h, r, c, cg, e, sigma, cf, cb, newcg, iso, wh, phi, i, testCFs, tmpcb,tmpsp,
      getWeyl, liealgs, whs, pos, L1, Lj, isoms, l, isos, j, inn, R, cfs, cft, bb, W, notmv, act,
      tmp,  applyReflection, hs, found, es, fs, h0, vals, pr, posr,  s, B, C, en, base, where,
      sps, mv, cf2, mat, newpr, ims, bbase, bcg, theta, sums, posK, posP, orb, bbas,
      prKind, prK, prP, dim, bas, ct, tt, rr, todo, todoE, todoL, getNewWeyl;
 
      if HasVoganDiagram(L) then return VoganDiagram(L); fi;

      Info(InfoCorelg,2,"   start Vogan Diagram for realification; get CartDecomp and CSA");
      cd    := CartanDecomposition(L);
      h     := MaximallyCompactCartanSubalgebra(L);
      rank  := Dimension(h);
      theta := cd.CartanInv;
      sigma := RealStructure(L);
      R    := RootsystemOfCartanSubalgebra(L,h);
      cb   := ChevalleyBasis(R);
      cg   := CanonicalGenerators(R);
      ct   := CartanType(CartanMatrix(R));
      if not ForAll(Basis(h),x->theta(x) in h) then
         Error("need a theta-stable CSA; Cartan Dec and CSA must be compatible!");
      fi;
      SetIsRealFormOfInnerType(L,false);
      Info(InfoCorelg,2,"   ... done; continue with Vogan Diagram for realification");
      

     #find h0 to define new root ordering; take CSA compatible with h
      tmp := Intersection(cd.K,h);
      hs  := ShallowCopy( CanonicalGenerators(RootsystemOfCartanSubalgebra(cd.K,tmp))[3]);
      
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
      if ct.types[1] = ["F",4] then
           en := en{[4,1,3,2,   8,5,7,6]};
      fi;             
      base := base{en};
       
     #now construct corresponding canonical generators
      newcg := corelg.makeCanGenByBase(pr,cb,base);
      es    := newcg[1];
      fs    := newcg[2];

      sps := List( es, x -> SubspaceNC( L, [x],"basis" ) );
      mv  := [];
      for i in [1..Length(es)] do
         j := PositionProperty( sps, U -> theta( es[i] ) in U );
         if j > i then
            Add(mv,[i,j]);
            es[j]:= theta( es[i] );
            fs[j]:= theta( fs[i] );
         fi;
     od;
     Sort(mv);
     notmv := Filtered([1..rank],x->not x in Flat(mv));
     newcg := [es,fs,List( [1..Length(es)], i -> es[i]*fs[i] ) ];
      
     #computes coefficients wrt sigma and theta 
      testCFs := function(newcg)
      local cft, cfs, i;
         cfs  := ListWithIdenticalEntries(rank,1);
         cft  := ListWithIdenticalEntries(rank,1);
         for i in notmv do
            cfs[i] := Coefficients(Basis(SubspaceNC(L,[newcg[2][i]],"basis"),[newcg[2][i]]),
                                   sigma(newcg[1][i]))[1];
         od;
         for i in mv do
            cfs[i[1]] := Coefficients(Basis(SubspaceNC(L,[newcg[2][i[2]]],"basis"),[newcg[2][i[2]]]), 
                                       sigma(newcg[1][i[1]]))[1] ;
            cfs[i[2]] := Coefficients(Basis(SubspaceNC(L,[newcg[2][i[1]]],"basis"),[newcg[2][i[1]]]), 
                                       sigma(newcg[1][i[2]] ))[1] ;
          od;
         for i in notmv do
            cft[i] := Coefficients(Basis(SubspaceNC(L,[newcg[1][i]],"basis"),[newcg[1][i]]), 
                                    theta(newcg[1][i] ))[1] ;
         od;
         for i in mv do
            cft[i[1]] := Coefficients(Basis(SubspaceNC(L,[newcg[1][i[2]]],"basis"),[newcg[1][i[2]]]),
                                      theta(newcg[1][i[1]]))[1] ;
            cft[i[2]] := Coefficients(Basis(SubspaceNC( L,[newcg[1][i[1]]],"basis"),[newcg[1][i[1]]]), 
                                      theta(newcg[1][i[2]] ))[1] ;
         od;
         if not ForAll([1..rank],x-> cfs[x]*cft[x]<cft[x]*0) then Error("mhmm..signs wrong"); fi; ##^0
         if not ForAll(Flat(mv),x->cft[x]=cft[x]^0) then Error("mhmm"); fi;
         return rec(cfs := cfs, cft := cft);
      end; 


 
      tmp := testCFs(newcg);
      cft := tmp.cft;
      cfs := tmp.cfs;
      tmp := corelg.VoganDiagramOfRealForm(L,
                 rec(cg:=newcg,  
                     base:=base, 
                     mv := mv,
                     signs:=ListWithIdenticalEntries(2*Length(mv),"?"),
                     cfsigma:=cfs, 
                     cftheta:=cft));
      SetVoganDiagram(L,tmp);
     
      if Length(tmp!.param)=1 then
         mv := IdRealForm(L);
         SetRealFormParameters(L,RealFormParameters(RealFormById(mv)));
      fi;
      Info(InfoCorelg,2,"   end Vogan Diagram for realification"); 
      tmp := VoganDiagram(L);
### added this
      if Length(ct.types)=2 and ct.types[1] = ct.types[2] then
         tmp!.sstypes :=  [Concatenation(ct.types[1],[0])];
        #Print("added ",tmp!.sstypes,"\n");
      else
         Display("did NOT add id to realification...");
      fi;
       return tmp;
end;


##############################################################################
# input:  simple real form,
#         defined over Gaussian rationals, with attached Cartan decompositions
#         (record with entries K, P and a function CartanInv)
# output: vogan diagram of Lie algebra
#
corelg.SingleVoganDiagram := function(L)
local res, rank, cd, h, r, c, cg, e, sigma, cf, cb, newcg, iso, wh, phi, i, testCFs, tmpcb,tmpsp,
      getWeyl, liealgs, whs, pos, L1, Lj, isoms, l, isos, j, inn, R, cfs, cft, bb, W, notmv, act,
      tmp,  applyReflection, hs, found, es, fs, h0, vals, pr, posr,  s, B, C, en, base, where,
      sps, mv, cf2, mat, newpr, ims, bbase, bcg, theta, sums, posK, posP, orb, bbas,
      prKind, prK, prP, dim, bas, ct, tt, rr, todo, todoE, todoL, getNewWeyl;
 
      if HasVoganDiagram(L) then return VoganDiagram(L); fi;

      Info(InfoCorelg,2,"   start Vogan Diagram for simple LA; get CartDecomp and CSA");
      cd    := CartanDecomposition(L);
      h     := MaximallyCompactCartanSubalgebra(L);
      rank  := Dimension(h);
      theta := cd.CartanInv;
      sigma := RealStructure(L);
      inn  := rank = Dimension(CartanSubalgebra(cd.K));
      R    := RootsystemOfCartanSubalgebra(L,h);
      cb   := ChevalleyBasis(R);
      cg   := CanonicalGenerators(R);
      ct   := CartanType(CartanMatrix(R));
      if not ForAll(Basis(h),x->theta(x) in h) then
         Error("need a theta-stable CSA; Cartan Dec and CSA must be compatible!");
      fi;
      SetIsRealFormOfInnerType(L,inn);
      Info(InfoCorelg,2,"    ... done; continue with Vogan Diagram for simple LA");

    ##compact form
      if Dimension(cd.P)=0 then
         base  := SimpleSystemNF(R);
         pr    := PositiveRootsNF(R);
         en    := Concatenation( CartanType(CartanMatrix(R)).enumeration );
         base  := base{en};
        #now enumeration is [1...r]; for F4 make it [2,4,3,1]:
         if ct.types[1]=["F",4] then base := base{[4,1,3,2]}; fi;
        #now construct corresponding canonical generators
         cg := corelg.makeCanGenByBase(pr,cb,base);

         cfs  := ListWithIdenticalEntries(rank,1);
         for i in [1..rank] do
            cfs[i] := Coefficients(BasisNC(SubspaceNC(L,[cg[2][i]],"basis"),[cg[2][i]]),
                                                    sigma(cg[1][i]))[1];
         od;  
         tmp  := corelg.VoganDiagramOfRealForm(L,rec(cg   := cg, 
                                                      base := base, 
                                                      mv:=[],
                                                      signs:=ListWithIdenticalEntries(rank,"K"), 
                                                      cfsigma:=cfs, 
                                                      cftheta:=ListWithIdenticalEntries(rank,1)));
         SetVoganDiagram(L,tmp);
         SetRealFormParameters(L,[ct.types[1][1],ct.types[1][2],ListWithIdenticalEntries(ct.types[1][2],1),()]);

         tmp := VoganDiagram(L);
         tmp!.sstypes := [IdRealForm(L)];
         L!.sstypes   := [IdRealForm(L)];
         return tmp;
      fi;
 
      

    ################################
    # SOME PRELIMINARY FUCTIONS    #
    ################################

    #input Cartan decomposition "cd" and can gens "cg"
    #returns list of "K" and "P" wrt cg[i] lying in cd.K or cd.P
      where := function(cd,cg)
      local i, wh; 
         wh := [];     
         for i in cg[1] do
            if i in cd.K then Add(wh,"K"); 
            elif i in cd.P then Add(wh,"P"); 
            else Add(wh,"?"); fi;
         od; 
      return wh;
      end;


     #computes coefficients wrt sigma and theta 
      testCFs := function(newcg)
      local cft, cfs, i;
         cfs  := ListWithIdenticalEntries(rank,1);
         cft  := ListWithIdenticalEntries(rank,1);
         for i in notmv do
            cfs[i] := Coefficients(Basis(SubspaceNC(L,[newcg[2][i]],"basis"),[newcg[2][i]]),
                                   sigma(newcg[1][i]))[1];
         od;
         for i in mv do
            cfs[i[1]] := Coefficients(Basis(SubspaceNC(L,[newcg[2][i[2]]],"basis"),[newcg[2][i[2]]]), 
                                       sigma(newcg[1][i[1]]))[1] ;
            cfs[i[2]] := Coefficients(Basis(SubspaceNC(L,[newcg[2][i[1]]],"basis"),[newcg[2][i[1]]]), 
                                       sigma(newcg[1][i[2]] ))[1] ;
          od;
         for i in notmv do
            cft[i] := Coefficients(Basis(SubspaceNC(L,[newcg[1][i]],"basis"),[newcg[1][i]]), 
                                    theta(newcg[1][i] ))[1] ;
         od;
         for i in mv do
            cft[i[1]] := Coefficients(Basis(SubspaceNC(L,[newcg[1][i[2]]],"basis"),[newcg[1][i[2]]]),
                                      theta(newcg[1][i[1]]))[1] ;
            cft[i[2]] := Coefficients(Basis(SubspaceNC( L,[newcg[1][i[1]]],"basis"),[newcg[1][i[1]]]), 
                                      theta(newcg[1][i[2]] ))[1] ;
         od;
         if not ForAll([1..rank],x-> cfs[x]*cft[x]<cft[x]*0) then Error("mhmm..signs wrong"); fi; ##^0
         if not ForAll(Flat(mv),x->cft[x]=cft[x]^0) then Error("mhmm"); fi;
         return rec(cfs := cfs, cft := cft);
      end; 
     
     #apply the reflection s_{base[j]}\in W to the can gen set newcg
     #return record with new can gens, new base, new Weyl group gens (wrt new base)
      applyReflection := function(newcg,j,base)
      local tmp, pos,ims,W;

        #get Weyl automorphism
          ims := List(base,x-> x-(2*(x*B* base[j])/( base[j]*B* base[j]))* base[j]);           
          W   := LieAlgebraIsomorphismByCanonicalGenerators(L,newcg,L,corelg.makeCanGenByBase(pr,cb,ims));

         newcg := List(newcg,x->List(x,y->Image(W,y)));
         for i in mv do
            #really need this, e.g. if L=RealFormById("E",6,2)
             newcg[1][i[2]] := theta(newcg[1][i[1]]);
             newcg[2][i[2]] := theta(newcg[2][i[1]]);
             newcg[3][i[2]] := newcg[1][i[2]]*newcg[2][i[2]];
         od;
         wh  := where(cd,newcg);         
         tmp := [];
         for i in newcg[1] do 
            pos := PositionProperty(cb[1],x->i in SubspaceNC(L,[x],"basis")); 
            if not pos = fail then
               Add(tmp,pr[pos]);
            else
               pos := PositionProperty(cb[2],x->i in SubspaceNC(L,[x],"basis")); 
               Add(tmp,-pr[pos]);
            fi;
         od;
         
         return rec(cg:=newcg, wh:=wh, base:=tmp);
      end;

 
     ################################################# 

     #consider form of INNER TYPE
 
     if inn then

         mv    :=[];
         notmv := [1..rank];
         base  := SimpleSystemNF(R); 
         pr    := PositiveRootsNF(R);
         en    := Concatenation( CartanType(CartanMatrix(R)).enumeration );
         base  := base{en};

        #now enumeration is [1...r]; for F4 make it [2,4,3,1]:
         if ct.types[1]=["F",4] then 
            base := base{[4,1,3,2]};
         fi;

        #now construct corresponding canonical generators
         newcg := corelg.makeCanGenByBase(pr,cb,base);
         B     := BilinearFormMatNF(R);
  
         wh    := where(cd,newcg);
        #Print("this is 1st wh ",wh,"\n");
         tt    := ct.types[1][1];  
         rr    := ct.types[1][2]; 
         
         if tt="D" and rr=4 then
             pos   := PositionsProperty(wh,x->x="P");
             todo := [];
             if Length(pos)=4 then todo := [2]; fi;
             if Length(pos)=3 then
                tmp := Filtered([1..4],x-> not x in pos)[1];
                if tmp = 2 then todo := [1,2]; else todo := [2,tmp]; fi;
             fi;
             if Length(pos)=2 then
                tmp := Filtered([1..4],x-> not x in pos);
                if 2 in tmp then 
                   i    := Filtered(tmp,x-> not x = 2)[1]; 
                   todo := [pos[1],2,i];
                else
                   todo := Filtered(pos,x->not x=2);
                fi;
             fi;
             for i in todo do
                newcg := applyReflection(newcg,i,base);
                wh    := newcg.wh; 
                base  := newcg.base;
                newcg := newcg.cg;
                pos := PositionsProperty(wh,x->x="P");
               #Print("this is new wh ",wh,"\n");
            od;
            if not Length(pos)=1 then Error("ups"); fi;
            if not pos[1] in [2,3] then
               tmp   := Filtered([1..4],x-> not x in [pos[1],2]);
               tmp   := [tmp[1],2,pos[1],tmp[2]];
               if not IsDuplicateFreeList(tmp) then Error("ups"); fi;
               base  := base{tmp};
               newcg := corelg.makeCanGenByBase(pr,cb,base);
               wh    := where(cd,newcg);
               pos   := PositionsProperty(wh,x->x="P");              
            fi;
            if not Length(pos)=1 or not pos[1] in [2,3] then Error("upsi"); fi;
         fi;

       ###
       # TYPE A and B:
       # find a base such that
       #   A: have at most one "P" in the first \lceil rank/2\rceil entries
       #   B: have at most one "P"
       ###
         if tt in ["A","B"] then
            pos := PositionsProperty(wh,x->x="P");
            while Length(pos)>1 do 
               newcg := applyReflection(newcg,pos[Length(pos)-1],base);
               wh    := newcg.wh; 
               base  := newcg.base;
               newcg := newcg.cg;
               pos := PositionsProperty(wh,x->x="P");
            od;
           #diag aut:
            if tt = "A" and pos[1] > rank/2 then
               base  := Reversed(base);
               newcg := corelg.makeCanGenByBase(pr,cb,base);
               wh    := where(cd,newcg);
               pos := PositionsProperty(wh,x->x="P");
            fi;
         fi;     

        ###
        # TYPE C and D:
        # find a base such that
        #   C: have at most one "P"
        #   D: have at most one "P" in first < (n+1)/2 entries, or 1..1-11
        ###
         if tt in ["C","D"] and not [tt,rr]=["D",4] then
            pos := PositionsProperty(wh,x->x="P");
            while Length(pos)>1 and not (tt="D" and pos = [rank-1,rank]) do 
               newcg := applyReflection(newcg,pos[2],base);
               wh    := newcg.wh; 
               base  := newcg.base;
               newcg := newcg.cg;
               pos := PositionsProperty(wh,x->x="P");
              #Print("this is new wh ",wh,"\n");
            od;
            todo := [];
            if tt="D" and pos = [rank-1,rank] then 
               todo := Reversed([1..rank-1]); 
               for i in todo do
                  newcg := applyReflection(newcg,i,base);
                  wh    := newcg.wh; 
                  base  := newcg.base;
                  newcg := newcg.cg;
                  pos := PositionsProperty(wh,x->x="P");
                 #Print("this is new wh ",wh,"\n");
               od;
            fi;
           #bring P in first half
            if tt="D" and Length(pos)=1 and pos[1]>=(rank+1)/2 and not pos[1] in [rank,rank-1] then
               tmp  := rank-2-pos[1]+2;
               todo := List(Reversed([1..pos[1]]),x->List([1..tmp],y->x+y-1));
               todo := Concatenation(todo); 
               for i in todo do
                  newcg := applyReflection(newcg,i,base);
                  wh    := newcg.wh; 
                  base  := newcg.base;
                  newcg := newcg.cg;
                  pos := PositionsProperty(wh,x->x="P");
                 #Print("this is new wh ",wh,"\n");
               od;
            fi;
           #apply diag aut for D
            if tt="D" and pos = [rank] then
               tmp   := base[rank]; base[rank] := base[rank-1]; base[rank-1] := tmp;
               newcg := corelg.makeCanGenByBase(pr,cb,base);           
               wh    := where(cd,newcg);
               pos := PositionsProperty(wh,x->x="P");
            fi;
            if tt = "C" and (pos[1]>(rank)/2 and not pos[1]=rank) then
               tmp  := rank-1-pos[1];
               todo := List(Reversed([1..pos[1]]),x->List([0..tmp],y->x+y));
               todo := Concatenation(todo); 
               for i in todo do
                  newcg := applyReflection(newcg,i,base);
                  wh    := newcg.wh; 
                  base  := newcg.base;
                  newcg := newcg.cg;
                  pos := PositionsProperty(wh,x->x="P");
                 #Print("this is new wh ",wh,"\n");
               od;
            fi;
         fi; 


         ####
         # TYPE G2, F4, E6, E7, E8:
         # find base such that in std numbering of simple roots:
         #   G2: -11
         #   F4: 1-111 or 11-11 
         #   E6: (-111111), (1-11111)
         #   E7: (-1111111), (1-111111), (111111-1)
         #   E8: (-11111111), (1111111-1)
         ####
         if tt in ["G","F","E"] then
            pos  := PositionsProperty(wh,x->x="P");
            todo := [];
           #case G2
            if rr = 2 then
               if Length(pos)=2 then todo := [2]; fi;
               if pos = [1] then todo := [1,2]; fi;
           #case F4
            elif rr = 4 then
              #this is for the case that base has can ord [1,2,3,4]:
              #todoL := [[[],[]],[[2],[]],[[1,2,3],[2]],[[1,2,3,4],[3,2]],
              #       [[1,2,4],[4,3,2]],[[1,3],[1,2]],
              #       [[1,3,4],[1,3,2]],[[1,4],[1,4,3,2]],
              #       [[2,4],[2,3,2]],[[2,3,4],[2,4,3,2]],
              #       [[2,3],[3,2,4,3,2]],[[1,2],[2,3,2,4,3,2]],
              #       [[1],[1,2,3,2,4,3,2]],[[3],[]],[[3,4],[3]],[[4],[4,3]]];
              # todo := todoL[Position(List(todoL,x->x[1]),pos)][2];

              #this is for the case that base has can ord [2,4,3,1]:
              todoL:=[[[],[]],[[4],[]],[[2,3,4],[4]],[[1,2,3,4],[3,4]],[[1,2,4],[1,3,4]],
                     [[2,3],[2,4]],[[1,2,3],[2,3,4]],[[1,2],[2,1,3,4]],[[1,4],[4,3,4]],
                     [[1,3,4],[4,1,3,4]],[[3,4],[3,4,1,3,4]],[[2,4],[4,3,4,1,3,4]],
                      [[2],[2,4,3,4,1,3,4]],[[3],[]],[[1,3],[3]],[[1],[1,3]]];
              #todoL := [ [[],[]], [[1],[1,3]], [[2],[2,4,3,4,1,3,4]]];
               todo := todoL[Position(List(todoL,x->x[1]),pos)][2];

           #case E
            else
                while Length(pos)>1 and not (Length(pos)=2 and 2 in pos) do
                  tmp := pos[Length(pos)-1];
                  if tmp=2 then tmp:=pos[Length(pos)-2]; fi;
                  newcg := applyReflection(newcg,tmp,base);
                  wh    := newcg.wh; 
                  base  := newcg.base;
                  newcg := newcg.cg;
                  pos := PositionsProperty(wh,x->x="P");
                 #Print("this is new wh ",wh,"\n");
               od;               
              #E: now either 1x(-1) or 2x(-1) where lambda_2=-1
               if rr = 6 then       
                  todoL := [[[],[]],[[1],[]],[[1,3],[1]],[[3,4],[3,1]],
                           [[2,6],[6,5,4,3,1]],[[2,5],[2,4,3,1]],
		    	   [[3,5],[5,4,2,6,5,4,3,1]],
			   [[1,6],[1,3,4,2,5,4,3,1]],
			   [[1,2],[2,4,3,5,4,2,6,5,4,3,1]],
			   [[2,3],[3,1,4,3,5,4,2,6,5,4,3,1]],
			   [[4,5],[4,3,2,1,4,3,5,4,2,6,5,4,3,1]],
			   [[5,6],[5,4,3,2,1,4,3,5,4,2,6,5,4,3,1]],
			   [[6],[6,5,4,3,2,1,4,3,5,4,2,6,5,4,3,1]],
			   [[2],[]],[[2,4],[2]],[[3,6],[6,5,4,2]],
			   [[1,5],[1,3,4,2]],[[1,4],[4,2,1,5,4,3,6,5,4,2]]
			   ,[[4],[4,3,1,5,4,3,6,5,4,2]],
			   [[4,6],[4,3,2,1,4,3,6,5,4,2]],
			   [[5],[5,4,3,2,1,4,3,5,4,2]],
			   [[3],[3,4,2,5,4,3,6,5,4,2]]];
               elif rr = 7 then
                  todoL :=[[[],[]],[[7],[]],[[6,7],[7]],[[5,6],[6,7]],
                           [[4,5],[5,6,7]],[[2,3],[2,4,5,6,7]],
			   [[1,7],[7,6,5,4,3,2,4,5,6,7]],
			   [[1,2],[1,3,4,5,6,7]],
			   [[3,5],[3,1,4,3,2,4,5,6,7]],
			   [[2,6],[6,5,4,3,1,7,6,5,4,3,2,4,5,6,7]],
			   [[2],[]],[[2,4],[2]],[[3,7],[7,6,5,4,2]],
			   [[1,5],[1,3,4,2]],[[4,7],[4,3,1,5,4,3,6,5,4,2]]
			   ,[[5],[5,4,3,1,6,5,4,3,7,6,5,4,2]],[[1],[]],
			   [[1,3],[1]],[[3,4],[3,1]],
			   [[2,7],[7,6,5,4,3,1]],[[2,5],[2,4,3,1]],
			   [[3,6],[6,5,4,2,7,6,5,4,3,1]],
			   [[1,6],[1,3,4,2,5,4,3,1]],
			   [[1,4],[4,2,1,5,4,3,6,5,4,2,7,6,5,4,3,1]],
			   [[4],[4,3,1,5,4,3,6,5,4,2,7,6,5,4,3,1]],
			   [[4,6],[4,3,2,1,4,3,6,5,4,2,7,6,5,4,3,1]],
			   [[5,7],[5,4,3,2,1,4,3,5,4,2,7,6,5,4,3,1]],
			   [[6],[6,5,4,3,2,1,4,3,5,4,2,6,5,4,3,1]],
			   [[3],[3,4,2,5,4,3,6,5,4,2,7,6,5,4,3,1]] ];
               elif rr=8 then
                  todoL := [[[],[]],[[8],[]],[[7,8],[8]],[[6,7],[7,8]],[[5,6],[6,7,8]],
                            [[4,5],[5,6,7,8]],
                            [[2,3],[2,4,5,6,7,8]],[[1,8],[8,7,6,5,4,3,2,4,5,6,7,8]],
			    [[1,2],[1,3,4,5,6,7,8]],[[3,5],[3,1,4,3,2,4,5,6,7,8]],
			    [[2,7],[7,6,5,4,3,1,8,7,6,5,4,3,2,4,5,6,7,8]],
			    [[2,6],[2,4,3,1,5,4,3,2,4,5,6,7,8]],
			    [[3,6],[6,5,4,2,7,6,5,4,3,1,8,7,6,5,4,3,2,4,5,6,7,8]],
			    [[1,7],[1,3,4,2,5,4,3,1,6,5,4,3,2,4,5,6,7,8]],
			    [[1,4],[4,2,1,5,4,3,6,5,4,2,7,6,5,4,3,1,8,7,6,5,4,3,2,4,5,6,7,8]],
			    [[4],[4,3,1,5,4,3,6,5,4,2,7,6,5,4,3,1,8,7,6,5,4,3,2,4,5,6,7,8]],
			    [[4,6],[4,3,2,1,4,3,6,5,4,2,7,6,5,4,3,1,8,7,6,5,4,3,2,4,5,6,7,8]],
			    [[5,7],[5,4,3,2,1,4,3,5,4,2,7,6,5,4,3,1,8,7,6,5,4,3,2,4,5,6,7,8]],
			   [[6,8],[6,5,4,3,2,1,4,3,5,4,2,6,5,4,3,1,8,7,6,5,4,3,2,4,5,6,7,8]],
			    [[7],[7,6,5,4,3,2,1,4,3,5,4,2,6,5,4,3,1,7,6,5,4,3,2,4,5,6,7,8]],
			    [[3],[3,4,2,5,4,3,6,5,4,2,7,6,5,4,3,1,8,7,6,5,4,3,2,4,5,6,7,8]],
			    [[1],[]],[[1,3],[1]],[[3,4],[3,1]],[[2,8],[8,7,6,5,4,3,1]],
			    [[2,5],[2,4,3,1]],[[3,7],[7,6,5,4,2,8,7,6,5,4,3,1]],
			    [[1,6],[1,3,4,2,5,4,3,1]],[[4,8],[4,3,1,5,4,3,6,5,4,2,7,6,5,4,3,1]],
			    [[4,7],[4,3,2,1,4,3,7,6,5,4,2,8,7,6,5,4,3,1]],
			    [[1,5],[5,4,2,1,6,5,4,3,7,6,5,4,2,8,7,6,5,4,3,1]]
			    ,[[5],[5,4,3,1,6,5,4,3,7,6,5,4,2,8,7,6,5,4,3,1]
			    ],[[5,8],[5,4,3,2,1,4,3,5,4,2,8,7,6,5,4,3,1]],
			    [[6],[6,5,4,3,2,1,4,3,5,4,2,6,5,4,3,1]],
			    [[3,8],[3,4,2,5,4,3,6,5,4,2,7,6,5,4,3,1]],
			    [[2,4],[4,3,5,4,2,6,5,4,3,7,6,5,4,2,8,7,6,5,4,3,1]],
			    [[2],[2,4,3,5,4,2,6,5,4,3,7,6,5,4,2,8,7,6,5,4,3,1]]];
               fi;
               todo := todoL[Position(List(todoL,x->x[1]),pos)][2];
            fi;

           #Print("this is new wh ",wh,"\n");
            for i in todo do 
              #Print("act with ",i,"\n");
               newcg := applyReflection(newcg,i,base);
               wh    := newcg.wh; 
               base  := newcg.base;
               newcg := newcg.cg;
               pos := PositionsProperty(wh,x->x="P");
              #Print("this is new wh ",wh,"\n");
            od; 
         fi;
                     

        #Print(CartanType(corelg.CartanMatrixOfCanonicalGeneratingSet(L,newcg)),"\n");
        #Print("this is new wh",wh,"\n");

         tmp := testCFs(newcg);
         cft := tmp.cft;
         cfs := tmp.cfs;

        
         tmp := corelg.VoganDiagramOfRealForm(L,
                    rec( cg   := newcg, 
                         base := base, 
                         mv   := mv,
                         signs:= wh,  
                         cfsigma:=cfs,
                         cftheta:=cft));


         SetVoganDiagram(L,tmp);
         
         if Length(tmp!.param)=1 then
            mv := IdRealForm(L);
            SetRealFormParameters(L,RealFormParameters(RealFormById(mv)));
         fi;

         Info(InfoCorelg,2,"   end Vogan Diagram for simple LA");   

         tmp := VoganDiagram(L);
         tmp!.sstypes := [IdRealForm(L)];
         return tmp;
        

     ###########################################
     #here consider form of OUTER TYPE
     ###########################################
      else

        #find h0 to define new root ordering; take CSA compatible with h

         tmp := Intersection(cd.K,h);
         hs  := ShallowCopy( CanonicalGenerators(RootsystemOfCartanSubalgebra(cd.K,tmp))[3]);
         
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

        #adjust D_4 so that root 3 and 4 are swapped
         if ct.types[1] = ["D",4] then
            wh   := where(cd,newcg);
           #Print("this is 1st wh ",wh,"\n");
            pos  := Filtered([1..4],x->wh[x]="?");
            tmp  := Filtered([1..4],x->not x in pos and not x=2)[1];
            tmp  := [tmp,2,pos[1],pos[2]];
            if not IsDuplicateFreeList(tmp) then Error("ups..."); fi;
            base := base{tmp};
            newcg := corelg.makeCanGenByBase(pr,cb,base);
            es    := newcg[1];
            fs    := newcg[2];
         fi;


         sps := List( es, x -> Subspace( L, [x],"basis" ) );
         mv  := [];
         for i in [1..Length(es)] do
            j := PositionProperty( sps, U -> theta( es[i] ) in U );
            if j > i then
               Add(mv,[i,j]);
               es[j]:= theta( es[i] );
               fs[j]:= theta( fs[i] );
            fi;
        od;
        Sort(mv);
        notmv := Filtered([1..rank],x->not x in Flat(mv));
        newcg := [es,fs,List( [1..Length(es)], i -> es[i]*fs[i] ) ];

    
        #for roots not moved by theta, determine whether root space 
        #lies in K or in P
         wh := where(cd,newcg);  
    
        #Print("out, this is first wh ",wh,"\n");

        #now consider the Weyl group action to adjust it; the first case is E_6
         if rank=6 and Size(Filtered(wh,x->not x="?"))=2 and not (wh[2]=1 and wh[4]=1) then
            if wh[2]="P" and wh[4]="K" then
               newcg := applyReflection(newcg,2,base);
               wh    := newcg.wh; 
               base  := newcg.base;
               newcg := newcg.cg;
             
            fi;
            if wh[2]="P" and wh[4]="P" then     
               newcg := applyReflection(newcg,4,base);
               wh    := newcg.wh; 
               base  := newcg.base;
               newcg := newcg.cg;
            fi;
            if not wh[2]="K" and wh[4]="P" then Error("E6 error"); fi;
          
    
       #now this is case D
         elif Length(Filtered(wh,x->x="?"))=2 and "P" in wh and not Length(wh)=3 then
           #Print("this is wh", wh,"\n");
           
            pos := PositionsProperty(wh,x->x="P");
            while Length(pos)>1 do 
               newcg := applyReflection(newcg,pos[Length(pos)-1],base);
               wh    := newcg.wh; 
               base  := newcg.base;
               newcg := newcg.cg;
               pos := PositionsProperty(wh,x->x="P");
              #Print("this is pos ",pos,"\n");
            od;
            if pos[1] > First([1..rank],x-> x>= rank /2)-1 then 
               bbase := [];
               for i in [1..rank-2] do
                  bbase[i] := base[rank-i-1];
               od;
               bbase[rank-1] := -Sum(base{[1..rank-1]});
               bbase[rank]   := -Sum(base{[1..rank-2]})-base[rank];
               bcg   := corelg.makeCanGenByBase(pr,cb,bbase);
               iso   := LieAlgebraIsomorphismByCanonicalGenerators(L,newcg,L,bcg);
               newcg := List(newcg,x->List(x,y->Image(iso,y)));
               base  := bbase;
             
               wh := [];
               for i in newcg[1] do 
                  if i in cd.K then Add(wh,"K");
                  elif i in cd.P then Add(wh,"P"); 
                  else Add(wh,"?"); 
                  fi;
               od;
               for i in mv do
                  newcg[1][i[2]] := theta(newcg[1][i[1]]);
                  newcg[2][i[2]] := theta(newcg[2][i[1]]);
                  newcg[3][i[2]] := newcg[1][i[2]]*newcg[2][i[2]];
               od;
               pos := PositionsProperty(wh,x->x="P");
              #Print("... swapped this to ",pos,"\n");

            fi;
         fi;
 
         tmp := testCFs(newcg);
         cft := tmp.cft;
         cfs := tmp.cfs;
        #Print("this is new wh ",wh,"\n");

         tmp := corelg.VoganDiagramOfRealForm(L,
                    rec(cg:=newcg,  
                        base:=base, 
                        mv := mv,
                        signs:=wh, 
                        cfsigma:=cfs, 
                        cftheta:=cft));
         SetVoganDiagram(L,tmp);
        
         if Length(tmp!.param)=1 then
            mv := IdRealForm(L);
            SetRealFormParameters(L,RealFormParameters(RealFormById(mv)));
         fi;
 

         Info(InfoCorelg,2,"   end Vogan Diagram"); 
         tmp := VoganDiagram(L);
         tmp!.sstypes := [IdRealForm(L)];
         return tmp;
        

        #Add(res,rec(L := L, cd:=cd, h:=h, cg:=newcg, base := base, inn:=inn, mv := mv, cftheta := cft,
        #            cfsigma:=cfs, where := wh, rank:=Dimension(h)));

      fi;

  
end;




##############################################################################
# M is a simple LA of type "type" with can gens "cg", 
# "realification" is false if its a real form, and true if its a realification
# returns a root system of M, with corresp. can gen "cg"
#
#############################################################################
corelg.getRootsystem := function(M,cg,type,realification)
   local F,K,iso,cb,RK,R ,pr, ss;

   Info(InfoCorelg,3,"    start getRootsystem by can gen");
   F   := LeftActingDomain(M);
   K   := RealFormById( type[1],type[2],2,F);
   if realification then
      K := DirectSumOfAlgebras(K,K);
   fi;
   iso := LieAlgebraIsomorphismByCanonicalGenerators(
                K,CanonicalGenerators(RootSystem(K)),M,cg);
   cb  := List( ChevalleyBasis(K), x -> List( x, y -> Image( iso, y ) ) );

   if not cb[3]=cg[3] then Error("wrong last part"); fi;
   if not List(CanonicalGenerators(RootSystem(K)),x->List(x,i->Image(iso,i)))=cg then
      Error("wrong cg");
   fi;

   RK  := RootSystem(K);
      
   R:= Objectify( NewType( NewFamily( "RootSystemFam", IsObject ),
               IsAttributeStoringRep and IsRootSystemFromLieAlgebra ), 
               rec() );
   SetCanonicalGenerators( R, [cg[1],cg[2], cg[3]]);
   SetUnderlyingLieAlgebra( R, M );
   SetPositiveRootVectors( R, List(PositiveRootVectors(RK),x->Image(iso,x)));
   SetNegativeRootVectors( R, List(NegativeRootVectors(RK),x->Image(iso,x)));
   SetCartanMatrix( R, CartanMatrix(RK) );

   pr := PositiveRoots(RK);
   if F=SqrtField then pr := pr*One(SqrtField); pr:=SqrtFieldMakeRational(pr); fi;
   SetPositiveRoots(R, pr);

   pr := NegativeRoots(RK);
   if F=SqrtField then pr := pr*One(SqrtField); pr:=SqrtFieldMakeRational(pr); fi;
   SetNegativeRoots(R, pr);

   pr := SimpleSystem(RK);
   if F=SqrtField then pr := pr*One(SqrtField); pr:=SqrtFieldMakeRational(pr); fi;
   SetSimpleSystem(R, pr);
   
   SetChevalleyBasis(R, cb);
   SetRootSystem(M,R);
   Info(InfoCorelg,3,"    end getRootsystem by can gen");
   return R;
end;


################################################
corelg.makeBlockDiagMat := function ( mats )
   local  n, M, m, d;
   n := Sum( mats, x->Length(x[1]));
   M := NullMat( n, n );
   n := 0;
   for m  in mats  do
       d := Length( m );
       M{[ 1 .. d ] + n}{[ 1 .. d ] + n} := m;
       n := n + d;
   od;
   return M;
end;

###########################################################################
InstallMethod( VoganDiagram,
   "for Lie algebras",
   true,
   [ IsLieAlgebra ], 0, function(L)
local cd, h, rank, theta, sigma, R, C, base, cg, pr, ct, cb, algs, cgs,
      rA, a, mysort, en, ranks, i, testIt, makeCartInv, ha, cda,
      aK, aP, mv, rewr, hh, halg, perm, j, ctrf, ctreal;

   if HasVoganDiagram(L) then return VoganDiagram(L); fi;

   Info(InfoCorelg,1,"start Vogan Diagram; get CartDecomp and CSA");
   cd    := CartanDecomposition(L);
   h     := MaximallyCompactCartanSubalgebra(L);
   rank  := Dimension(h);
   theta := cd.CartanInv;
   sigma := RealStructure(L);
   if not ForAll(Basis(h),x->theta(x) in h) then
      Error("need a theta-stable CSA; Cartan Dec and CSA must be compatible!");
   fi;

   R    := RootsystemOfCartanSubalgebra(L,h);
   C    := CartanMatrix(R);
   base := SimpleSystem(R);
   cg   := CanonicalGenerators(R);
   pr   := PositiveRoots(R);  
   ct   := CartanType(C);
   cb   := ChevalleyBasis(R);
   en   := Concatenation(ct.enumeration);

   if Length(ct.types)=1 then
      Info(InfoCorelg,1,"call Vogan diagram for simple LA");
      return corelg.SingleVoganDiagram(L);
   fi;

  #identify realifications
   halg := List(ct.enumeration,x -> SubalgebraNC(L,Concatenation(List(cg,i->i{x}))));
   hh   := List(halg,x->Basis(x)[1]);
   perm := [];
   for i in [1..Length(hh)] do
      Add(perm,[i,First([1..Length(hh)],j-> theta(hh[i]) in halg[j])]);
   od;
   perm := AsSet(List(perm,AsSet));

  #have a single realification?
   if Length(perm)=1 and Length(perm[1])= 2 then
     #Print("start VD for real\n");
      return corelg.VoganDiagramRealification(L);
   fi;


  #now adjust cartantypes so that enums corresponding to
  #realifications are merged
   ctrf   := rec(types:=[],enumeration:=[]);
   for i in perm do
      if Length(i)=1 then 
         Add(ctrf.types,Concatenation(ct.types[i[1]],[1]));
         Add(ctrf.enumeration, ct.enumeration[i[1]]);
      else
         Add(ctrf.types,Concatenation(ct.types[i[1]],[2]));
         Add(ctrf.enumeration, 
             Concatenation(ct.enumeration[i[1]],ct.enumeration[i[2]]));
      fi;
   od;
   ct := ctrf;

  #adjust type F4
   for i in [1..Length(ct.types)] do
      if ct.types[i][1] = "F" then
         if ct.types[i][3] = 1 then
            ct.enumeration[i] := ct.enumeration[i]{[4,1,3,2]};
         elif ct.types[i][3] = 2 then
            ct.enumeration[i] := ct.enumeration[i]{[4,1,3,2,8,5,7,6]};
         fi;
      fi;
   od;

  #first sort by types
   mysort := function(a,b)
      local typ, pa, pb;
      typ := ["A","B","C","D","E","F","G"];
      pa  := Position(typ,a[1]);
      pb  := Position(typ,b[1]);
      if pa<pb then return true; fi;
      if pa>pb then return false; fi;
      if a[2]<b[2] then return true; fi;
      if a[2]>b[2] then return false; fi;
      return a[3]<b[3];
   end;
   SortParallel(ct.types,ct.enumeration,mysort);
   for i in [1..Length(ct.types)] do ct.types[i] := ct.types[i]{[1..2]}; od;

 
  #now sort by vectors
   en    := Concatenation(ct.enumeration);
   base  := List(base{en});
   rank  := List(ct.enumeration,Length);
   ranks := List([1..Length(rank)],x->[Sum(rank{[1..x-1]})+1..Sum(rank{[1..x]})]);
   cg    := corelg.makeCanGenByBase(pr,cb,base);

  #ct    := CartanType(corelg.CartanMatrixOfCanonicalGeneratingSet(L,cg));

   cgs  := List(ranks,x-> List(cg,i->i{x})); 
   algs := List(cgs,x->SubalgebraNC(L,corelg.myflat(x))); 
   

   #### small test if all is compatible
   testIt := function(L)
   local R,H,cb,cg,C,basH,h,i,r,cf,pr;

      H := CartanSubalgebra(L);
      R := RootSystem(L);
      pr := PositiveRoots(R);
      cb:= ChevalleyBasis(R);
      cg := CanonicalGenerators(R);
      C   := LieCentraliser(L,H);
      if not IsAbelian(C) or not Dimension(C)=Dimension(H) then
         Error("not a CSA");
      fi;
      if not cg[3] = cb[3]{[1..Length(cg[3])]} then Error("error 1"); fi;
      basH := Basis(H,cg[3]);
      for h in cg[3] do
         for i in [1..Length(pr)] do      
            r  := pr[i]*Coefficients(basH,h);
            cf := Coefficients(Basis(SubspaceNC(L,[cb[1][i]],"basis"),[cb[1][i]]),h*cb[1][i])[1];
            if not r=cf then Error("error cf"); fi;
         od;
      od;
     #Print("test ok\n");
   end;


   makeCartInv := function(a,aK,aP)
   local bas;
     bas := BasisNC(a,Concatenation(Basis(aK),Basis(aP)));
     return function(v)
     local k, p, cf, i;   
        k   := Length(Basis(aK));
        p   := Length(Basis(aP));
        cf  := List(Coefficients(bas,v),x->x);
        for i in [k+1..k+p] do cf[i] := -cf[i]; od;
           return cf*bas;
        end;
   end; 

  #here we reduce to each direct factor; compute its RS and Vogan diagram
   for i in [1..Length(algs)] do
      Info(InfoCorelg,2,"  start direct factor ",i," of type ",ct.types[i]);
      a   := algs[i];
      ha  := Intersection(a,h);
      SetCartanSubalgebra(a,ha);
      SetMaximallyCompactCartanSubalgebra(a,ha);
      SetRealStructure(a,sigma); ##!! or define wrt basis of L?
      aK  := Intersection(cd.K,a);
      SetCartanSubalgebra(aK,Intersection(aK,CartanSubalgebra(cd.K)));
      aP  := Intersection(cd.P,a);
      if not Dimension(aK) + Dimension(aP) = Dimension(a) then Error("dim"); fi;
      cda := rec(K:=aK, P:=aP, CartanInv:=makeCartInv(a,aK,aP));
      SetCartanDecomposition(a,cda);
      if not ForAll(Basis(aK),x->cda.CartanInv(x)=x) then Error("k"); fi;
      if not ForAll(Basis(aP),x->cda.CartanInv(x)=-x) then Error("p"); fi;
      if ct.types[i][2] = Length(ct.enumeration[i]) then
         rA := corelg.getRootsystem(a,cgs[i],ct.types[i],false);
      else
         rA := corelg.getRootsystem(a,cgs[i],ct.types[i],true);
      fi;
      SetRootSystem(ha,rA);
      Info(InfoCorelg,2,"  end direct factor ",i);
      VoganDiagram(a);
   od;

   mv := List(algs,x->MovedPoints(VoganDiagram(x)));
   for i in [2..Length(rank)] do
      mv[i] := mv[i]+Sum(List([1..i-1],j->rank[j]));
   od;
 
  #for rewriting signs-vector
   rewr := function(k) if k=1 then return "K"; elif k=-1 then return "P"; fi; end;

   a  := rec(cg   := List([1..3],i->
                     Concatenation(List(algs,x->CanonicalGenerators(VoganDiagram(x))[i]))),
             base := corelg.makeBlockDiagMat(List(algs,x->BasisOfSimpleRoots(VoganDiagram(x)))), 
             mv   := Concatenation(mv),
             signs:= List(Concatenation(List(algs,x->Signs(VoganDiagram(x)))),rewr),
             cfsigma:=Concatenation(List(algs,x->CoefficientsOfSigmaAndTheta(VoganDiagram(x)).cfsigma)), 
             cftheta:=Concatenation(List(algs,x->CoefficientsOfSigmaAndTheta(VoganDiagram(x)).cftheta)));
  #Print("this is a",a,"\n");
   a := corelg.VoganDiagramOfRealForm(L,a);

   if ForAll(algs,x->IsBound(VoganDiagram(x)!.sstypes)) then
      a!.sstypes := Concatenation(List(algs,a->VoganDiagram(a)!.sstypes));
   fi;
   SetVoganDiagram(L,a);
   if IsBound(a!.sstypes) then L!.sstypes := StructuralCopy(a!.sstypes); fi;


   Info(InfoCorelg,1,"end Vogan diagram");
   return a;
end);




############################################################################
#input is lie alg L and list l=[type,rank,pos of -1, outer?]
corelg.computeIdRealForm := function(L,l)
local id, nr, nr1, nr2;

   if l[1]="A" then
      if l[2]=1 and l[3]=1 then return [l[1],l[2],2]; fi;
      nr := First([1..l[2]],x-> x>= l[2]/2);
      if l[4] and l[3]=0 then return [l[1],l[2],1+nr+1];
      elif l[4] and l[3]>0 then return  [l[1],l[2],1+nr+2];
      elif not l[4] then return  [l[1],l[2],1+l[3]];
      fi;
   elif l[1] = "B" then
      return [l[1],l[2],l[3]+1];
   elif l[1] = "C" then
      if l[3] < NumberRealForms(l[1],l[2]) then return [l[1],l[2],l[3]+1]; fi;
      return [l[1],l[2], NumberRealForms(l[1],l[2]) ];
   elif l[1] = "D" then
      if l[2]>4 then
         nr1 := First([1..l[2]+1],x-> x> l[2]/2)-1;
         nr2 := First([1..l[2]+1],x-> x> (l[2]-1)/2)-1;
         if not l[4] then
            if l[3]=0 then return [l[1],l[2],1]; fi;
            if l[3]=l[2]-1 then return [l[1],l[2],2+nr1]; fi;
            return [l[1],l[2],1+l[3]];
         else
            if l[3]=0 then return [l[1],l[2],3+nr1]; fi;
            return [l[1],l[2],3+nr1+l[3]];
         fi;
      else
         if not l[4] then
            if l[3]=3 then return [l[1],l[2],2]; fi;
            if l[3]=2 then return [l[1],l[2],3]; fi;
         else
            if l[3] = 0 then return [l[1],l[2],5]; fi;
            if l[3] = 1 then return [l[1],l[2],4]; fi;
         fi;
      fi;
   elif l[1] = "G" then
      if l[3] = 0 then return [l[1],l[2],1]; fi;
      if l[3] = 2 then return [l[1],l[2],2]; fi;
   elif l[1] = "F" then
      if l[3] = 0 then return [l[1],l[2],1]; fi;
      if l[3] = 4 then return [l[1],l[2],2]; fi;
      if l[3] = 3 then return [l[1],l[2],3]; fi;
   elif l[1] = "E" then
      if l[2] = 8 then
         if l[3]=0 then return [l[1],l[2],1]; fi;
         if l[3]=8 then return [l[1],l[2],2]; else return [l[1],l[2],3]; fi;
      elif l[2] = 7 then
         if l[3]=0 then return [l[1],l[2],1]; fi;
         if l[3]=2 then return [l[1],l[2],2]; fi;
         if l[3]=7 then return [l[1],l[2],3]; fi;
         if l[3]=1 then return [l[1],l[2],4]; fi;
      elif l[2] = 6 then     
         if not l[4] then
            if l[3]=0 then return [l[1],l[2],1]; fi;
            if l[3]=2 then return [l[1],l[2],3]; fi;
            if l[3]=1 then return [l[1],l[2],4]; fi;
         else 
            if l[3]=4 then return [l[1],l[2],2]; fi;
            if l[3]=0 then return [l[1],l[2],5]; fi;
         fi;
     
      fi;
   fi;  
end;




##############################################################################
#
#
InstallGlobalFunction( IdRealForm, function(L)
local id,vd,pos,tmp,j;

   if IsBound(L!.id) then return L!.id; fi;
   if IsBound(L!.sstypes) then return L!.sstypes; fi;
  
   if (HasIsCompactForm(L) and IsCompactForm(L)) or Dimension(CartanDecomposition(L).P)=0 then
      id := CartanType(CartanMatrix(VoganDiagram(L))).types;
      if Length(id) = 1 then
         id := id[1];
         Add(id,1);
      else
         id := ShallowCopy(id);
         for j in id do Add(j,1); od;
      fi;
      L!.id := id;
      return id;
   fi;

   if HasRealFormParameters(L) then
      tmp := RealFormParameters(L);
      pos := Position(tmp[3],-1);
      if pos = fail then pos := 0; fi;
      id  := corelg.computeIdRealForm(L,[tmp[1],tmp[2],pos,
             Length(Filtered(Orbits(Group(tmp[4]),[1..tmp[2]]),
                    x->Length(x)=2))>0]);  
   else
      vd  := VoganDiagram(L);
      if IsBound(L!.sstypes) then return L!.sstypes; fi;
      tmp := vd!.param;
      if Length(tmp)>1 then
         Error("id functionality only for simple LAs");
      else
         tmp := tmp[1];
      fi;
      pos := Position(Signs(vd),-1);
      if pos = fail then pos := 0; fi;
      id := corelg.computeIdRealForm(L,[tmp[1],tmp[2],pos,
            Length(MovedPoints(vd))>0]);  
   fi;
   L!.id := id;
   return id;

end);


###############################################################################
InstallGlobalFunction( NumberRealForms, function(t,r)
local nr, mv, nr1, nr2;
    if t = "A" then
      #if not r > 1 then Error("rank must be at least 2"); fi;
       if r = 1 then return 2; fi;
       nr := First([1..r],x-> x>= r/2);
       mv := 1; if IsOddInt(r) then mv := 2; fi;   
       return 1+nr+mv;
    elif t = "B" then
       if not r > 1 then Error("rank must be at least 2"); fi;
       return r+1;
    elif t = "C" then
       if not r > 2 then Error("rank must be at least 3"); fi;
       return First([1..r],x-> x> r/2)+1;
    elif t = "D" then
       if not r > 1 then Error("rank must be at least 4"); fi;
       nr1 := First([1..r+1],x-> x > r/2)-1;
       nr2 := First([1..r+1],x-> x > (r-1)/2)-1;
       if r >4 then return 1+nr1+1+1+nr2; else return 5; fi;
    elif t = "G" then
       if not r =2 then Error("rank must be 2"); fi;
       return 2;
    elif t = "F" then
       if not r =4 then Error("rank must be 4"); fi;
       return 3;
    elif t = "E" then
       if r=6 then return 5; elif r=7 then return 4; elif r=8 then return 3; fi;
       Error("rank must be 6,7, or 8"); 
    fi;
end);



####################################################################################
InstallGlobalFunction( RealFormsInformation, function(t,r)
local nr, mv, nr1, nr2, en;
    Print("\n");
    if t = "A" and r=1 then
        Print("  There are 2 simple real forms with complexification ",t,r,"\n");
        Print("    1 is of type su( 2 ), compact form\n");
        Print("    2 is of type su(1,1)=sl(2,R)\n");
    elif t = "A" and r>1 then
        nr := First([1..r],x-> x>= r/2);
        mv := 1; if IsOddInt(r) then mv := 2; fi;       
        Print("  There are ",1+nr+mv, " simple real forms with complexification ",t,r,"\n");
        Print("    1 is of type su(",r+1,"), compact form\n");
        Print("    2 - ",nr+1," are of type su(p,",r+1,"-p) with 1 <= p <= ",nr,"\n");
        if mv = 1 then
           Print("    ",nr+2," is of type sl(",r+1,",R)\n");
        elif mv = 2 then
           Print("    ",nr+2," is of type sl(",(r+1)/2,",H)\n");
           Print("    ",nr+3," is of type sl(",r+1,",R)\n"); 
        fi;
    elif t = "B" then
        if not r > 1 then Error("rank must be at least 2"); fi;
        Print("  There are ",r+1, " simple real forms with complexification ",t,r,"\n");
        Print("    1 is of type so(",2*r+1,"), compact form\n");
        Print("    2 - ",r+1," are of type so(2*p,",2*r,"-(2*p)+1) with 1 <= p <= ",r,"\n");
    elif t = "C" then
       if not r > 2 then Error("rank must be at least 3"); fi;
       nr := First([1..r],x-> x> r/2)+1;
       Print("  There are ",nr, " simple real forms with complexification ",t,r,"\n");
       Print("    1 is of type sp(",r,"), compact form\n");
       Print("    2 - ",nr-1," are of type sp(p,",r,"-p) with 1 <= p <= ",nr-2,"\n");
       Print("    ",nr," is of type sp(",r,",R)\n");
    elif t = "D" then
       if not r > 3 then Error("rank must be at least 4"); fi;
       nr1 := First([1..r+1],x-> x> r/2)-1;
       nr2 := First([1..r+1],x-> x> (r-1)/2)-1;
       if r = 4 then nr1 := nr1-1; fi;
       Print("  There are ",1+nr1+1+1+nr2, " simple real forms with complexification ",t,r,"\n");
       Print("    1 is of type so(",2*r,"), compact form\n");
       if r > 4 then
          Print("    2 - ",nr1+1," are of type so(2p,",2*r,"-2p) with 1 <= p <= ",nr1,"\n");
          Print("    ",nr1+2," is of type so*(",2*r,")\n");
          Print("    ",nr1+3," is of type so(",2*r-1,",1)\n");
          Print("    ",nr1+4," - ",nr1+3+nr2," are of type so(2p+1,",2*r,"-2p-1) with 1 <= p <= ",nr2,"\n");
       else
          Print("    2 is of type so*(8)\n");
          Print("    3 is of type so(4,4)\n");
          Print("    4 is of type so(3,5)\n");
          Print("    5 is of type so(1,7)\n");    
        # Print("    4 is of type so(",2*r-1,",1)\n");
#       #  Print("    5 is of type so(3,5)\n");
## corrected (3,5) and (1,7) (swap and typo)
       fi;
       
    elif t = "G" then
       if not r=2 then Error("rank must be 2"); fi;
       Print("  There are ",2, " simple real forms with complexification ",t,r,"\n");
       Print("    1 is the compact form\n");
       Print("    2 is G2(2) with k_0 of type su(2)+su(2) (A1+A1)\n");
    elif t = "F" then
       if not r =4  then Error("rank must be 4"); fi;
       Print("  There are ",3, " simple real forms with complexification ",t,r,"\n");
       Print("    1 is the compact form\n");
       Print("    2 is F4(4) with k_0 of form sp(3)+su(2) (C3+C1)\n"); #signs Params [11-11]
       Print("    3 is F4(-20) with k_0 of form so(9) (B4)\n");       #signs Params [1-111]
    elif t = "E" then
       Print("  There are ",NumberRealForms(t,r)," simple real forms with complexification ",t,r,"\n");
       Print("    1 is the compact form\n");
       if r = 6 then
          Print("    2 is EI   = E6(6), with k_0 of type sp(4) (C4)\n");
          Print("    3 is EII  = E6(2), with k_0 of type su(6)+su(2) (A5+A1)\n");
          Print("    4 is EIII = E6(-14), with k_0 of type so(10)+R (D5+R)\n");
          Print("    5 is EIV  = E6(-26), with k_0 of type f_4 (F4)\n");
       elif r=7 then
          Print("    2 is EV   = E7(7), with k_0 of type su(8) (A7)\n");
          Print("    3 is EVII  = E7(-25), with k_0 of type e_6+R (E6+R)\n"); #so(12)+su(2)\n");
### notation swap: EVII and EVI changed
          Print("    4 is EVI = E7(-5), with k_0 of type so(12)+su(2) (D6+A1)\n");#e_6+R (\n");
       elif r=8 then
          Print("    2 is EVIII = E8(8), with k_0 of type so(16) (D8)\n");
          Print("    3 is EIX   = E8(-24), with k_0 of type e_7+su(2) (E7+A1)\n");
       fi;
       if not r in [6,7,8] then Error("rank must be 6,7, or 8"); fi;
    fi;
    
    Print("  Index '0' returns the realification of ",t,r,"\n\n");
end);


#####################################################################################

InstallGlobalFunction( RealFormById, function(arg)
local r,t,id, par,sign,perm,mv,nr,tmp, nr1,nr2, F, cf,testCF, vd, rsc, en, sigma, cg, ct, base, L;

    if IsField(arg[Length(arg)]) then 
        F := arg[Length(arg)]; 
        arg := arg{[1..Length(arg)-1]};
    else 
        F := SqrtField; 
    fi;
  
    if IsList(arg[1]) and IsList(arg[1][1]) and IsList(arg[1][1][1]) then
        L := RealFormById(arg[1][1],F);
        for r in [2..Size(arg[1])] do
           L := DirectSumOfAlgebras(L,RealFormById(arg[1][r],F));
        od;
        L!.sstypes := arg[1];
        return L;
    fi;

    
    if IsList(arg[1]) and Length(arg[1])>1 then arg := arg[1]; fi;
    t:=arg[1]; r:=arg[2]; id:=arg[3];
    if not id in [0..NumberRealForms(t,r)] then 
        Error("there are only ",NumberRealForms(t,r), " real forms"); 
    fi;
    par  := [t,r];
    sign := ListWithIdenticalEntries(r,One(F));
    perm := ();

   #realification of simple
    if id = 0 then
       tmp := corelg.realification(t,r,F);
       sigma := RealStructure(tmp);
       tmp!.id  := [t,r,0];
       if F = SqrtField then tmp!.std := true; fi;
       SetIsRealFormOfInnerType(tmp,false);
       SetIsCompactForm(tmp,false);
       SetIsRealification(tmp,true);
       rsc := RootsystemOfCartanSubalgebra(tmp);

       ct   := CartanType(CartanMatrix(rsc));
       en   := Concatenation( ct.enumeration );
       if ct.types[1] = ["F",4] then
          en := en{[4,1,3,2,   8,5,7,6]};
       fi;             
       base := SimpleSystem(rsc){en};
       cg := corelg.makeCanGenByBase(PositiveRoots(rsc),ChevalleyBasis(rsc),base);

       vd  := corelg.VoganDiagramOfRealForm(tmp,
                rec(cg   := cg,
                    base := base,
                    mv   := List([1..r],x->[x,x+r]),
                    signs:= List([1..2*r],x->"?"),
                    cfsigma:= List([1..2*r],x->-One(F)),
                    cftheta:= List([1..2*r],x->One(F))));
  
       SetVoganDiagram(tmp,vd);
       SetcorelgCompactDimOfCSA(CartanSubalgebra(tmp),r);
       return tmp;
    fi;

   #compact form
    if id = 1 then 
       tmp     := corelg.Sub3( t, r, F );
       tmp!.id := [t,r,id];
       SetIsRealFormOfInnerType(tmp,true);
       SetRealFormParameters(tmp,[t,r,ListWithIdenticalEntries(r,1),()]);  
       SetIsRealification(tmp,false);
       rsc := RootsystemOfCartanSubalgebra(tmp);
       vd  := corelg.VoganDiagramOfRealForm(tmp,
                rec(cg   := CanonicalGenerators(rsc), 
                    base := SimpleSystem(rsc),
                    mv   := Filtered(Orbits(Group(perm),[1..r]),x->Length(x)=2),
                    signs:= List([1..r], function(i)
                              if sign[i]=-1 then return "P"; fi;
                              if not i^perm=i then return "?"; fi;
                              if i^perm=i and sign[i]=1 then return "K"; fi; end), 
                    cfsigma:= -sign, 
                    cftheta:= sign));
       vd!.sstypes := [ [t,r,id] ];
       tmp!.sstypes := [[t,r,id]];
       SetVoganDiagram(tmp,vd);
       if F=SqrtField then tmp!.std := true; fi;
       SetcorelgCompactDimOfCSA(CartanSubalgebra(tmp),Dimension(CartanSubalgebra(tmp)));
       return tmp;
    fi;

    if t = "A" and r=1 and id=2 then
       sign := [-One(F)];
       perm := ();
    elif t = "A" and r > 1 then
       nr := First([1..r],x-> x>= r/2);
       if id in [2..nr+1] then 
          sign[id-1] := -One(F);
          perm := ();
       fi; 
       if id >= nr+2 then perm := PermList(Reversed([1..r])); fi;
       if id = nr+3 then sign[(r+1)/2] := -One(F); fi;
       Add(par,sign); Add(par,perm);
    elif t = "B" then
       sign[id-1] := -One(F);
       perm := ();
    elif t = "C" then
       if id < NumberRealForms(t,r) then sign[id-1] := -One(F); else sign[r] := -One(F); fi;
       perm := ();
    elif t = "D" then
       if r>4 then
          nr1 := First([1..r],x-> x> r/2)-1;
          nr2 := First([1..r+1],x-> x>(r-1)/2)-1;
          perm := ();
          if id-1 =nr1+1 then sign[r-1] := -One(F); fi;
          if id-1< nr1+1 then sign[id-1] := -One(F); fi;
          if id-1 >nr1+1 then perm := (r-1,r); fi;
          if id-1 >nr1+2 then sign[id-nr1-2-1]:=-One(F); fi;
      else
         if id=2 then sign[3]:=-One(F); fi;
         if id=3 then sign[2]:=-One(F); fi;
         if id=4 then sign[1]:=-One(F); perm:=(3,4); fi;
         if id=5 then perm:=(3,4); fi;
      fi;
    elif t = "G" then
       perm    := ();
       sign[2] := -One(F);
    elif t = "F" then
       perm     := ();
       if id = 2 then sign := [1,1,1,-1]*One(F); fi;
       if id = 3 then sign := [1,1,-1,1]*One(F); fi;
    elif t = "E" then
       perm := ();
       if r=7 and id=2 then sign[2] := -One(F); fi;
       if r=7 and id=3 then sign[7] := -One(F); fi;
       if r=7 and id=4 then sign[1] := -One(F); fi;
       if r=8 and id=2 then sign[1] := -One(F); fi;
       if r=8 and id=3 then sign[8] := -One(F); fi;
       if r=6 then
          if id=2 then perm:=(1,6)(3,5); sign[4]:=-One(F); fi;
          if id=3 then sign[2] := -One(F); fi;
          if id=4 then sign[1] := -One(F); fi;
          if id=5 then perm:=(1,6)(3,5); fi;
       fi;
    fi;

   #tmp := corelg.NonCompactRealFormsOfSimpleLieAlgebra([t,r,sign,perm],F);
    tmp := corelg.SuperLie( t, r, sign, perm, F );
    SetRealFormParameters(tmp, [t,r,sign,perm]);
    if E(4) in F or IsSqrtField(F) then sigma := RealStructure(tmp); fi;  


    if perm=() then SetIsRealFormOfInnerType(tmp,true); else SetIsRealFormOfInnerType(tmp,false); fi;
    SetIsRealification(tmp,false);
   #attach Vogan diagram if not compact form
    rsc := RootsystemOfCartanSubalgebra(tmp);
    vd  := corelg.VoganDiagramOfRealForm(tmp,
                rec(cg   := CanonicalGenerators(rsc), 
                    base := SimpleSystem(rsc),
                    mv   := Filtered(Orbits(Group(perm),[1..r]),x->Length(x)=2),
                    signs:= List([1..r], function(i)
                              if sign[i]=-1 then return "P"; fi;
                              if not i^perm=i then return "?"; fi;
                              if i^perm=i and sign[i]=1 then return "K"; fi; end), 
                    cfsigma:= -sign, 
                    cftheta:= sign));
    vd!.sstypes := [ [t,r,id] ]; 
    SetVoganDiagram(tmp,vd);
 
    tmp!.id  := [t,r,id];
    tmp!.sstypes := [[t,r,id]];

    if F = SqrtField then  tmp!.std := true; fi;
    SetcorelgCompactDimOfCSA(CartanSubalgebra(tmp),
             Dimension(Intersection(CartanDecomposition(tmp).K,CartanSubalgebra(tmp))));

    return tmp;
 
end);


#############################################################################
corelg.getPureLA := function(arg)
local L, M,F,t,r,id;
   if IsField(arg[Length(arg)]) then 
       F := arg[Length(arg)]; 
       arg := arg{[1..Length(arg)-1]};
   else 
       F := SqrtField; 
   fi;
   if IsList(arg[1]) and Length(arg[1])>1 then arg := arg[1]; fi;
   t:=arg[1]; r:=arg[2]; id:=arg[3];
   L:=RealFormById(t,r,id,F);;
   L:=LieAlgebraByStructureConstants(F,StructureConstantsTable(Basis(L)));;
   return L;
end;

corelg.getDirectSumOfPureLA := function(arg)
local l,L,v,F;
   if Length(arg)=1 then F:=SqrtField; else F:=arg[2]; fi;
   v := arg[1];
   l := List(v,x->corelg.getPureLA(x,F));
   L := DirectSumOfAlgebras(l);
   return rec(alg:=L, algs:=l);
end;


###############################################################################
InstallGlobalFunction( AllRealForms, function(arg)
local F, t,r;
   t := arg[1];
   r := arg[2];
   if Length(arg) = 3 then F := arg[3]; else F := SqrtField; fi;
   return List([1..NumberRealForms(t,r)],x->RealFormById(t,r,x,F));
end);



##############################################################################
# 
#
InstallGlobalFunction( IsomorphismOfRealSemisimpleLieAlgebras, function(L1,L2)
local H1, H2, cd1, cd2, mkWhere, L, cd, whs, K1, K2, phi, l, cf,rank, res, vd,  
      cm1, cm2, cm3;

   if not IsLieAlgebra(L1) or not IsLieAlgebra(L2) then return false; fi;
   if not Dimension(L1)=Dimension(L2) then return false; fi;
   if not LeftActingDomain(L1)=LeftActingDomain(L2) then return false; fi;
   H1  := MaximallyCompactCartanSubalgebra(L1);
   H2  := MaximallyCompactCartanSubalgebra(L2);
   if not Dimension(H1)=Dimension(H2) then return false; fi;
   cd1 := CartanDecomposition(L1);
   cd2 := CartanDecomposition(L2);
   if not Dimension(cd1.P)=Dimension(cd2.P) then return false; fi;
   if not Dimension(cd1.K)=Dimension(cd2.K) then return false; fi;

 #now do basically the same as for IsomorphismsOfRealForms
 
   res := [];

   mkWhere := function(signs,mv)
   local i, new;
      new :=[];
      for i in [1..Length(signs)] do
         if signs[i]=-1 then Add(new,"P");
         elif i in Flat(mv) then Add(new,"?");
         else Add(new,"K");
         fi;
     od;
     return new;
   end;

   for L in [L1,L2] do 
      vd := VoganDiagram(L); 
      Add(res, rec(L:=L, cd := CartanDecomposition(L), 
                   cg := List(CanonicalGenerators(vd),x->List(x,y->y)),
                   base := BasisOfSimpleRoots(vd),  
                   mv := MovedPoints(vd),  
                   where := mkWhere(Signs(vd),MovedPoints(vd)),
                   type := vd!.param,    cftheta :=  CoefficientsOfSigmaAndTheta(vd).cftheta, 
                                         cfsigma := CoefficientsOfSigmaAndTheta(vd).cfsigma));
   od;


   whs  := List(Collected(List(res,x->[x.type,x.where])),x->x[1]);
   if Length(whs)=2 then return false; fi;

   K1   := res[1];
   rank := Sum(List(K1.type,x->x[2]));
   K2   := res[2];
   cf   := List([1..rank],x->Sqrt(K1.cfsigma[x]/K2.cfsigma[x]));
   if not ForAll(cf,x-> x in LeftActingDomain(L1)) then Error("cannot do this over ",LeftActingDomain(L1)); fi;
   if not ForAll(cf,x-> x = ComplexConjugate(x)) then Error("ups, not real"); fi;

   if not ForAll(cf,x-> x in Rationals or not SqrtFieldMakeRational(x)=false) 
         and not cf[1] in SqrtField then
      Print("(isom would have to be defined over SqrtField!)\n");
      return fail;
   else
      for l in [1..rank] do 
         K2.cg[1][l] := cf[l]*K2.cg[1][l];
         K2.cg[2][l] := cf[l]^-1*K2.cg[2][l];
      od;
      phi := LieAlgebraIsomorphismByCanonicalGenerators(K1.L,K1.cg,K2.L,K2.cg);
   fi;
      
  #Info(InfoCorelg,1,"  now test isom");
  #cd1 := CanonicalGenerators(VoganDiagram(L1));
  #cd2 := CanonicalGenerators(VoganDiagram(L2));;
  #cm1 := corelg.CartanMatrixOfCanonicalGeneratingSet(L1,cd1);
  #cm2 := corelg.CartanMatrixOfCanonicalGeneratingSet(L2,cd2);
  #if not ForAll(Flat(cd1),x->x in L1) or not ForAll(Flat(cd2),x->x in L2) then
  #   Error("wrong can gen sets in VoganDiag");
  #fi;
  #cm3 := corelg.CartanMatrixOfCanonicalGeneratingSet(L2,List(cd1,x->List(x,i->Image(phi,i))));
  #if not cm1=cm2 or not cm1=cm3 then
  #   Error("isom wrong! (at least CM different...)");
  #fi;
  #Info(InfoCorelg,1,"  test ok");

   return phi;
end);






#########################################################################################################
#########################################################################################################
#########################################################################################################


corelg.splitRealFormOfSL := function(rank)
local M, S, K, P, i, j, l, tmp, diag, bas, h, makeCartInv, n, eij, R, T, SS, writeToSS, KK, PP, RR, iso, F;

    n    := rank+1;
    F    := SqrtField;
    M    := MatrixLieAlgebra(F,n);
    bas  := BasisVectors(Basis(M));
    eij  := function(i,j) return bas[(i-1)*n+j]; end;
    tmp  := Filtered(Basis(M),x->Trace(x)=Zero(F)); 
    h    := List([1..rank],i-> eij(i,i)-eij(i+1,i+1));
    tmp  := Concatenation(tmp,h);
    S    := SubalgebraNC(M,tmp,"basis");
    SetCartanSubalgebra(S,Subalgebra(S,h));
    R    := RootsystemOfCartanSubalgebra(S);
    SetRootSystem(S,R);
    
   #have K = X in S with X^\intercal = -X
   #and P = X in S with X^\intercal = X
    K := [];
    P := ShallowCopy(h);
    for i in [1..n-1] do
       for j in [i+1..n] do
           Add(K,eij(i,j)-eij(j,i));
           Add(P,eij(i,j)+eij(j,i));
       od;
    od;
    K := SubalgebraNC(S,K,"basis");
    P := SubspaceNC(S,P,"basis");

   
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

   SetCartanDecomposition( S, rec( K:= K, P:= P,
                          CartanInv := makeCartInv(S,K,P)));

   T   := StructureConstantsTable(Basis(S));;
   SS  := LieAlgebraByStructureConstants(F,T);
   bas := Basis(SS);

   writeToSS := v->Coefficients(Basis(S),v)*Basis(SS);
   SetCartanSubalgebra(SS,SubalgebraNC(SS,List(Basis(CartanSubalgebra(S)),writeToSS)));
   RR := RootsystemOfCartanSubalgebra(SS);
   SetRootSystem(SS,RR);
   KK  := SubalgebraNC(SS,List(Basis(K),writeToSS));
   PP  := SubspaceNC(SS, List(Basis(P),writeToSS));
   iso := AlgebraHomomorphismByImagesNC(SS,S,Basis(SS),Basis(S));
   SetCartanDecomposition(SS, rec(K := KK,
                                  P := PP,
                                  CartanInv := makeCartInv(SS,KK,PP)));

   return rec(matalg := S, liealg := SS, iso := iso);
end;







#########################################
corelg.STKDTA := function(L)

    local H, R, cd, c, a, prv, nrv, pr, posR, cfs, i, j, cfc, cfa, t, sums, D, B, C, en, cc,
          pos, cpt, cf, sp, pairs, posres, rts, resbas, posRv, negRv,
          posresRv, negresRv, zerosp, C0, tp, hts, p, R0, eqns, eq, d, bh, mat, sol;


    H:= MaximallyNonCompactCartanSubalgebra(L);
    R:= RootsystemOfCartanSubalgebra( L, H );
    cd:= CartanDecomposition(L);

    bh:= Basis(H,CanonicalGenerators(R)[3]);
    mat:= List( bh, h -> Coefficients( bh, cd.CartanInv(h) ) );
    sol:= NullspaceMat( mat + One(LeftActingDomain(L))*IdentityMat(Length(mat)) );
    c:= List( sol, x -> x*bh );

    #c:= BasisVectors( Basis( Intersection( cd.P, H ) ) );
    a:= BasisVectors( Basis( Intersection( cd.K, H ) ) );

    prv:= PositiveRootVectors( R );
    nrv:= NegativeRootVectors( R );

    pr:= PositiveRootsNF(R);
    posR:= [ ];
    cfs:= [ ];
    posRv:= [ ];
    negRv:= [ ]; 
    for i in [1..Length(pr)] do
        cfc:= List( c, h -> Coefficients( Basis( SubspaceNC(L,[prv[i]],"basis"),[prv[i]]), h*prv[i])[1]);
        cfa:= List( a, h -> 
                  E(4)*Coefficients( Basis( SubspaceNC(L,[prv[i]],"basis"),[prv[i]]), h*prv[i])[1]);
        t:= First( cfc, x -> x<>0 );
        if t <> fail then
           if t > 0 then
              Add( posR, pr[i] );
              Add( posRv, prv[i] );
              Add( negRv, nrv[i] );
              Add( cfs, [cfc,cfa] );
           else
              Add( posR, -pr[i] );
              Add( posRv, nrv[i] );
              Add( negRv, prv[i] );
              Add( cfs, [-cfc,-cfa] );
           fi;
        else
           #cf:= List( a, h -> 
           #       E(4)*Coefficients( Basis( SubspaceNC(L,[prv[i]],"basis"),[prv[i]]), h*prv[i])[1]);
           #t:= First( cfa, x -> x<>0 );
           Add( posR, pr[i] );
           Add( posRv, prv[i] );
           Add( negRv, nrv[i] );
           Add( cfs, [cfc,cfa] );
        fi;
    od;

    sums:= [ ];
    for i in posR do for j in posR do Add( sums, i+j ); od; od;

    D:= Filtered( posR, x -> not x in sums );
    B:= BilinearFormMatNF(R);
    C:= List( D, x -> List( D, y -> 2*x*B*y/(y*B*y) ) );
    en:= Concatenation( CartanType(C).enumeration );
    D:= D{en};
    C:= C{en}{en};

    cpt:= [ ];
    for i in [1..Length(D)] do
        pos:= Position( posR, D[i] );
        cc:= cfs[pos];
        if IsZero( cc[1] ) then Add( cpt, i ); fi;
    od;

    sp:= Basis( VectorSpace( Rationals, D ), D );
    pairs:= [ ];
    for i in [1..Length(D)] do
        if not i in cpt then
           pos:= Position( posR, D[i] );
           cc:= List( cfs[pos], ShallowCopy );
           cc[2]:= -cc[2];
           pos:= Position( cfs, cc );
           cc:= Coefficients( sp, -posR[pos] );
           for j in [1..Length(D)] do
               if cc[j] = -1 and not j in cpt and i < j then
                  Add( pairs, [i,j] );
               fi;
           od;
        fi;
    od;

    return rec( CM:= C, cpt:= cpt, sym:= pairs, bas:= D, posR:= posR,
                cfs:= cfs, posRv:= posRv, negRv:= negRv );

end;


#############################################################################
InstallMethod( ViewObj,
   "for Satake diagram",
   true,
   [ IsSatakeDiagramOfRealForm ], 0,
   function( o )
   local tmp, i;
      tmp := ""; #Concatenation(o!.type,String(o!.rank));
      for i in [1..Length(o!.type)] do
          Append( tmp, o!.type[i][1] );
          Append( tmp, String(o!.type[i][2]) );
          if i < Length(o!.type) then Append( tmp, "x" ); fi;
      od;
      Print(Concatenation(["<Satake diagram in Lie algebra of type ",tmp,">"]));
end );


#############################################################################
#
# for printing Satake and Vogan diagrams
#
corelg.prntdg:= function( C, blc )

   local t, en, type, i, b, s, rank, offset, bound,lv;

   t:= CartanType(C);
   for lv in [1..Length(t.enumeration)] do
      if not lv = 1 then Print("\n"); fi;
      en:= t.enumeration[lv];
      rank:= Length(en);

      type:= t.types[lv][1];

      if type ="D" then
         b:= Length( Intersection( blc, en{[1..Length(en)-2]} ));
         offset:= 3+3*b+(Length(en)-2-b) + 3*(rank-3)+2;
         s:= ""; for i in [1..offset] do Append( s, " "); od;
         Append(s," ");
         if en[rank-1] in blc then
            Append( s, "("); Append( s, String(en[rank-1]) ); Append( s, ")");
         else
            Append( s, String(en[rank-1]) );
         fi;
         Append(s,"\n");
         Print(s);
         s:= ""; for i in [1..offset] do Append( s, " "); od;
         Append( s, "/\n" );
         Print(s);
      fi;

      if type ="E" then
         b:= 3;
         if en[1] in blc then
            b:= b+3;
         else
            b:= b+1;
         fi;
         if en[3] in blc then
            b:= b+3;
         else
            b:= b+1;
         fi;
         if en[4] in blc then
            b:= b+2;
         else
            b:= b+1;
         fi;
         b:= b+6;

         offset:= b;
         s:= ""; for i in [1..offset] do Append( s, " "); od;
         if not en[2] in blc then Append( s, " " ); fi;
         if en[2] in blc then
            Append( s, "("); Append( s, String(en[2]) ); Append( s, ")");
         else
            Append( s, String(en[2]) );
         fi;
         Append(s,"\n");
         Print(s);
         s:= " "; for i in [1..offset] do Append( s, " "); od;
         Append( s, "|\n" );
         Print(s);
      fi;  

      Print( t.types[lv][1], t.types[lv][2], ":  " );
      if type in ["A","B","C","F","G","E"] then
         bound:= rank;
      elif type = "D" then
         bound:= rank-2;
      fi; 
         for i in [1..bound] do
          if type <> "E" or i <> 2 then
             if en[i] in blc then
                Print("(",en[i],")");
             else
                Print(en[i]);
             fi;
          fi;

          if i < Length(en) then
      
             if type = "A" or (type in ["B","C"] and i < Length(en)-1) 
                     or (type = "D" and i < Length(en)-2) or (type="E" and i <> 2) 
                     or (type = "F" and i in [1,3] ) then
                Print( "---");
             elif type in ["B","C"] and i = Length(en)-1 then
                if type = "B" then
                   Print("=>=");
                else
                   Print("=<=");
                fi;
             elif type = "G" then
                Print("#>#");
             elif type ="F" and i=2 then
                Print("=>=");
             fi;
          fi;
      od;

      if type ="D" then
         s:= "\n"; for i in [1..offset] do Append( s, " "); od;
         Append( s, "\\\n" );
         Print(s);
         b:= Length( Intersection( blc, en{[1..Length(en)-2]} ));
         offset:= 3+3*b+(Length(en)-2-b) + 3*(rank-3)+2;
         s:= ""; for i in [1..offset] do Append( s, " "); od;
         Append(s," ");
         if en[rank] in blc then
            Append( s, "("); Append( s, String(en[rank]) ); Append( s, ")");
         else
            Append( s, String(en[rank]) );
         fi;
         Append(s,"\n");
         Print(s);
      fi;

   od;

end;


############################################################################################
InstallMethod( PrintObj,
   "for Satake diagram",
   true,
   [ IsSatakeDiagramOfRealForm ], 0,
function(d)

  #Print("Dynkin diagram:\n\n");
   corelg.prntdg( CartanMatrix(d), CompactSimpleRoots(d) );
   Print("\n");
   Print("Involution:  ",ThetaInvolution(d));
end);


#############################################################################################
InstallMethod( SatakeDiagram,
   "for a Lie algebra",
   true,
   [ IsLieAlgebra ], 0, 
function(L)

   local fam, tip, s, tp, d, p, u;

   fam:= NewFamily( "SatakeFam", IsSatakeDiagramOfRealForm );
   tip:= NewType( fam, IsSatakeDiagramOfRealForm and IsAttributeStoringRep );

   s  := corelg.STKDTA( L );
   tp := CartanType( s.CM );
   d  := Objectify( tip, rec(type:= tp.types) );

   SetCartanMatrix( d, s.CM );
   SetBasisOfSimpleRoots( d, s.bas );
   p:= ();
   for u in s.sym do
       p:= p*(u[1],u[2]);
   od;
   SetThetaInvolution( d, p );
   SetCompactSimpleRoots( d, s.cpt );
   return d;
end );





#################################################################################################
#
# Now the functions for computing all CSA up to conjugacy
#
#

corelg.so_sets:= function( type, n )

  local sim, k, Omega, Omega1, i, j, pls, bas, b0, b1, b2, rt, sets;

  sim:= IdentityMat( n );
  if type = "A" then
     if IsEvenInt(n) then 
        k:= n-1;
     else
        k:= n;
     fi;

     return Concatenation( [[]], List([1,3..k], i -> sim{[1,3..i]} ) );

  elif type = "B" then

     if IsEvenInt(n) then
        Omega:= [1,3..n-1];
	    Omega1:= [1,3..n-3];
     else
        Omega:= [1,3..n-2];
	    Omega1:= [1,3..n-2];
     fi;

     pls:= [ ]; # this will be the roots of the form v_i+v_{i+1}
     for i in [1..n-1] do
         rt:= List( [1..n], x -> 0 );
         rt[i]:= 1;
         for j in [i+1..n] do
             rt[j]:= 2;
         od;
         Add( pls, rt );
     od;

     sets:= [ ];
     for i in [0..Length(Omega)] do # ie construct the so sets with full indices 1,3,..,Omega[i]
         bas:= [ ];  # contains the roots v_1 - v_2, v_1+v_2,..,v_r-v_{r+1},v_r+v_{r+1},
                     # where r = Omega[i]
         for j in [1..i] do
             Add( bas, sim[Omega[j]] );
             Add( bas, pls[Omega[j]] );
         od;

         Add( sets, bas );
         for j in [i+1..Length(Omega)] do
             b0:= ShallowCopy(bas);
             for k in [i+1..j] do
                 Add( b0, sim[Omega[k]] );
             od;
             Add( sets, b0 );
         od;
     od;

     for i in [0..Length(Omega1)] do # ie construct the so sets with full indices 1,3,..,Omega[i]
         bas:= [sim[n] ];  # contains the roots v_1 - v_2, v_1+v_2,..,v_r-v_{r+1},v_r+v_{r+1},
                     # where r = Omega[i]
         for j in [1..i] do
             Add( bas, sim[Omega1[j]] );
             Add( bas, pls[Omega1[j]] );
	     od;

         Add( sets, bas );
         for j in [i+1..Length(Omega1)] do
             b0:= ShallowCopy(bas);
             for k in [i+1..j] do
                 Add( b0, sim[Omega1[k]] );
             od;
             Add( sets, b0 );
         od;
     od;
     Sort( sets, function(s,t)  return Length(s) < Length(t); end );
     return sets;
 
  elif type = "C" then

     if IsEvenInt(n) then
        Omega:= [1,3..n-1];
	b1:= n/2;
     else
        Omega:= [1,3..n-2];
	b1:= (n-1)/2;
     fi;

     pls:= [ ]; # this will be the roots of the form 2v_i
     for i in [1..n-1] do
          rt:= List( [1..n], x -> 0 );
        for j in [i..n-1] do
             rt[j]:= 2;
         od;
	rt[n]:= 1;
        Add( pls, rt );
     od;
     rt:= List( [1..n], x -> 0 );
     rt[n]:= 1;
     Add( pls, rt );
 
    sets:= [ ];
     for i in [0..b1] do # ie construct the so sets with not bad not full indices 1,3,..,Omega[i]
         bas:= [ ];  # contains the roots v_1 - v_2, v_r-v_{r+1},
                     # where r = Omega[i]
         for j in [1..i] do			# we add the roots  v_1-v_2... v_k-v_{k+1}
             Add( bas, sim[Omega[j]] );
         od;

   	     Add( sets, bas );
         for j in [2*i+1..n] do   # we add the roots  2v_r... 2 v_s
             b0:= ShallowCopy(bas);
             for k in [2*i+1..j] do
                 Add( b0, pls[k] );
             od;
             Add( sets, b0 );
         od;
     od;

     Sort( sets, function(s,t)  return Length(s) < Length(t); end );
     return sets;

  elif type = "D" then
     
     if IsEvenInt(n) then
        Omega:= [1,3..n-1];
     else
        Omega:= [1,3..n-2];
     fi;

     pls:= [ ]; # this will be the roots of the form v_i+v_{i+1}
     for i in [1..n-2] do
         rt:= List( [1..n], x -> 0 );
         rt[i]:= 1;
         for j in [i+1..n-2] do
             rt[j]:= 2;
         od;
         rt[n-1]:= 1;
         rt[n]:= 1;
         Add( pls, rt );
     od;
     rt:= List( [1..n], x -> 0 );
     rt[n]:= 1;
     Add( pls, rt );

     sets:= [ ];
     for i in [0..Length(Omega)] do # ie construct the so sets with full indices 1,3,..,Omega[i]
         bas:= [ ];  # contains the roots v_1 - v_2, v_1+v_2,..,v_r-v_{r+1},v_r+v_{r+1},
                     # where r = Omega[i]
         for j in [1..i] do
             Add( bas, sim[Omega[j]] );
             Add( bas, pls[Omega[j]] );
         od;

         Add( sets, bas );
         for j in [i+1..Length(Omega)] do
             b0:= ShallowCopy(bas);
             for k in [i+1..j] do
                 Add( b0, sim[Omega[k]] );
             od;
             Add( sets, b0 );
         od;
     od;

     # add the special one...
     if IsEvenInt(n) then
        b0:= sim{Omega};
        b0[Length(b0)]:= sim[Length(sim)];
        Add( sets, b0 );
     fi;
     Sort( sets, function(s,t)  return Length(s) < Length(t); end );
     return sets;

  elif type = "E" and n = 6 then

     return [ [  ], [ [ 1, 0, 0, 0, 0, 0 ] ], [ [ 1, 0, 0, 0, 0, 0 ], [ 0, 1, 0, 0, 0, 0 ] ], 
  [ [ 1, 0, 0, 0, 0, 0 ], [ 0, 1, 0, 0, 0, 0 ], [ 0, 0, 0, 0, 1, 0 ] ], 
  [ [ 1, 0, 0, 0, 0, 0 ], [ 0, 1, 0, 0, 0, 0 ], [ 0, 0, 0, 0, 1, 0 ], [ 1, 1, 2, 2, 1, 0 ] ] ];


  elif type = "E" and n = 7 then

      return [ [  ], [ [ 1, 0, 0, 0, 0, 0, 0 ] ], [ [ 1, 0, 0, 0, 0, 0, 0 ], 
               [ 0, 1, 0, 0, 0, 0, 0 ] ], 
     [ [ 1, 0, 0, 0, 0, 0, 0 ], [ 0, 1, 0, 0, 0, 0, 0 ], [ 0, 0, 0, 0, 1, 0, 0 ] ], 
  [ [ 1, 0, 0, 0, 0, 0, 0 ], [ 0, 1, 0, 0, 0, 0, 0 ], [ 1, 2, 2, 4, 3, 2, 1 ] ], 
  [ [ 1, 0, 0, 0, 0, 0, 0 ], [ 0, 1, 0, 0, 0, 0, 0 ], [ 0, 0, 0, 0, 1, 0, 0 ], 
[ 0, 0, 0, 0, 0, 0, 1 ] ], 
  [ [ 1, 0, 0, 0, 0, 0, 0 ], [ 0, 1, 0, 0, 0, 0, 0 ], [ 0, 0, 0, 0, 1, 0, 0 ], 
[ 1, 1, 2, 2, 1, 0, 0 ] ], 
  [ [ 1, 0, 0, 0, 0, 0, 0 ], [ 0, 1, 0, 0, 0, 0, 0 ], [ 0, 0, 0, 0, 1, 0, 0 ], 
[ 0, 0, 0, 0, 0, 0, 1 ], 
      [ 1, 1, 2, 2, 1, 0, 0 ] ], [ [ 1, 0, 0, 0, 0, 0, 0 ], [ 0, 1, 0, 0, 0, 0, 0 ], 
[ 0, 0, 0, 0, 1, 0, 0 ], 
      [ 0, 0, 0, 0, 0, 0, 1 ], [ 1, 1, 2, 2, 1, 0, 0 ], [ 1, 1, 2, 2, 2, 2, 1 ] ], 
  [ [ 1, 0, 0, 0, 0, 0, 0 ], [ 0, 1, 0, 0, 0, 0, 0 ], [ 0, 0, 0, 0, 1, 0, 0 ], 
[ 0, 0, 0, 0, 0, 0, 1 ], 
      [ 1, 1, 2, 2, 1, 0, 0 ], [ 1, 1, 2, 2, 2, 2, 1 ], [ 1, 2, 2, 4, 3, 2, 1 ] ] ];

  elif type = "E" and n = 8 then

    return
[ [  ], [ [ 1, 0, 0, 0, 0, 0, 0, 0 ] ], [ [ 1, 0, 0, 0, 0, 0, 0, 0 ], 
[ 0, 0, 0, 0, 0, 0, 0, 1 ] ], 
  [ [ 1, 0, 0, 0, 0, 0, 0, 0 ], [ 0, 0, 0, 0, 0, 0, 0, 1 ], [ 1, 1, 2, 2, 1, 0, 0, 0 ] ], 
  [ [ 1, 0, 0, 0, 0, 0, 0, 0 ], [ 0, 0, 0, 0, 0, 0, 0, 1 ], [ 1, 1, 2, 2, 1, 0, 0, 0 ], 
[ 0, 1, 0, 0, 0, 0, 0, 0 ] ], 
  [ [ 1, 0, 0, 0, 0, 0, 0, 0 ], [ 0, 0, 0, 0, 0, 0, 0, 1 ], [ 1, 1, 2, 2, 1, 0, 0, 0 ], 
[ 2, 3, 4, 6, 5, 4, 2, 1 ] ], 
  [ [ 1, 0, 0, 0, 0, 0, 0, 0 ], [ 0, 0, 0, 0, 0, 0, 0, 1 ], [ 1, 1, 2, 2, 1, 0, 0, 0 ], 
[ 0, 1, 0, 0, 0, 0, 0, 0 ], 
      [ 0, 0, 0, 0, 1, 0, 0, 0 ] ], [ [ 1, 0, 0, 0, 0, 0, 0, 0 ], [ 0, 0, 0, 0, 0, 0, 0, 1 ], 
[ 1, 1, 2, 2, 1, 0, 0, 0 ], 
      [ 0, 1, 0, 0, 0, 0, 0, 0 ], [ 0, 0, 0, 0, 1, 0, 0, 0 ], [ 1, 1, 2, 2, 2, 2, 2, 1 ] ], 
  [ [ 1, 0, 0, 0, 0, 0, 0, 0 ], [ 0, 0, 0, 0, 0, 0, 0, 1 ], [ 1, 1, 2, 2, 1, 0, 0, 0 ], 
[ 0, 1, 0, 0, 0, 0, 0, 0 ], 
      [ 0, 0, 0, 0, 1, 0, 0, 0 ], [ 1, 1, 2, 2, 2, 2, 2, 1 ], [ 1, 2, 2, 4, 3, 2, 2, 1 ] ], 
  [ [ 1, 0, 0, 0, 0, 0, 0, 0 ], [ 0, 0, 0, 0, 0, 0, 0, 1 ], [ 1, 1, 2, 2, 1, 0, 0, 0 ], 
[ 0, 1, 0, 0, 0, 0, 0, 0 ], 
      [ 0, 0, 0, 0, 1, 0, 0, 0 ], [ 1, 1, 2, 2, 2, 2, 2, 1 ], [ 1, 2, 2, 4, 3, 2, 2, 1 ], 
[ 2, 3, 4, 6, 5, 4, 2, 1 ] ] ];

  elif type = "F" and n = 4 then

     return [ [  ], [ [ 0, 0, 1, 0 ] ], [ [ 1, 0, 0, 0 ] ], [ [ 0, 0, 1, 0 ], [ 0, 1, 2, 2 ] ], 
[ [ 1, 0, 0, 0 ], [ 1, 2, 2, 0 ] ], 
  [ [ 0, 0, 1, 0 ], [ 0, 1, 2, 2 ], [ 2, 3, 4, 2 ] ], [ [ 1, 0, 0, 0 ], [ 1, 2, 2, 0 ], 
[ 1, 2, 2, 2 ] ], 
  [ [ 1, 0, 0, 0 ], [ 1, 2, 2, 0 ], [ 1, 2, 2, 2 ], [ 1, 2, 4, 2 ] ] ];

  elif type = "G" and n = 2 then

     return [ [  ], [ [ 1, 0 ] ], [ [ 0, 1 ] ], [ [ 1, 0 ], [ 3, 2 ] ] ];

  else
    Error("no such root system");
  fi;


end;





############################################################################################
corelg.conj_func:= function( type, n )

  local R, Bil, W, wt, wts, wts1, wts2, o, conj;

  if type = "A" then

     conj:= function( A, B ) return Length(A) = Length(B); end;
     return conj;

  elif type = "B" then

     R:= RootSystem( type, n );
     W:= WeylGroup( R );
     wt:= List( [1..n], x -> 0 ); wt[1]:= 1;
     o:= WeylOrbitIterator( W, wt );
     wts:= [ ]; while not IsDoneIterator( o ) do Add( wts, NextIterator(o) ); od;
     Bil:= BilinearFormMatNF(R);
     conj:= function(A,B)

         local nmA, nmB, spA, spB;

         if Length(A) <> Length(B) then return false; fi;
         nmA:= List( A, x -> x*Bil*x ); Sort( nmA );
         nmB:= List( B, x -> x*Bil*x ); Sort( nmB );
         if nmA <> nmB then return false; fi;

         spA:= VectorSpace( Rationals, List(A,x->x*SimpleRootsAsWeights(R)) );
         spB:= VectorSpace( Rationals, List(B,x->x*SimpleRootsAsWeights(R)) );

         if Length( Filtered( wts, x -> x in spA ) ) <> 
            Length( Filtered( wts, x -> x in spB ) ) then
            return false;
         fi;
         return true;

     end;

     return conj;
 
  elif type = "C" then

     R:= RootSystem( type, n );
     Bil:= BilinearFormMatNF(R);
     conj:= function(A,B)

         local nmA, nmB, spA, spB;

         if Length(A) <> Length(B) then return false; fi;
         nmA:= List( A, x -> x*Bil*x ); Sort( nmA );
         nmB:= List( B, x -> x*Bil*x ); Sort( nmB );
         if nmA <> nmB then return false; fi;

         return true;

     end;

     return conj;

  elif type = "D" then
     
     R:= RootSystem( type, n );
     W:= WeylGroup( R );
     wt:= List( [1..n], x -> 0 ); wt[1]:= 1;
     o:= WeylOrbitIterator( W, wt );
     wts1:= [ ]; while not IsDoneIterator( o ) do Add( wts1, NextIterator(o) ); od;
     if IsEvenInt(n) then
        wt:= List( [1..n], x -> 0 ); wt[n]:= 1;
        o:= WeylOrbitIterator( W, wt );
        wts2:= [ ]; while not IsDoneIterator( o ) do Add( wts2, NextIterator(o) ); od;
     else
        wts2:= [ ];
     fi;

     conj:= function(A,B)

         local nmA, nmB, spA, spB;

         if Length(A) <> Length(B) then return false; fi;
         spA:= VectorSpace( Rationals, List(A,x->x*SimpleRootsAsWeights(R)) );
         spB:= VectorSpace( Rationals, List(B,x->x*SimpleRootsAsWeights(R)) );
         if Length( Filtered( wts1, x -> x in spA ) ) <> 
            Length( Filtered( wts1, x -> x in spB ) ) then
            return false;
         fi;
         if IsEvenInt(n) and Length(A) = n/2 then
            if Length( Filtered( wts2, x -> x in spA ) ) <> 
               Length( Filtered( wts2, x -> x in spB ) ) then
               return false;
            fi;
         fi;
         return true;

     end;

     return conj;

  elif type = "E" and n = 6 then

     conj:= function( A, B ) return Length(A) = Length(B); end;
     return conj;

  elif type = "E" and n = 7 then

     R:= RootSystem( type, n );
     W:= WeylGroup( R );
     wt:= List( [1..n], x -> 0 ); wt[n]:= 1;
     o:= WeylOrbitIterator( W, wt );
     wts1:= [ ]; while not IsDoneIterator( o ) do Add( wts1, NextIterator(o) ); od;

     conj:= function(A,B)

         local nmA, nmB, spA, spB;

         if Length(A) <> Length(B) then return false; fi;

         spA:= VectorSpace( Rationals, List(A,x->x*SimpleRootsAsWeights(R)) );
         spB:= VectorSpace( Rationals, List(B,x->x*SimpleRootsAsWeights(R)) );
         if Length( Filtered( wts1, x -> x in spA ) ) <> 
            Length( Filtered( wts1, x -> x in spB ) ) then
            return false;
         fi;
         return true;

     end;

     return conj;


  elif type = "E" and n = 8 then

     R:= RootSystem( type, n );
     wts1:= PositiveRootsNF(R);

     conj:= function(A,B)

         local nmA, nmB, spA, spB;

         if Length(A) <> Length(B) then return false; fi;
         spA:= VectorSpace( Rationals, A );
         spB:= VectorSpace( Rationals, B );
         if Length( Filtered( wts1, x -> x in spA ) ) <> 
            Length( Filtered( wts1, x -> x in spB ) ) then
            return false;
         fi;
         return true;

     end;

     return conj;

  elif type = "F" and n = 4 then

     R:= RootSystem( type, n );
     Bil:= BilinearFormMatNF(R);
     conj:= function(A,B)

         local nmA, nmB, spA, spB;

         if Length(A) <> Length(B) then return false; fi;
         nmA:= List( A, x -> x*Bil*x ); Sort( nmA );
         nmB:= List( B, x -> x*Bil*x ); Sort( nmB );
         if nmA <> nmB then return false; fi;

         return true;

     end;

     return conj;

  elif type = "G" and n = 2 then

     R:= RootSystem( type, n );
     Bil:= BilinearFormMatNF(R);
     conj:= function(A,B)

         local nmA, nmB, spA, spB;

         if Length(A) <> Length(B) then return false; fi;
         nmA:= List( A, x -> x*Bil*x ); Sort( nmA );
         nmB:= List( B, x -> x*Bil*x ); Sort( nmB );
         if nmA <> nmB then return false; fi;

         return true;

     end;

     return conj;

  else
    Error("no such root system");
  fi;


end;




#########################################################################################
corelg.SOSets:= function( C )

   local t, sim, pieces, k, l, simk, sts, sets, s, inds, done;

   t:= CartanType(C);
   sim:= C^0;

   pieces:= [ ];  # pieces[k] is a list of so sets of the k-th component
   for k in [1..Length(t.types)] do
       simk:= sim{t.enumeration[k]};
       sts:= corelg.so_sets( t.types[k][1], t.types[k][2] );
       sets:= [ ];
       for s in sts do
           Add( sets, List( s, u -> LinearCombination( u, simk ) ) );
       od;
       Add( pieces, sets );
   od;
                
   # now an so set is a union of one set from each piece...

   inds:= List( pieces, x -> 1 );
   sets:= [ ];
   done:= false;
   while not done do
      Add( sets, Union( List( [1..Length(inds)], x -> pieces[x][ inds[x] ] ) ) );
      l:= Length( inds );
      while l >= 1 and inds[l] = Length(pieces[l]) do
          inds[l]:= 1; l:= l-1;
      od;
      if l = 0 then
         done:= true;
      else
         inds[l]:= inds[l]+1;
      fi;
   od;

   Sort( sets, function(s,t)  return Length(s) < Length(t); end );
   return sets;

end;



##########################################################################################
corelg.ConjugationFct:= function( C )

   local t, k, conj, conjs;

   t:= CartanType(C);
   conjs:= [ ];
   for k in [1..Length(t.types)] do
       Add( conjs, corelg.conj_func( t.types[k][1], t.types[k][2] ) );
   od;
                
   conj:= function( A, B )

      local k, A0, B0;

      if Length(A) <> Length(B) then return false; fi;
      for k in [1..Length(t.enumeration)] do
          A0:= List( A, x -> x{ t.enumeration[k] } ); 
          A0:= Filtered( A0, x -> not IsZero(x) );
          B0:= List( B, x -> x{ t.enumeration[k] } ); 
          B0:= Filtered( B0, x -> not IsZero(x) );

          if Length(A0) = 0 then 
             if Length(B0) <> 0 then 
                return false; 
             fi;
          elif Length(B0) = 0 then
             return false;
          else
             if not conjs[k]( A0, B0 ) then return false; fi;
          fi;
      od;
      return true;
   end;

   return conj;

end;



##############################################################################################
InstallMethod( CartanSubalgebrasOfRealForm,
   "for a Lie algebra",
   true,
   [ IsLieAlgebra ], 0, function( L )

   local cd, H, C, R, ch, p, pr, i, j, s, sets, so, sums, sim, B, CM, f, bc, kappa, cfi, cfj, 
         CSAs, hh, mat, sol, pos, M, r, Kappa, sp, b;
  
   cd:= CartanDecomposition(L);
   H:= MaximallyNonCompactCartanSubalgebra( L );
   C:= Intersection( H, cd.P );
   R:= RootsystemOfCartanSubalgebra( L, H );
   ch:= ChevalleyBasis(R);
   p:= PositiveRootsNF(R);
   pr:= [ ];
   for i in [1..Length(p)] do
       if ch[1][i]*ch[2][i] in C then
          Add( pr, p[i] );
       fi;
   od;
   if Length(pr) = 0 then return [ H ]; fi;

   sums:= [ ];
   for i in [1..Length(pr)] do
       for j in [i+1..Length(pr)] do 
           Add( sums, pr[i]+pr[j] );
       od;
   od;
   sim:= Filtered( pr, x -> not x in sums );  

   B:= BilinearFormMatNF(R);
   CM:= List( sim, x -> List( sim, y -> 2*x*B*y/(y*B*y) ) );
   
   so:= corelg.SOSets(CM);
   so:= List( so, x -> List( x, y -> y*sim ) );

   f:= corelg.ConjugationFct( CartanMatrix(R) );
   sets:= [so[1]];
   for i in [2..Length(so)] do
       if ForAll( sets, x -> not f(x,so[i]) ) then
          Add( sets, so[i] );
       fi;
   od;

   bc:= BasisVectors( Basis(C) );
   kappa:= List( bc, x -> [ ] );
   Kappa:= KillingMatrix( Basis(L) );
   for i in [1..Length(bc)] do
       for j in [i..Length(bc)] do
           cfi:= Coefficients( Basis(L), bc[i] );
           if i = j then
              kappa[i][i]:= cfi*Kappa*cfi;
           else
              cfj:= Coefficients( Basis(L), bc[j] );
              kappa[i][j]:= cfi*Kappa*cfj;
              kappa[j][i]:= kappa[i][j];
           fi;
       od;
   od; 

   CSAs:= [ ];
   for s in sets do 

       if Length(s) = 0 then
          Add( CSAs, H );
       else
          hh:= [ ];
          for r in s do
              pos:= Position( p, r );
              Add( hh, Coefficients( Basis(C), ch[1][pos]*ch[2][pos] ) );
          od;
          mat:= TransposedMat( hh*kappa );
          sol:= NullspaceMat( mat );
          hh:= List( sol, x -> x*Basis(C) );

          while Length(hh) < Dimension(H) do
               sp:= MutableBasis( LeftActingDomain(L), hh, Zero(L) );
               b:= Filtered( Basis(cd.K), x -> ForAll( hh, y -> IsZero(x*y) ) );
               pos:= PositionProperty( b, x -> not IsContainedInSpan( sp, x ) );
               if pos <> fail then
                  Add( hh, b[pos] );
                  CloseMutableBasis( sp, b[pos] );
               else
                  M:= Intersection( cd.K, LieCentralizer( L, Subalgebra( L, hh ) ) );
                  b:= BasisVectors( Basis( CartanSubalgebra( M ) ) );
                  i:= 1;
                  while Length(hh) < Dimension(H) do
                      if not IsContainedInSpan(sp,b[i]) then
                         Add(hh,b[i]);
                         CloseMutableBasis( sp, b[i] );
                      fi;
                      i:= i+1;
                  od;

               fi;
          od;
          Add( CSAs, SubalgebraNC( L, hh ) );
       fi;
   od;
   return CSAs;       

end );

corelg.namesimple:= function( id )    
    
local nr, mv, nr1, nr2, en, q, p, t, r;

    t:= id[1]; r:= id[2]; q:= id[3];

    if t = "A" and r=1 then
       if q=0 then return "sl(2,C)"; fi;
       if q=1 then return "su(2)"; fi;
       if q=2 then return "sl(2,R)"; fi;
    elif t = "A" and r>1 then
       if q=0 then return Concatenation( "sl(", String(r+1),",C)"); fi;
       if q=1 then return Concatenation( "su(", String(r+1),")"); fi;
        nr := First([1..r],x-> x>= r/2);
        mv := 1; if IsOddInt(r) then mv := 2; fi;
        p:= Position( [2..nr+1], q );
        if p <> fail then 
           return Concatenation( "su(", String(p),",",String(r+1-p),")");
        fi;
        if mv = 1 and q = nr+2 then
           return Concatenation( "sl(", String(r+1),",R)");
        fi;
        if mv = 2 and q = nr+2 then 
           return Concatenation("sl(",String((r+1)/2),",H)");
        fi;
        if mv=2 and q = nr+3 then
           return Concatenation( "sl(",String(r+1),",R)"); 
        fi;
    elif t = "B" then
        if q=0 then return Concatenation( "so(", String(2*r+1),",C)"); fi;
        if q=1 then return Concatenation( "so(", String(2*r+1),")"); fi;
        p:= Position( [2..r+1], q );
        return Concatenation("so(",String(2*p),",",String(2*r-(2*p)+1),")"); 
    elif t = "C" then
       if q=0 then return Concatenation( "sp(", String(2*r),",C)"); fi;
       if q=1 then return Concatenation( "sp(", String(r),")"); fi;    
       nr := First([1..r],x-> x> r/2)+1;
       p:= Position( [2..nr-1], q );
       if p <> fail then 
          return Concatenation( "sp(",String(p),",",String(r-p),")");
       fi;
       if q = nr then 
          return Concatenation("sp(",String(r),",R)");
       fi;
    elif t = "D" then
       nr1 := First([1..r+1],x-> x> r/2)-1;
       nr2 := First([1..r+1],x-> x> (r-1)/2)-1;
       if r = 4 then nr1 := nr1-1; fi;
       if q=0 then return Concatenation( "so(", String(2*r),",C)"); fi;
       if q=1 then return Concatenation( "so(", String(2*r),")"); fi;
       if r > 4 then
          p:= Position( [2..nr1+1], q );
          if p <> fail then
             return Concatenation("so(",String(2*p),",",String(2*r-2*p),")");
          fi;
          if q = nr1+2 then
             return Concatenation("so*(",String(2*r),")");
          fi;
          if q = nr1+3 then
             return Concatenation( "so(",String(2*r-1),",1)");
          fi;
          p:= Position( [nr1+4..nr1+r+nr2], q );
          if p <> fail then
             return Concatenation("so(",String(2*p+1),",",String(2*r-2*p-1),")");
          fi;
       else
          if q=2 then
             return "so*(8)";
          elif q=3 then
             return "so(4,4)";
          elif q=4 then
             return "so(3,5)";
          elif q=5 then
             return "so(1,7)";
          fi;
       fi;
       
    elif t = "G" then
       if q=0 then return "G2(C)"; fi;
       if q=1 then return "G2c"; fi;
       if q=2 then return "G2(2)"; fi;    
    elif t = "F" then
       if q=0 then return "F4(C)"; fi;
       if q=1 then return "F4c"; fi;
       if q=2 then return "F4(4)"; fi;
       if q=3 then return "F4(-20)"; fi;
    elif t = "E" then
       if q=0 then return Concatenation("E",String(r),"(C)"); fi;
       if q=1 then return Concatenation("E",String(r),"c"); fi;
       if r = 6 then
          if q=2 then return "E6(6)"; fi;
          if q=3 then return "E6(2)"; fi;
          if q=4 then return "E6(-14)"; fi;
          if q=5 then return "E6(-26)"; fi;
       elif r=7 then
          if q=2 then return "E7(7)"; fi;
          if q=3 then return "E7(-25)"; fi;
          if q=4 then return "E7(-5)"; fi;
       elif r=8 then
          if q=2 then return "E8(8)"; fi;
          if q=3 then return "E8(-24)"; fi;
       fi;
    fi;
    
end;

InstallMethod( NameRealForm,
   "for a Lie algebra",
   true,
   [ IsLieAlgebra ], 0, function( L )


        local C, L0, cd, H, v, id, s, i, k, p;

# we assume that L is reductive!

        if IsBound(L!.sstypes) then
	   id:= L!.sstypes;
	   s:= "";
           for i in [1..Length(id)] do
               s:= Concatenation( s, corelg.namesimple( id[i] ) );
               if i < Length( id ) then
                  s:= Concatenation( s, "+" );
               fi;
           od;
	   return s;
        fi;	   

        C:= LieCentre(L);
        if Dimension(C) = 0 then
           L0:= L;
        else
           L0:= LieDerivedSubalgebra(L);
           cd:= CartanDecomposition(L);
           SetCartanDecomposition( L0, rec( CartanInv:= cd.CartanInv,
                K:= Intersection( L0, cd.K), P:= Intersection( L0, cd.P ) ) );              if HasMaximallyCompactCartanSubalgebra( L ) then
               H:= MaximallyCompactCartanSubalgebra( L );
               SetMaximallyCompactCartanSubalgebra( L0,
                    Intersection( L0, H ) ); 
            fi;
        fi;
        v:= VoganDiagram(L0);
        id:= v!.sstypes;
        s:= "";
        for i in [1..Length(id)] do
            s:= Concatenation( s, corelg.namesimple( id[i] ) );
            if i < Length( id ) then
               s:= Concatenation( s, "+" );
            fi;
        od;

        if Dimension(C) > 0 then
           k:= Dimension( Intersection( C, cd.K ) );
           p:= Dimension( Intersection( C, cd.P ) );
           s:= Concatenation( s, " + a torus of ");
           if k > 0 then 
              s:= Concatenation( s, String( k ), " compact dimensions" );
           fi;
           if p > 0 then
              if k > 0 then s:= Concatenation( s, " and" ); fi;
              s:= Concatenation( s, " ", String(p), " non-compact dimensions");
           fi;
        fi;
        return s;
end );

InstallMethod( MaximalReductiveSubalgebras,
 "for type, rank and number", 
 true, [ IsString, IsInt, IsInt ], 0,   

function( type, rk, no ) 

      local L, b1, b2, c, u, i, F, K, cd, C, H, list, subs, l;

      if not rk in [1..8] then
         Error("The maximal reductive subalgebras are available for simple real Lie algebras of rank up to 8");
      fi;

      if not no in [1..NumberRealForms(type,rk)] then
         Error("There is no real form with the input parameters");
      fi;

      F:= SqrtField;
      L:= RealFormById( type, rk, no, F );

      list:= Filtered( corelg.Linc, x -> x[1] = type and x[2] = rk and x[3] = no );
      subs:= [ ];
      for l in list do 

          b1:= [ ];
          for c in l[4] do
              u:= Zero(L);
              for i in [1,3..Length(c)-1] do
                  u:= u + (c[i+1]*One(F))*Basis(L)[c[i]];
              od;
              Add( b1, u );
          od;

          b2:= [ ];
          for c in l[5] do
              u:= Zero(L);
              for i in [1,3..Length(c)-1] do
                  u:= u + (c[i+1]*One(F))*Basis(L)[c[i]];
              od;
              Add( b2, u );
          od;

          K:= Subalgebra( L, b1, "basis" );
          if Length(b2) < rk then
             C:= LieCentralizer( L, K );
             if Dimension(C) > 0 then
                Append( b1, BasisVectors( Basis(C) ) );
                Append( b2, BasisVectors( Basis(C) ) );
                K:= Subalgebra( L, b1, "basis" );
             fi;
          fi;

          cd:= CartanDecomposition(L);
          SetCartanDecomposition( K, rec( CartanInv:= cd.CartanInv,
                K:= Intersection( K, cd.K), P:= Intersection( K, cd.P ) ) );              H:= Subalgebra( L, b2, "basis" );
          SetMaximallyCompactCartanSubalgebra( K, H );      

          Add( subs, K );
      od;

      return rec( liealg:= L, subalgs:= subs );

end );
