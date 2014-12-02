########################################################################################
#  
#  this file contains the functions to read the nilpotent orbits from the database;
#  also the semi-automated functions for creating the database are listed here (but not
#  officially documented).
#
#  this file contains the following functions
#     NilpotentOrbitsOfRealForm
#     CarrierAlgebraOfNilpotentOrbit
#     corelg.readDBCA
#     corelg.readDBTriples
#     corelg.SL2tripleOfNilpotentElement
#     corelg.SL2tripleOfCharacteristic
#     corelg.SqrtEltMySign
#     corelg.CayleyTransform
#     corelg.CayleyTransformInverse
#     corelg.ChevalleySystemInnerType
#     corelg.checkTriples
#     corelg.mySLAfctCanBas
#     corelg.principalOrbitsOfRealForm
#     corelg.lookupRealCayleyTriple
#     corelg.RealCayleyTriplesOfRealForm
#     corelg.ConvertRealCayleyTriplesToNilpotentOrbits
#     corelg.viewReducedEquationsAndAttach
#     corelg.attachSolution
#     corelg.TryToFindComplexCayleyTriple
#     corelg.calgDBentries
#     corelg.WriteRealNilpotentOrbitsToDB
#     corelg.RealNilpotentOrbitsInDatabase
#     corelg.RealNilpotentOrbitsFromDatabase


#######################################################################################
## NOT OFFICIALLY DOCUMENTED:
##
## How to construct real nilpotent orbits using only the database of carrier algs:
##
## 1) Let form be an entry of corelg.NonCompactRealFormsOfSimpleLieAlgebra
##
## 2) res := corelg.RealCayleyTriplesOfRealForm(form)
##    now res is a record with entries
##    - form: the input form
##    - triples: a record with entries
##          - principal: contains records with entries
##                 - realsl2: real Cayley triple [f,h,e] with principal carrier alg
##                 - cdims  : dimensions of gradations of carrier algebra
##          - nonprincipal: contains records with entries
##                 - realsl2: real Cayley triple [f,h,e] with non-principal carrier alg
##                 - cdims  : dimensions of gradations of carrier algebra
##            OR
##                 - oldsl2: the original homogeneous real triple [f,h,e]
##                 - carrier: the corresponding carrier alg with entries g0, gp, gn
##                 - cdims : dimensions of gradation of carrier algebra
##    - tobedone: the indices i such that res.triples.nonprincipal[i] does NOT have
##                a realsl2 attached. These entries have to be solved as below
##
## 3) if res.tobedone=[], then all real triples were constructed, go to 5)
##
## 4) for an entry i in res.tobedone do
##      corelg.TryToFindComplexCayleyTriple(res, i, nrvars, nrtries)
##    here nrvars and nrtries are lists of integers.
##    IMPORTANT: name of variable 'res' MUST BE "res"
##               in order to save solution automatically
##    This function sets up equations to find a complex Cayley triple; 
##    Then all but nrvars[j] variables are set to zero and this is tried nrtries[j]
##    times until we find a system of equations which has a non-trivial Groebner basis.
##    If a solution can be generated automatically (which should almost always be the case),
##    then it is attached. Otherwise, one has to manually follow the instructions; or to do
##    the same again. Start with nrvars := [6..14] (or so).
##    nrtries can also be a single integers; then for every entry nrvars[j] the same number
##    of tries is used.
##    If a solution has been found, then the other open cases in res.triples.nonprincipal will
##    be checked and, if possible, solved automatically. Also, the solution will be written
##    to the database calg_db in the file carrierAlg.db
##    Afterwards, res.triples.tobedone will be updated.
##
## 5) If res.triples.tobedone = [], then do
##      orbs := corelg.ConvertRealCayleyTriplesToNilpotentOrbits(res)
##    This function converts the real triples into objects "NilpotentOrbit"
##    and computes the corresponding WDDs of the characteristics
##
##
## Use "corelg.WriteRealNilpotentOrbitsToDB(type,rank)" for the semi-automated version of this algorithm
## 
##
############################################################################################


#################################################################################
#################################################################################
#
# READ DATABASES
#
#################################################################################
#################################################################################

## DATABASE OF CARRIER ALGEBRAS
corelg.readDBCA := function()
   Print("#I CoReLG: read database of carrier algebras");
   ReadPackage( "corelg", "gap/carrierAlg.db" );
   Print(" ... done\n");
end;
if not IsBound(corelg.carrierAlgDB) then corelg.carrierAlgDB := []; fi;


## DATABASE OF REAL TRIPLES
corelg.readDBTriples := function()
   Print("#I CoReLG: read database of real triples");
   ReadPackage( "corelg", "gap/realTriples.db" );
   Print(" ... done\n");
end;

if not IsBound(corelg.realtriplesDB) then corelg.realtriplesDB := []; fi;

## IS SINGULAR LOADED?
if not IsBound(HasTrivialGroebnerBasis) then HasTrivialGroebnerBasis:=function()end; fi;



#################################################################################
#################################################################################
#
# SOME PRELIMINARY FUNCTIONS
#
#################################################################################
#################################################################################


################################################
#input:  liealg L and nilpotent element x in L
#output: an SL2-triple (f,h,x)
#remark: just a modification of SL2Triple,
#        avoid nilpotency test
################################################
corelg.SL2tripleOfNilpotentElement := function ( L, x )
local  n, F, B, xc, eqs, T, i, j, k, l, cij, b, v, z, h, R, BR, Rvecs, 
       H, e0, e1, y;
   n := Dimension( L );
   F := LeftActingDomain( L );
   B := Basis( L );
   T := StructureConstantsTable(B);
   xc := Coefficients( B, x );
   eqs := NullMat( 2 * n, 2 * n, F );
   for i  in [ 1 .. n ]  do
      for j  in [ 1 .. n ]  do
         cij := T[i][j];
         for k  in [ 1 .. Length( cij[1] ) ]  do
            l := cij[1][k];
            eqs[i][l] := eqs[i][l] + xc[j] * cij[2][k];
            eqs[n + i][n + l] := eqs[n + i][n + l] + xc[j] * cij[2][k];
         od;
      od;
      eqs[n + i][i] := One( F );
   od;
   b := [  ];
   for i  in [ 1 .. n ]  do
      b[i] := Zero( F );
      b[n + i] := 2 * One( F ) * xc[i];
   od;
   v := SolutionMat( eqs, b );
   if v = fail  then
      return fail;
   fi;
   z := LinearCombination( B, v{[ 1 .. n ]} );
   h := LinearCombination( B, v{[ n + 1 .. 2 * n ]} );
   R := LieCentralizer( L, SubalgebraNC( L, [ x ],"basis" ) );
   BR := Basis( R );
   Rvecs := BasisVectors( BR );
   H := List( Rvecs, function ( v )
           return Coefficients( BR, h * v );
       end );
   H := H + 2 * IdentityMat( Dimension( R ), F );
   e0 := Coefficients( BR, h * z + 2 * z );
   e1 := SolutionMat( H, e0 );
   if e1 = fail  then
      return fail;
   fi;
   y := z - LinearCombination( Rvecs, e1 );
   return [ y, h, x ];
end;


###################################################################
#input:  liealg, grading record gr with entries g0, gp, gn, and 
#        a characteristic h in span of g0
#output: an SL2-triple (f,h,x) with x in span of gp[1]
###################################################################
corelg.SL2tripleOfCharacteristic := function( L, gr, h )
local e, f, co, x, sp, mat, sol;
   e := gr.gp[1];
   f := gr.gn[1];
   while true do
      co := List( e, x -> Random([-2..2]) );
      x  := co*e;
      sp := SubspaceNC( L, List( f, y -> x*y) );
      if Dimension(sp) = Length(e) and h in sp then
         mat := List( f, u -> Coefficients( Basis(sp), x*u ) );
         sol := SolutionMat( mat, Coefficients( Basis(sp), h ) );
         return [sol*f,h,x];   
      fi;
   od;
end;


################################################
#input:  rational or real SqrtFieldElt with one monom
#output: the sign of v
################################################
corelg.SqrtEltMySign :=function(v) 
   if v=0*v then return 0; fi;
   if IsSqrtFieldElement(v) then
      if IsPosSqrtFieldElt(v) then return 1; else return -1; fi; 
   fi;
   if v>0 then return 1; else return -1; fi;
end;


################################################
#Input:  tr: complex Cayley triple [f,h,e]
#Output: real Cayley triple (cayley transfom)
################################################
corelg.CayleyTransform := function(tr,F)
local e, f, h,i;
   f := tr[1]; 
   h := tr[2]; 
   e := tr[3]; 
   i := E(4)*One(F);
   return [(i/2)*(e-f+h),e+f,(i/2)*(e-f-h)];
end;

################################################
#Input:  tr: real Cayley triple [f,h,e]
#Output: complex Cayley triple (cayley transfom)
################################################
corelg.CayleyTransformInverse := function(tr,F)
local e, f, h,i;
   f := tr[1]; 
   h := tr[2]; 
   e := tr[3]; 
   i := E(4)*One(F);
   return [(-1/2*One(F))*(-h-i*f-i*e),i*(e-f),(1/2*One(F))*(h-i*f-i*e)];
end;


################################################
#Input:  calg:  carrier  subalgebra (of L)
#        tr  :  record with entries g0, gp, gn
#        sigma: complex conjugation in L  
#Output: Chevalley system in form of list
#          [ [x_alpha], [x_{-alpha}], [h_\alpha]]
#Remark: Need that CSA of L lie in K where L=K+P
################################################
corelg.ChevalleySystemInnerType := function(calg, tr, sigma)
local rs, cgen, cf1, cf2, cf, im, pv, nv, hs, i, fstCB, sndCB, cgen2;

  Info(InfoCorelg,5,"    start ChevalleySystemInnerType in dimension ",Dimension(calg));
  #construct carrier algebra with roots system and 
  #canonical generators
   rs := RootSystemOfZGradedLieAlgebra(calg,tr);
   SetRootSystem(calg,rs);

  #get a (non-special) Chevalley system
   cgen    := List(CanonicalGenerators(rs),x->List(x,y->y));
   cgen[2] := -cgen[2];

  #find scalars to construct a special system
   cf1 := List(List(cgen[1],sigma), 
               x->First(Coefficients(Basis(calg),x),y->not y=0*y));
   cf2 := List(cgen[2],             
               x->First(Coefficients(Basis(calg),x),y->not y=0*y));
   
   if IsSqrtField(LeftActingDomain(calg)) then
      cf  := List([1..Length(cf1)],i->Sqroot(AbsoluteValue(cf1[i]^-1*cf2[i])));
   else
      cf  := List([1..Length(cf1)],i->Sqrt(AbsoluteValue(cf1[i]^-1*cf2[i])));
   fi;   
   
  #now set up automorphism and construct special Chevalley system
   im  := StructuralCopy(cgen);
   for i in [1..Length(cgen[1])] do
      im[1][i] := im[1][i]*cf[i];
      im[2][i] := im[2][i]*cf[i]^-1;
   od;
 
  #mapping fstCB to sndCB is an automorphism
   fstCB := List(Flat(SLAfcts.canbas( calg, [cgen[1],-cgen[2],cgen[3]])),
                      x-> Coefficients(Basis(calg),x));
   sndCB := Flat(SLAfcts.canbas( calg, [im[1],-im[2],im[3]]));
   pv      := List(ChevalleyBasis(calg)[1],x->   
              SolutionMat(fstCB,Coefficients(Basis(calg),x))*sndCB);
   nv      := List(ChevalleyBasis(calg)[2],x->   
              SolutionMat(fstCB,Coefficients(Basis(calg),-x))*sndCB);
   hs      := List([1..Length(pv)],x-> -pv[x]*nv[x]);
   hs      := Concatenation(hs,List([1..Length(nv)],x-> -nv[x]*pv[x]));

   Info(InfoCorelg,5,"    end ChevalleySystemInnerType");
   return rec(chevSys := [pv,nv,hs],
              rank    := Length(SimpleSystem(rs)));
end;




################################################
#Output: true, if all rsl2 are real Cayley triples
################################################
corelg.checkTriples := function(form, triples)
local L, K, P, theta, sigma, tr, f, e, h;
   Info(InfoCorelg,5,"    start corelg.checkTriples");
   L     := form.liealgSF;
   theta := CartanDecomposition(L).CartanInv;
   sigma := RealStructure(L);
   for tr in triples do
      if IsBound(tr.realsl2) then
         f := tr.realsl2[1]; 
         h := tr.realsl2[2]; 
         e := tr.realsl2[3];
         if not (h*f = -2*f and h*e=2*e and e*f=h) or
            not theta(e)=-f or not tr.realsl2 = List(tr.realsl2,sigma) then
               Error("not a real sl2 triple");
         fi;
      fi;   
   od; 
   if Length(Filtered(triples,x->not IsBound(x.realsl2)))>0 then
      Print("corelg.checkTriples: there are entries withouth attached realsl2\n");
   fi;
   Info(InfoCorelg,5,"    end corelg.checkTriples");
   return true;
end;



################################################
#Input:  slightly mod version of SLAfct.canbas
################################################
corelg.mySLAfctCanBas := function ( L, c )
    local  x, y, x1, y1, done, levelx, levely, newlevx, newlevy, sp, i, j, u, tmp;
    x := c[1];
    y := c[2];
    x1 := ShallowCopy( x );
    y1 := ShallowCopy( y );
    done := false;
    levelx := ShallowCopy( x );
    levely := ShallowCopy( y );
    while not done  do
        newlevx := [  ];
        newlevy := [  ];
        sp  := MutableBasis( SqrtField, [  ], Zero(SqrtField)*c[1][1]);
        for i  in [ 1 .. Length( x ) ]  do
            for j  in [ 1 .. Length( levelx ) ]  do
                u := x[i] * levelx[j];         
                if not IsZero( u ) and not u in sp then 
                   # corelg.eltInSubspace(L,BasisVectors(sp), u) then
                    Add( newlevx, u );
                    CloseMutableBasis( sp, u );
                    u   := y[i] * levely[j];
                    Add( newlevy, u );
                fi;
            od;
        od;
        if newlevx <> [  ]  then
            Append( x1, newlevx );
            Append( y1, newlevy );
            levelx := newlevx;
            levely := newlevy;
        else
            done := true;
        fi;
    od;
    return [ x1, y1, c[3] ];
end;




#################################################################################
#################################################################################
#
# FUNCTIONS TO CONSTRUCT REAL TRIPLES FROM SCRATCH (ONLY PRINC CAlgs)
#
#################################################################################
#################################################################################


################################################
#Input:  form: a real form of a lie algebra
#              record containing liealg and grading
#Output: record with the following entries
#            - principal: record with
#                   -oldsl2    (old triple of K in p) 
#                   -cayleysl2 (complex Cayley triple)
#                   -realsl2   (corresp. real Cayley triple)
#            - nonprincipal: record with
#                    -oldsl2    (old triple of K in p)
#                    -carrier    rec with g0, gp, gn
################################################
corelg.principalOrbitsOfRealForm := function(form)
local res, L, K, P, sl2, sigma, tr, f, h, e, calg, cs, writeToSF,
      r, mat, ns, cf, x, rsl2, ll, tt1, tt2, F, r2, tmp, t, g,
      LSF, sigmaSF, csl2, T;  

   Info(InfoCorelg,4,"start corelg.principalOrbitsOfRealForm");
   res := rec(principal :=[], nonprincipal :=[]); 
   L   := form.liealg;
   F   := LeftActingDomain(L);
   K   := Basis(CartanDecomposition(L).K);
   P   := Basis(CartanDecomposition(L).P);
   g   := [K,P];
   writeToSF := form.writeToSF;
   LSF := form.liealgSF;

   if Dimension(CartanSubalgebra(CartanDecomposition(L).K)) 
      = Dimension(CartanSubalgebra(L)) then
      sl2 := SLAfcts.nil_orbs_inner( L, g[1], g[2], g[2] );;
   else
    
      sl2 := corelg.nil_orbs_outer(L, g[1], g[2], g[2] );;
   fi;
   for t in [1..Length(sl2.sl2)] do
      tmp        := RegularCarrierAlgebraOfSL2Triple( L, sl2.sl2[t] );
      sl2.sl2[t] := rec(sl2  := sl2.sl2[t], 
                        g0   := tmp.g0, 
                        gp   := tmp.gp, 
                        gn   := tmp.gn);
   od;

  #complex conjugation
   sigma   := RealStructure(L);
   sigmaSF := RealStructure(LSF);

   for t in [1..Length(sl2.sl2)] do
      tr := sl2.sl2[t]; 
     #tr has the form rec(sl2=rec(f,h,e), g0,gp,gn) 
     #where the g0, gp, gn describe a carrier alg
      Info(InfoCorelg,4,"   consider triple ",t," of ",Length(sl2.sl2));
      f    := tr.sl2[1]; 
      h    := tr.sl2[2];
      e    := tr.sl2[3];
      calg := SubalgebraNC(L,corelg.myflat(Concatenation(tr.g0,tr.gp,tr.gn)),"basis");
     tr := sl2.sl2[t]; 
     #tr has the form rec(sl2=rec(f,h,e), g0,gp,gn) 
     #where the g0, gp, gn describe a carrier alg
      Info(InfoCorelg,4,"   consider triple ",t," of ",Length(sl2.sl2));
      f    := tr.sl2[1]; 
      h    := tr.sl2[2];
      e    := tr.sl2[3];
      calg := SubalgebraNC(L,corelg.myflat(Concatenation(tr.g0,tr.gp,tr.gn)),"basis");

     #do we have principal case?
      if IsAbelian(SubalgebraNC(L,tr.g0)) then
        #set CSA
        #basis of rootsystem lies in tr.gp[1]    
         SetCartanSubalgebra(calg,SubalgebraNC(calg,tr.g0));
         cs   := corelg.ChevalleySystemInnerType(calg,tr,sigma);
         r    := [1..Length(SimpleSystem(RootSystem(calg)))];
      
         mat  := List(Concatenation(cs.chevSys[3]{r},[-h]),
                      x->Coefficients(Basis(L),x));
         ns   := NullspaceMat(mat)[1];
         cf   := List(ns{[1..Length(r)]}, x->Sqroot(x));
         x    := Sum(List([1..Length(r)], i-> cf[i]*(writeToSF(Flat(cs.chevSys)[r[i]])))); 
         csl2 := [sigmaSF(x),writeToSF(h),x];
         rsl2 := corelg.CayleyTransform(csl2,SqrtField); 

         Add(res.principal, rec(cdims     := [[Length(tr.g0)],List(tr.gp,Length),
                                                            List(tr.gn,Length)],  
                               #cayleysl2 := csl2,
                                realsl2   := rsl2));
      else

         Add(res.nonprincipal, rec(oldsl2  := tr.sl2,
                               cdims   := [[Length(tr.g0)],List(tr.gp,Length),
                                                            List(tr.gn,Length)],
                               carrier := rec(g0  := tr.g0,
                                              gp  := tr.gp,
                                              gn  := tr.gn)));                  
      fi;
   od; 
   if not corelg.checkTriples(form,res.principal) then Print("dammit!"); fi;
   Info(InfoCorelg,4,"end corelg.principalOrbitsOfRealForm");
   return res;
end;



################################################
#Input:  form and an entry of nonprincipal out 
#        containing entries oldsl2, cdims, carrier
#Output: try to look-up a real Cayley triple for oldsl2
#        in the database and stores it if possible;
#        otherwise returns false
################################################
corelg.lookupRealCayleyTriple := function(form, out)
local L, sigma, ca, grad, rs, cm, pos, ct, cg, salgs, sl2, i, j, t,enum, F,
      newcg, new, ngrad, ndims, cf, csl2, canbas, my, cb, myenum, dbenum, l1, l2,
      mycg, mygr, mycf, l, mys1, mys1c, newtr, dbcf, db, cand, tmp, perm, ins1, mycf2,
      bas, dim, cfh, bas2, dim2,tr,cs,r,mat,ns,x, h, y, xy, K, P,s, my2, mycg2, signs,
      writeToSF, LSF, sigmaSF, mycgSF, getDBpositions;

   if Length(corelg.carrierAlgDB)=0 then corelg.readDBCA(); fi;

  #find entry in DB corelg.carrierAlgDB
   getDBpositions := function(type, rank, dims,ins1,signs)
   return 
       Filtered([1..Length(corelg.carrierAlgDB)],x->
             corelg.carrierAlgDB[x].type = type and 
             corelg.carrierAlgDB[x].rank = rank and 
             corelg.carrierAlgDB[x].dims = dims and 
             corelg.carrierAlgDB[x].ins1 = ins1 and
             corelg.carrierAlgDB[x].ordsigns =signs);
   end;

   L     := form.liealg;
   LSF   := form.liealgSF;
   writeToSF := form.writeToSF;
   sigma   := RealStructure(L);
   sigmaSF := RealStructure(LSF);

   ca    := SubalgebraNC(L,corelg.myflat(Concatenation(out.carrier.g0,
                                    out.carrier.gp,out.carrier.gn)),"basis"); 
   SetCartanSubalgebra(ca,Intersection(ca,CartanSubalgebra(L)));
   grad  := List([ [out.carrier.g0], out.carrier.gp, out.carrier.gn],
                 x-> List(x,y->SubspaceNC(L,y,"basis")));

  #split into simple parts
   rs    := RootSystemOfZGradedLieAlgebra(ca,out.carrier);
   SetRootSystem(ca,rs);
   cm    := CartanMatrix(rs);
   ct    := CartanType(cm);
   cg    := CanonicalGenerators(rs);
   cb    := corelg.myflat(SLAfcts.canbas( ca, cg));
   salgs := [];
   sl2   := out.oldsl2;
   for i in [1..Length(ct.types)] do
      t     := ct.enumeration[i];
      newcg := List(cg,x->x{t});
      enum  := [1..Length(t)];
      new   := rec(type := ct.types[i],
                   enum := enum,
                   alg  := SubalgebraNC(ca,corelg.myflat(newcg),"basis"));
      ngrad := List(grad, x-> List(x,y->Intersection(y,new.alg)));
      ngrad := List(ngrad,x-> Filtered(x,y->Dimension(y)>0));
      ndims := List(ngrad,x-> List(x,Dimension));
      new.grad := ngrad;
      new.dims := [ndims[1],ndims[2],ndims[3]];
      new.cgen := newcg;
      new.isprincipal := IsAbelian(SubalgebraNC(new.alg,Basis(ngrad[1][1]),"basis"));
      Add(salgs,new);
   od;

  #have to split characteristic for principle orbits?
   if Length(Filtered(salgs,x->x.isprincipal))>0 then
      bas  := List(salgs,x->BasisVectors(CanonicalBasis(x.alg)));
      dim  := List(bas,Length);
      bas2 := Basis(VectorSpace(CF(4),corelg.myflat(bas),"basis"),corelg.myflat(bas));
      cfh  := Coefficients(bas2,sl2[2]);
      dim2 := [[1..dim[1]]];
      for i in [2..Length(dim)] do
         dim2[i] := dim2[i-1][Length(dim2[i-1])]+[1..dim[i]];
      od;
      cfh  := List(dim2,x->cfh{x});
      for i in [1..Length(salgs)] do salgs[i].h := cfh[i]*bas[i]; od;
   fi;

   if Length(ct.types) > 1 then
      Info(InfoCorelg,4,"  Carrier algebra has decomposition",ct.types);
   fi;

   newtr    := [];
   for i in [1..Length(salgs)] do
      my     := salgs[i].alg;
      mycg   := salgs[i].cgen;
      mygr   := salgs[i].grad;

      if not salgs[i].isprincipal then
         myenum := salgs[i].enum;
        #order gens st first gens lie in s0
         l      := Length(mycg[1]);
        #mys1   := Filtered([1..l], x-> corelg.eltInSubspace(L,Basis(mygr[2][1]),mycg[1][x]));
         mys1   := Filtered([1..l], x-> mycg[1][x] in mygr[2][1] );

         mys1c  := Filtered([1..l], x-> not x in mys1);
         mycg   := List(mycg, x-> Concatenation(x{mys1},x{mys1c}));
         perm   := PermList(Concatenation(mys1,mys1c))^-1;
         myenum := List(myenum, x-> x^perm);
         ins1   := AsSortedList(List([1..Length(mys1)],x->Position(myenum,x)));
         mycf   := List([1..l], x->
                   Coefficients(Basis(VectorSpace(CF(4),[mycg[2][x]],"basis"),
                               [mycg[2][x]]),sigma(mycg[1][x]))[1]);
         signs  := List(mycf{myenum},corelg.SqrtEltMySign);
         pos    := getDBpositions(salgs[i].type[1],salgs[i].type[2], 
                                  salgs[i].dims,ins1,signs);
  
         if pos = [] then 
            return false;
         fi;
        #now find isomorphism from calg in db to my
         cand   := corelg.carrierAlgDB[pos[1]];
         dbenum := cand.enum;
         dbcf   := cand.cfImage;   #coefficient of sigma(xi) wrt yi   
         csl2   := cand.cfCsl2*Sqroot(1);    #cayley sl2 triple
        
        #adjust ordering wrt enum
         perm := MappingPermListList(myenum,dbenum);
         l1   := [1..Length(mys1)];
         l2   := [Length(mys1)+1..l];
         if not (ForAll(l1,x->x^perm in l1) and ForAll(l2,x->x^perm in l2)) then  
            Error("cannot use this db entry!");
         fi;
         mycg2  := [[],[],[]];
         mycf2  := ShallowCopy(mycf);
         for j in [1..l] do
            mycg2[1][j^perm] := mycg[1][j];
            mycg2[2][j^perm] := mycg[2][j];
            mycg2[3][j^perm] := mycg[3][j];
            mycf2[j^perm]     := mycf[j];
         od;
         myenum := List(myenum,x->x^perm);
         mycg   := mycg2;
         mycf   := mycf2;
         cf     := List([1..l],x->Sqroot(dbcf[x]/mycf[x]));
         mycgSF := [[],[],[]];
         for j in [1..l] do
             mycgSF[1][j] := cf[j]*writeToSF(mycg[1][j]);
             mycgSF[2][j] := cf[j]^-1*writeToSF(mycg[2][j]);
             mycgSF[3][j] := writeToSF(mycg[3][j]);
         od;
         canbas := corelg.myflat(corelg.mySLAfctCanBas(LSF, mycgSF));
         csl2   := List(csl2,x->x*canbas);
         if i = 1 then
            newtr := csl2;
         else
            for j in [1..3] do newtr[j] := newtr[j]+csl2[j]; od;
         fi;
   
     #else we have a principal carrier algebra
      else

         if salgs[i].type = ["A",1] then
            x    := Basis(mygr[2][1])[1];
            h    := salgs[i].h;
            cf   := Coefficients(Basis(VectorSpace(CF(4),[h],"basis"),
                               [h]),x*sigma(x))[1];
            x    := Sqroot(1/cf)*writeToSF(x);
            csl2 := [sigmaSF(x),writeToSF(h),x];
         else
            tr := rec(g0 := BasisVectors(Basis(mygr[1][1])),
                      gp := List(mygr[2],x->BasisVectors(Basis(x))),
                      gn := List(mygr[3],x->BasisVectors(Basis(x))));

            my   := SubalgebraNC(L,corelg.myflat(Concatenation(tr.g0,tr.gp,tr.gn)),"basis");
            cs   := corelg.ChevalleySystemInnerType(my,tr,sigma).chevSys;
            r    := Concatenation(cs[1],cs[2]);
           #r    := Filtered([1..Length(r)], x-> corelg.eltInSubspace(L,tr.gp[1],r[x])); 
            r    := Filtered([1..Length(r)], x-> r[x] in tr.gp[1] );
            if r = [] then Error("ups.. torus as s0; case A1?"); fi; 
            mat  := List(Concatenation(cs[3]{r},[-salgs[i].h]),
                     x->Coefficients(Basis(L),x));
            ns   := NullspaceMat(mat)[1];
            cf   := List(ns{[1..Length(r)]}, x->Sqroot(x));
            x    := Sum(List([1..Length(r)], i-> cf[i]*(writeToSF(Flat(cs)[r[i]])))); 
            csl2 :=  [sigmaSF(x),writeToSF(salgs[i].h),x]; 
         fi;
         if i = 1 then
            newtr := csl2;
         else
            for j in [1..3] do newtr[j] := newtr[j]+csl2[j]; od;
         fi;
      fi;
   od;
   out.realsl2   := corelg.CayleyTransform(newtr,SqrtField);
   Unbind(out.carrier);
   Unbind(out.oldsl2);
   corelg.checkTriples(form,[out]);
   return true;
end;


################################################
#Input:  real form
#Output: record with entries
#           - form: the input form
#           - triples: record with entries:
#               -principal:    cdims, realsl2
#               -nonprincipal: cdims, realsl2
#           - tobedone: those indices i where
#                 nonprincipal[i] has NO realsl2;
#                 here we need to solve an equation!
#################################################
corelg.RealCayleyTriplesOfRealForm := function(form)
local triples, nonpr, lookup, tr, ctr, wdd;
   Info(InfoCorelg,4,"start RealSL2Triples");
   triples   := corelg.principalOrbitsOfRealForm(form);
   nonpr     := triples.nonprincipal;
   Info(InfoCorelg,1,"there are ",Length(nonpr)," non-principal triples");
   lookup := List([1..Length(nonpr)],x->corelg.lookupRealCayleyTriple(form,nonpr[x]));
   lookup := Filtered([1..Length(nonpr)],x->lookup[x]=false);
   if Length(lookup)>0 then
      Print("TO BE SOLVED (getEquationsForX): ",Length(lookup)," carrier algs,\n");
      Print("their dims are",List(nonpr{lookup},x->Sum(Flat(x.cdims))),"\n");
   elif Length(nonpr)>0 then
      Print("could solve all non-principal carrier algebras\n");
   fi;
   Info(InfoCorelg,4,"end RealSL2Triples");
   return rec(form     := form,
              triples  := triples, 
              tobedone := lookup);
end;


################################################
#Input:  output of RealCayleyTriplesOfRealForm with tobedone=[]
#Output: corresponding nilpotent orbits
################################################
corelg.ConvertRealCayleyTriplesToNilpotentOrbits := function( res )
local T, wdd, cb, L, LSF, form, new, i, o;

   form := res.form;
   if not res.tobedone = [] then
      Error("tobedone is not emtpy");
   fi;
   new := [];
   L   := form.liealg;
   LSF := form.liealgSF;
   T   := SignatureTable(L);
   cb  := BasisNC(LSF,corelg.myflat(ChevalleyBasis(LSF)));
   Info(InfoCorelg,4,"compute WDDs and coefficients for constructing orbits");
   for i in res.triples.principal do
      wdd := corelg.WDD(L, List(Coefficients(Basis(LSF),i.realsl2[2]),
                         SqrtFieldEltToCyclotomic)*Basis(L),T);
      o   := NilpotentOrbit( L, wdd );
      SetRealCayleyTriple(o, i.realsl2);
      SetInvariants(o, rec(wdd             := wdd,
                           carrierAlgebra  := rec(dims   := i.cdims,
                                                  principal := true)));
        
     #SetCoefficientsWRTChevBasis( o, List(i.realsl2,x->Coefficients(cb,x)));
      Add(new,o);
   od;
   for i in res.triples.nonprincipal do
      wdd := corelg.WDD(L, List(Coefficients(Basis(LSF),i.realsl2[2]),
                         SqrtFieldEltToCyclotomic)*Basis(L),T);
      o   := NilpotentOrbit( L, wdd );
      SetRealCayleyTriple(o, i.realsl2);
      SetInvariants(o, rec(wdd             := wdd,
                           carrierAlgebra  := rec(dims   := i.cdims,
                                                  principal := true)));        
     #SetCoefficientsWRTChevBasis( o, List(i.realsl2,x->Coefficients(cb,x)));
      Add(new,o);
   od;
   Info(InfoCorelg,4,"done");
   return rec(form := form, nilpotentOrbits := new);
end;

#################################################################################
#################################################################################
#
# FUNCTIONS TO GET, VIEW AND ATTACH SOLUTIONS FOR NON-PRINCIPAL CAlgs
#
#################################################################################
#################################################################################


################################################
#Input:  this is only called from TryToFindComplexCayleyTriple
#        res is output of RealCayleyTriplesOfRealForm
#        res.triples.nonprincipal[i] had attached eqs, GR etc
#Output: tries to attach a solution of the equations GR;
#        otherwise prints GR so that equations can maybe
#        solved manually
################################################
corelg.viewReducedEquationsAndAttach := function(res,i)
local out, eqs, ok, notok, isgood,j,jj ,isgood2, ok2, notok2, 
      goodvars, k, var, str;
  str := "";
  out := res.triples.nonprincipal[i];
  if not IsBound(out.eqs) then
     Error("no eqs attached");
  fi;
  isgood2 := function(l)
     return Length(l)=4 and Length(l[1])=2 and l[1][2] = 1 
            and l[2]=1 and IsGaussRat(l[4]) and
            ForAll(List([1..Length(l[3])/2],j->l[3][2*j-1]),x-> x in goodvars);
  end;
  isgood := function(l)
     return Length(l)=4 and l[1]=[] and IsGaussRat(l[2])
            and Length(l[3])=2 and l[3][2]=2 and l[4]=1;
  end;
  if not out.GR = fail then
     Print("these are equations:\n",out.GR,"\n");
     eqs   := out.GR;
     ok    := [];
     notok := [];
     for j in eqs do
        if isgood(ExtRepPolynomialRatFun(j)) then
           Add(ok,ExtRepPolynomialRatFun(j));
        else
           Add(notok,j);
        fi;
     od;
     if not ok=[] then
        goodvars := List(ok,x->x[3][1]);
        Print("\ncan copy-paste the following part (modify var name 'res'):\n");
        Append(str,Concatenation( "corelg.attachSolution(res,",String(i),",["));
        Print("corelg.attachSolution(res,",String(i),",[");
        for jj in [1..Length(ok)] do
           j := ok[jj];
           Append(str,Concatenation("[",String(j[3][1]),",Sqroot(",String(-j[2]),")]"));
           Print("[",String(j[3][1]),",Sqroot(",String(-j[2]),")]");
           if not jj = Length(ok) or not notok=[] then Append(str,","); Print(","); fi;
        od;
        if notok=[] then
           Append(str,"]);");
           Print("]);\n\n");
           Print("now try to attach this solution:\n");
           EvalString(str);
        else
           ok2    := [];
           notok2 := [];
           for j in notok do
              if isgood2(ExtRepPolynomialRatFun(j)) then
                 Add(ok2,ExtRepPolynomialRatFun(j));
              else
                 Add(notok2,j);
              fi;
           od;
           for jj in [1..Length(ok2)] do
              j := ok2[jj];
              Append(str, Concatenation("[",String(j[1][1]),",",String(-j[4]),"*Sqroot(" ));
              Print("[",j[1][1],",",-j[4],"*Sqroot(");
              for k in [1..Length(j[3])/2] do
                 var := j[3][2*k-1];
                 var := Filtered(ok,x->x[3][1] = var)[1];
                 Append(str,Concatenation("(",String(-var[2]),")"));
                 Print("(",-var[2],")");
                 if not k = Length(j[3])/2 then 
                    Print("*"); 
                    Append(str, "*");
                 fi;
              od;
              Append(str,")]");
              Print(")]");
              if not jj = Length(ok2) or not notok2 = [] then 
                 Append(str,","); 
                 Print(","); 
              else
                
              fi;
           od;          
           if notok2=[] then
              Append(str,"]);");
              Print("]);\n\n");
              Print("now try to attach this solution:\n");
              EvalString(str);
           else
              Print(",\n\n");
              Print("still to take into account:\n",notok2,"\n\n"); 
           fi;
        fi;
     fi;
  else
     Print("no reduced equations attached\n");
  fi;
end;





################################################
#Input:  this is only called from  corelg.viewReducedEquationsAndAttach    
#Output: attach a solution (complex Cayley triple)
#        to res.triples.nonprincipal[i] and write it to corelg.carrierAlgDB
################################################
corelg.attachSolution := function(res,i,v)
local csl2,j,w,tmp, out, bas, realsl2, L, sigmaSF, calg,  enum, form,
      rs, cg, gr, l, s1, s1c, cf, CB, cfcsl2, data, pos, perm, F, LSF,
      cgSF, path;

   if Length(corelg.carrierAlgDB)=0 then corelg.readDBCA(); fi;

   form := res.form;
   out  := res.triples.nonprincipal[i];
  
   if not IsBound(out.makeSol) then 
      Error("cannot attach solution, no makeSol entry");
   fi;
  
  #make sure everything is over SqrtField
   for i in [1..Length(v)] do v[i][2] := v[i][2]*Sqroot(1); od;

  #make solution
   csl2          := out.makeSol(v);
   realsl2       := corelg.CayleyTransform(csl2,SqrtField); 
   tmp           := rec(realsl2 := realsl2);
   corelg.checkTriples(res.form,[tmp]);
   if IsBound(out.realsl2) then 
      Error("Warning: realsl2 was alread bound!");
   fi;
   if not Sum(Flat(out.cdims)) = Dimension(res.form.liealg) then
      Error("dim of calg smaller than dim of alg; no need to store");
   fi;
   out.realsl2   := realsl2; 
   out.cayleysl2 := csl2;
   corelg.checkTriples(res.form,[out]); 

  #now rearrange can gens and store everything wrt can bas
   L       := res.form.liealg;
   LSF     := res.form.liealgSF;
   sigmaSF := RealStructure(LSF);
   calg := SubalgebraNC(L,corelg.myflat(Concatenation(out.carrier.g0, 
                                            out.carrier.gp,out.carrier.gn)),"basis");
   SetCartanSubalgebra(calg,Intersection(calg,CartanSubalgebra(L)));
   rs   := RootSystemOfZGradedLieAlgebra(calg,out.carrier);
   SetRootSystem(calg,rs);
   cg   := CanonicalGenerators(RootSystem(calg)); 
   gr   := List([[out.carrier.g0], out.carrier.gp, out.carrier.gn],
                 x-> List(x,y->SubspaceNC(L,y,"basis")));
   l    := Length(cg[1]);
  #s1   := Filtered([1..l], x -> corelg.eltInSubspace(calg,Basis(gr[2][1]),cg[1][x]));
   s1   := Filtered([1..l], x -> cg[1][x] in gr[2][1] );
   s1c  := Filtered([1..l], x -> not x in s1);
   cg   := List(cg, x-> Concatenation(x{s1},x{s1c}));
   cgSF := List(cg,x->List(x,res.form.writeToSF));
   cf   := List([1..l],x->
                Coefficients(Basis(VectorSpace(SqrtField,[cgSF[2][x]],"basis"),
                             [cgSF[2][x]]),sigmaSF(cgSF[1][x]))[1]);

   CB     := List(corelg.myflat(SLAfcts.canbas( calg, cg )),res.form.writeToSF);
   CB     := Basis(VectorSpace(SqrtField,CB,"basis"),CB);
   cfcsl2 := List(csl2,x->Coefficients(CB,x));

  #adjust enumeration according to rearranged cangens
   enum   := CartanType(CartanMatrix(rs)).enumeration[1];
   perm   := PermList(Concatenation(s1,s1c))^-1;
   enum   := List(enum,x->x^perm);
   tmp    := rec(type := res.form.type, 
              rank := res.form.rank,
              enum := enum,
              dims := List(gr,x->List(x,Dimension)),
              ins1 := AsSortedList(List([1..Length(s1)],x->Position(enum,x))),
              cfImage  := cf,          #sigma(cg[1][i]) = cf[i]*cg[2][i] 
              ordsigns := List(cf{enum},corelg.SqrtEltMySign),
              cfCsl2   := cfcsl2);      #cf of csl2 wrt can bas wrt cg
  #ins1 are the indicies of std ordered roots which lie in s1
  #ordsigns are the signs of cf, in order of std ordered roots
   
   Info(InfoCorelg,4,"   read database before writing");
   ReadPackage( "corelg", "gap/carrierAlg.db" );
   pos := Filtered(corelg.carrierAlgDB,x->  x.rank = tmp.rank and
                                x.type = tmp.type and 
                                x.dims = tmp.dims and
                                x.ins1 = tmp.ins1);
 
   if pos=[] then
      Add(corelg.carrierAlgDB, tmp);
    ##path := Concatenation(LOADED_PACKAGES.corelg[1],"/gap/carrierAlg.db");
      path := Filename(DirectoriesPackageLibrary("corelg","gap"),"carrierAlg.db");

      PrintTo(path,"corelg.carrierAlgDB:=");
      AppendTo(path,corelg.carrierAlgDB);
      AppendTo(path,";");
      Info(InfoCorelg,4,"  wrote new entry to database");
   else
      Info(InfoCorelg,4,"  entry was already contained in database");
   fi;
   Unbind(out.GR);
   Unbind(out.eqs);
   Unbind(out.var);
   Unbind(out.makeSol);
   Info(InfoCorelg,4,"update to-do list; check other unfinished nonprinciple triples");
   for i in res.triples.nonprincipal{res.tobedone} do 
      corelg.lookupRealCayleyTriple(form,i); 
   od; 
   res.tobedone := Filtered([1..Length(res.triples.nonprincipal)],
                   x->not IsBound(res.triples.nonprincipal[x].realsl2));
   Info(InfoCorelg,4,"this is new to-do list",res.tobedone);
   return true;
end;



################################################
#Input:  res is output of RealCayleyTriplesOfRealForm 
#        ind should lie in res.tobedone
#        nrvars list of nr of non-zero variables for
#               Groebner basis
#        nrtries list of tries for each number in nrvars
#                or just one integer (nr of tries)
#        IMPORTANT: name of variable 'res' MUST BE "res"
#                   in order to save solution automatically
#Output: computes equations to construct complex
#        Cayley triple and (hopefully) attaches it
#        automatically;
#################################################
corelg.TryToFindComplexCayleyTriple := function(res, ind, nrvars,nrtries)
local L, sigma, ca, s1, sl2, h, bas, s1b, new, n, PR, prb, lhs, rhs,
      eqs, rateqs, makeSol, complexConjugate, eq, i, I, GR, tmp, j,
      eqs1, done, ii, k, nzvars, form, out, sigmaSF, LSF, s1bSF, hSF; 

   if ind>Length(res.triples.nonprincipal) then
      Error("there are only ",Length(res.triples.nonprincipal)," nonprincipals\n");
   fi;
   form := res.form;
   out  := res.triples.nonprincipal[ind];
   if IsInt(nrtries) then 
      nrtries := ListWithIdenticalEntries(Length(nrvars),nrtries); 
   fi;

   L   := form.liealg;
   LSF := form.liealgSF;
   if IsBound(out.realsl2) then
      Print("real sl2 triple already attached, don't do anything.\n");
      return true;
   fi;
   if IsBound(out.eqs) then
      Print("equations already attached; re-compute them.\n");
   fi;
 
   sigma   := RealStructure(L);
   sigmaSF := RealStructure(LSF);

   ca    := SubalgebraNC(L,corelg.myflat(Concatenation(out.carrier.g0,
                                    out.carrier.gp,out.carrier.gn)),"basis"); 
   s1    := SubspaceNC(ca,out.carrier.gp[1],"basis");
   sl2   := out.oldsl2;
   h     := sl2[2];
   bas   := Basis(ca);
   s1b   := Basis(s1);
   s1bSF := List(s1b,form.writeToSF);    
   hSF   := form.writeToSF(h);

  #set up equations
  #new[i][j] = s1b(i) * sigma(s1b(i));
   new   := List(s1b, g-> List(s1b,h-> Coefficients(bas,g*sigma(h))));

  #if X=sum a_i s1b(i) then the system of equation is
  # (a_1,..,a_n) * new * sigma(a_1,..,a_n)^T = h
  #write ai = ui + E(4)vi; create indeterminates:
   n   := Length(s1b);
   PR  := PolynomialRing(CF(4),[1..2*n]);
   prb := IndeterminatesOfPolynomialRing(PR);
   lhs := List([1..n],i->prb[i]+prb[i+n]*(E(4)));
   rhs := List([1..n],i->prb[i]-prb[i+n]*(E(4)));
   new := lhs*new;
   new := Sum(List([1..n],i-> new[i]*rhs[i]));
   new := new - Coefficients(bas,h);
   eqs := Filtered(new,x->not x = 0*x);
   
  #write complex equation over rationals
   complexConjugate := function(eq)
   local tmp, fam, i;
      fam := FamilyObj(eq);
      tmp := MutableCopyMat(ExtRepPolynomialRatFun(eq));
      for i in [1..Length(tmp)] do
         tmp[i] := ComplexConjugate(tmp[i]);
      od;
      return PolynomialByExtRep(fam,tmp);  
   end;
   rateqs := [];
   for eq in eqs do
       tmp := complexConjugate(eq);
       Add(rateqs,(1/2*E(4))*(eq-tmp));
       Add(rateqs,(1/2)*(eq+tmp));
   od;
   eqs := Filtered(rateqs,x-> not x = x*0);
  
   Info(InfoCorelg,4,"there are ",Length(prb)," variables.");
   for ii in [1..Length(nrvars)] do
      i := nrvars[ii];
      j := nrtries[ii];
      Info(InfoCorelg,4,"start checking ",i," nonzero variables");
      for k in [1..j] do
         Print("start try nr ",k,"\n");
         nzvars := [];
         while Length(nzvars)<i do
            tmp := Random(prb);
            if not tmp in nzvars then Add(nzvars,tmp); fi;
         od;        
         eqs1 := Concatenation(eqs,Filtered(prb,x-> not x in nzvars)); 
         I    := Ideal(PR,eqs1);
         GR   := HasTrivialGroebnerBasis(I);
         done := not GR;
         if done then break; fi;
      od; 
      if done then break; fi;
   od;
   if done then 
      Print("found nontriv GR, now compute it\n");
      GR := ReducedGroebnerBasis(I,MonomialLexOrdering());
      GR := Filtered(GR,x-> not x in prb);
   else
      GR := fail;
      Print("Error:could not find nice Groebner basis\n");
   fi;
  
  #given solution vector l, creates solution triple
   makeSol := function(v)
   local e,f,j,l;
      l := ListWithIdenticalEntries(2*n,0);
      for j in v do l[j[1]] := j[2]; od;
      e := List([1..n],i-> l[i]+E(4)*l[i+n]);
      e := Sum(List([1..n],i->e[i]*s1bSF[i]));
      f := sigmaSF(e);
      return [f,hSF,e];
   end;

   out.eqs     := eqs;
   out.GR      := GR;
   out.var     := prb;
   out.makeSol := makeSol; 
   if not GR = fail then
      corelg.viewReducedEquationsAndAttach(res,ind);
   fi;
   return;
end;


#################################################
#Input:  ()
#Output: shows entries in database corelg_carrierAlgDB
#################################################
corelg.calgDBentries := function()
local i, tmp;
   if Length(corelg.carrierAlgDB)=0 then corelg.readDBCA(); fi;
   for i in ["A","B","C","D","G","E","F"] do 
      tmp := Collected(List(Filtered(corelg.carrierAlgDB,x->x.type=i),y->y.rank)); 
      Print("for ",i," have ",tmp,"\n"); 
   od;
end;



#################################################################################
#################################################################################
#
# THE FUNCTIONS FOR WRITING / READING / THE DATABASE realTriples.de
# (all real forms of simple complex LAs up to rank 10)
#
#################################################################################
#################################################################################


#################################################
#Input:  type and rank
#Output: adds all nilpotent orbits of given LA to
#        database realTriples.db
#################################################
corelg.WriteRealNilpotentOrbitsToDB := function(type, rank)
local form, res, triples, f, new, tr, T, newtr, L, LSF, tmp, cf, i, cnt, cb, path;

   if Length(corelg.realtriplesDB)=0 then corelg.readDBTriples(); fi;

   form := corelg.NonCompactRealFormsOfSimpleLieAlgebra(type,rank);
   for f in form do
      if ForAny(corelg.realtriplesDB, x-> 
                x.form = RealFormParameters(f.liealg)) then
         Info(InfoCorelg,4,"form ",RealFormParameters(f.liealg)," already in DB");
      else
         L    := f.liealg; 
         LSF  := f.liealgSF;
         cb   := BasisNC(LSF,corelg.myflat(ChevalleyBasis(LSF)));
         new  := rec( form := RealFormParameters(f.liealg), triples :=[]);
         tr   := corelg.RealCayleyTriplesOfRealForm(f);
         if not tr.tobedone=[] then
            Error("sth wrong, there are entries in tobedone!");
         fi;
         T   := SignatureTable(L);
         cnt := 1;
         for i in tr.triples.principal do
            Info(InfoCorelg,4," principal triple ",cnt," of ",Length(tr.triples.principal));
            cnt         := cnt + 1; 
            newtr       := rec();
            cf          := List(i.realsl2,x->Coefficients(Basis(LSF),x));
            tmp         := List(cf,x->Filtered([1..Length(x)],j->not x[j]=Zero(SqrtField)));
            newtr.rct   := List([1..3],x->Flat(List(tmp[x],y->[y,cf[x][y]])));
            cf          := List(i.realsl2,x->Coefficients(cb,x));
            tmp         := List(cf,x->Filtered([1..Length(x)],j->not x[j]=Zero(SqrtField)));
            newtr.rctcb := List([1..3],x->Flat(List(tmp[x],y->[y,cf[x][y]])));
            newtr.cdims := i.cdims;
            newtr.princ := true;
            newtr.wdd   := corelg.WDD(L, 
                               List(Coefficients(Basis(LSF),i.realsl2[2]),
                                    SqrtFieldEltToCyclotomic)*Basis(L),       
                               T);
            Add(new.triples,newtr);
         od;
         cnt := 1;
         for i in tr.triples.nonprincipal do
            Info(InfoCorelg,4," nonprincipal triple ",cnt," of ",Length(tr.triples.nonprincipal));
            cnt         := cnt + 1; 
            newtr       := rec();
            cf          := List(i.realsl2,x->Coefficients(Basis(LSF),x));
            tmp         := List(cf,x->Filtered([1..Length(x)],j->not x[j]=Zero(SqrtField)));
            newtr.rct   := List([1..3],x->Flat(List(tmp[x],y->[y,cf[x][y]])));
            cf          := List(i.realsl2,x->Coefficients(cb,x));
            tmp         := List(cf,x->Filtered([1..Length(x)],j->not x[j]=Zero(SqrtField)));
            newtr.rctcb := List([1..3],x->Flat(List(tmp[x],y->[y,cf[x][y]])));
            newtr.cdims := i.cdims;
            newtr.princ := false;
            newtr.wdd   := corelg.WDD(L, 
                               List(Coefficients(Basis(LSF),i.realsl2[2]),
                                    SqrtFieldEltToCyclotomic)*Basis(L),       
                               T);
            Add(new.triples,newtr);
         od;
         ReadPackage( "corelg", "gap/realTriples.db" ); 
         Info(InfoCorelg,4,"  re-read corelg.realtriplesDB");
         Add(corelg.realtriplesDB,new);  
       ##path := Concatenation(LOADED_PACKAGES.corelg[1],"/gap/realTriples.db");  
         path := Filename(DirectoriesPackageLibrary("corelg","gap"),"realTriples.db"); 
         PrintTo(path,"corelg.realtriplesDB:=");
         AppendTo(path,corelg.realtriplesDB);
         AppendTo(path,";");
         Info(InfoCorelg,4,"  wrote new entry to corelg.realtriplesDB");
      fi;
   od;
   return true;
end;






##############################################################################
##
##  displays the parameters of the real forms (of type <type>) for which the 
##  database contains contains its real nilpotent orbits
##
corelg.RealNilpotentOrbitsInDatabase := function(arg)
local i, tmp, t;
   if Length(corelg.realtriplesDB)=0 then corelg.readDBTriples(); fi;
   if Length(arg) = 1 then 
      t := [arg[1]];
   else 
      t := ["A","B","C","D","G","E","F"];
   fi;
   for i in t do
      tmp := Filtered(corelg.realtriplesDB,x->x.form[1]=i);
      Print("Triples of type ",i,"\n");
      for i in List(tmp,x->x.form) do Display(i); od;
   od;
end;



##############################################################################
##  returns all real nilpotent orbits in real form L, reads orbits reps from database:
##  at the moment, the database contains A2-A8, B2-B10, C2-C10, D4-D8, F4, G2, E6-E8
##  
corelg.RealNilpotentOrbitsFromDatabase := function(LL)
local type, rank, pars,param, i,  res, kacs, L, LSF, form, new, dim,orb, neworb,
      K, P, ff, ee, hh, tmp, sigma, theta, n, k, makeVec, o, forms, iso, cd,cg, H,
      h,e,f,cf, db;

   Info(InfoCorelg,1,"start RealNilpotentOrbitsFromDatabase");
   if Length(corelg.realtriplesDB)=0 then
      corelg.readDBTriples();
   fi;

   if not LeftActingDomain(LL)=SqrtField then
     Error("need LA over SqrtField");
   fi;

   tmp   := VoganDiagram(LL);
   param := tmp!.param;
   if Length(param)>1 then 
      Error("nilpotent orbits only for simple LAs");
   fi;
   param := [param[1][1],param[1][2]];
   Add(param,ShallowCopy(Signs(tmp)));
   Add(param,PermInvolution(tmp));

  #compact form?
   if IdRealForm(LL)[3] = 1 then
      return [];
   fi;

  # deal with case A1?
  if param[1]="A" and param[2]=1 then
     new :=[];
     H   := MaximallyNonCompactCartanSubalgebra(LL);
     cd  := CartanDecomposition(LL);
     cg  := CanonicalGenerators(RootsystemOfCartanSubalgebra(LL,H));
     h   := cg[3][1];
     e   := cg[1][1];
     f   := cg[2][1];
     cf  := Coefficients(Basis(SubspaceNC(LL,[h],"basis"),[h]),-e*cd.CartanInv(e))[1];
     e   := (1/Sqrt(cf))*e;
     f   := -cd.CartanInv(e); 
     if not e*f=h or not h*e=2*e or not h*f=-2*f or not cd.CartanInv(e)=-f then
        Error("wrong real triple");
     fi;  
     o := NilpotentOrbit( LL, [2] );
     SetRealCayleyTriple(o, [e,h,f]);
     SetInvariants(o, rec(wdd := [2]));
     Add(new,o);
     e := -e;
     f := -cd.CartanInv(e);
     if not e*f=h or not h*e=2*e or not h*f=-2*f or not cd.CartanInv(e)=-f then
        Error("wrong real triple");
     fi;  
     o := NilpotentOrbit( LL, [2] );
     SetRealCayleyTriple(o, [e,h,f]);
     SetInvariants(o, rec(wdd := [2]));
     Add(new,o);
     return new;
  fi;



   new := rec();
   db  := First(corelg.realtriplesDB,x-> x.form = param);
   if db = fail then
      Error("cannot find entry with these parameters",param);
   fi;

   if IsBound(LL!.std) and LL!.std then
      Info(InfoCorelg,2,"  don't have to construct isomorphism...");
      L := LL;
      n := Dimension(L);
      iso := IdentityMapping(L);
  else
      L := RealFormById(param[1],param[2],IdRealForm(LL)[3],SqrtField); 
      n := Dimension(L);
      VoganDiagram(L);
      Info(InfoCorelg,2,"  construct isomorphism...");
      iso    := IsomorphismOfRealSemisimpleLieAlgebras(L,LL);
      Info(InfoCorelg,2,"  ...done");
   fi;

   new.form   := rec(type:=param[1], rank:=param[2], liealgSF := LL);
   Info(InfoCorelg,2,"  read triples...");     

  #writes compressed coef vector to coef vector
   makeVec := function(v)
   local vec,i;
      vec := ListWithIdenticalEntries(n,Zero(SqrtField));
      for i in [1..Length(v)/2] do
         vec[v[2*i-1]] := v[2*i]*One(SqrtField);
      od;
      return vec;
   end;
    
   new.nilpotentOrbits := List(db.triples, x->
                       rec(rct   := List(x.rct,i->Image(iso,makeVec(i)*Basis(L))),
                           rctcb := List(x.rctcb,makeVec),
                           cdims := x.cdims,
                           princ := x.princ,
                           wdd   := x.wdd));

   for i in [1..Length(new.nilpotentOrbits)] do
      tmp := new.nilpotentOrbits[i];
      o   := NilpotentOrbit( LL, tmp.wdd );
      SetRealCayleyTriple(o, tmp.rct);
      SetInvariants(o, rec(wdd             := tmp.wdd,
                           carrierAlgebra  := rec(dims   := tmp.cdims,
                                               principal := tmp.princ)));
     #SetCoefficientsWRTChevBasis( o, tmp.rctcb);
      new.nilpotentOrbits[i] := o;
   od;
   sigma := RealStructure(LL);
   theta := CartanDecomposition(LL).CartanInv;
 
   Info(InfoCorelg,2,"  all triples constructed, now test them");
   for tmp in new.nilpotentOrbits do
      ff := RealCayleyTriple(tmp)[1];
      hh := RealCayleyTriple(tmp)[2];
      ee := RealCayleyTriple(tmp)[3];
      if not (hh*ff = -2*ff and hh*ee=2*ee and ee*ff=hh) or
         not theta(ee)=-ff or 
         not  RealCayleyTriple(tmp) = List(RealCayleyTriple(tmp),sigma) then
            Error("not a real sl2 triple");
      fi;
   Info(InfoCorelg,2,"  all triples OK");
   od;
   Info(InfoCorelg,1,"end RealNilpotentOrbitsFromDatabase");

   return new.nilpotentOrbits;
   
end;


#################################################
InstallMethod( NilpotentOrbitsOfRealForm,
   "for Lie algebras",
   true,
   [ IsLieAlgebra ], 0, function(L)

   if HasNilpotentOrbitsOfRealForm(L) then 
      return NilpotentOrbitsOfRealForm(L); 
   fi;
   
   return corelg.RealNilpotentOrbitsFromDatabase(L);
end);




##############################################################################
InstallGlobalFunction(CarrierAlgebraOfNilpotentOrbit, function(L, orb)
local sl2, i, j, g0, h, ca, tmp, esp, hm, K, z, t, grad, gr, tt, zz, old;
 
   sl2 := RealCayleyTriple(orb);
   h   := sl2[2];
   z   := LieCentraliser(L,SubalgebraNC(L,sl2,"basis"));
   zz  := LieDerivedSubalgebra(z);
   if Dimension(zz)>0 then
      t   := MaximallyNonCompactCartanSubalgebra(zz);
   else
      t := zz;
   fi;
  
   tmp := Concatenation(Basis(t),Basis(LieCentre(z)));
   if tmp =[] then t:=SubalgebraNC(L,[],"basis"); else  t   := SubalgebraNC(L,tmp);fi;

   z   := LieCentraliser(L,t);
   hm  := TransposedMat( AdjointMatrix(Basis(z),h) );
   esp := [[], [],[]];
   i   := 0;
   repeat
      old := Length(corelg.myflat(esp));
      if i=0 then 
         esp[1][1] := List(NullspaceMat(hm),x->x*Basis(z)); 
      else
         esp[2][i/2] := List(NullspaceMat(hm-i*hm^0),x->x*Basis(z));
         esp[3][i/2] := List(NullspaceMat(hm+i*hm^0),x->x*Basis(z));
      fi;
      i := i+2;
   until Length(corelg.myflat(esp))=old;

   K          := LieDerivedSubalgebra( SubalgebraNC( L, corelg.myflat(esp),"basis"));
   grad       := [[],[],[]];
   grad[1][1] := BasisVectors( Basis( Intersection(K,SubspaceNC(L,esp[1][1],"basis")) ) );
   for i in [1..Length(esp[2])] do
      grad[2][i] := BasisVectors( Basis( Intersection(K,SubspaceNC(L,esp[2][i],"basis")) ) );
   od;
   for i in [1..Length(esp[3])] do
      grad[3][i] := BasisVectors( Basis( Intersection(K,SubspaceNC(L,esp[3][i],"basis")) ) );
   od;

   return rec(liealg := K, grading := grad, wdd := Invariants(orb).wdd, dims := List(esp,x->List(x,Length)));
end);



#################################################
corelg.sortOrbitsByWDD := function(t,r,nr)
local L, orbs, wdds, res, i, j, cdims, tot,ok;

   L    := RealFormById(t,r,nr);
   orbs := corelg.RealNilpotentOrbitsFromDatabase(L);
   orbs := List(orbs,x->CarrierAlgebraOfNilpotentOrbit(L,x));
   wdds := List(Collected(List(orbs,x->x.wdd)),x->x[1]);
   res  := [];
   for i in [1..Length(wdds)] do
      res[i] := rec(wdd:=wdds[i], 
      orbs   := List(Filtered(orbs,x->x.wdd=wdds[i]),
                     y->rec(ca:=y.dims, orb :=y)));
   od; 
   Print("display only those wdds which at least two orbs with the same calg\n\n");
   ok  := 0;
   tot := 0;
   for i in res do
      if Length(i.orbs)>1 then
         tot := tot+1;
         cdims := List(i.orbs,x->x.ca);
         if not IsDuplicateFreeList(cdims) then
            Print("orbits with wdd ",i.wdd,":\n");
            for j in i.orbs do
               Print("   calg with dim ",j.ca,"\n");
            od;
            Print("\n");
         else
            ok := ok+1;
         fi;
      fi;
   od;
   Print("there are ",tot, " WDDs which have more than one orbits attached\n");
   Print(ok," of these have pairwise distinct carrier alg dims\n");
   if tot-ok=0 then
      Print("Hence ALL orbits can be distinguished by their wdd and calg dims\n");
   else
      Print("Hence: ",tot-ok," wdds have orbits which cannot be distinguised by their calg dims\n");
   fi;

   return rec(alg := L, res:=res, data:=[t,r,nr,tot,tot-ok]);
end;


#################################################
corelg.sortOrbitsByDim := function(t,r,nr)
local L, orbs, dims, res, i, ok, tot, wdds,j; 

   L    := RealFormById(t,r,nr);
   orbs := corelg.RealNilpotentOrbitsFromDatabase(L);
   orbs := List(orbs,x->CarrierAlgebraOfNilpotentOrbit(L,x));
 
   dims := List(Collected(List(orbs,x->x.dims)),x->x[1]);
   res  := [];
   for i in [1..Length(dims)] do
      res[i] := rec(dim  :=dims[i], 
                    orbs := List(Filtered(orbs,x->x.dims=dims[i]), y->rec(wdd:=y.wdd, orb :=y)));
   od;

   ok  := 0;
   tot := 0;
   for i in res do
      wdds := Collected(List(i.orbs,x->x.wdd));
      Print("dim ",i.dim, " has the following wdds\n");
      for i in wdds do Print("    ",i[1],"\n"); od;
   od;

end;
