#########################################################################
#
# this file contains the methods for defining the field SqrtField
# and for working with its elements, see sqrt.gd
#
#

###############################################
# Do deal with Flat for lists of matrices
#
corelg.myflat := function(L)
   return Concatenation(L);
end;


###############################################
# DUE TO A BUG IN BasisNC
InstallMethod( NiceBasisNC, 
               [IsBasisByNiceBasis and HasNiceBasis], NiceBasis );


#########################################################################

BindGlobal( "SqrtField_objectify", function( elm )
local u,isrt;
   isrt:= Length(elm) = 1 and Length(elm[1][2])=0;
   u:= Objectify( SqrtFieldType, [ Immutable(elm) ] );
   u![2]:= isrt;
   return u;
end ); 

#########################################################################

InstallMethod( \in,
   "for an object and a SqrtField",
   [ IsObject, IsSqrtField ], 1000000,
function( p, qf )
   return IsSqrtFieldElement(p);
end );

#########################################################################


InstallGlobalFunction(Sqroot, function(q)
local d, n, fc, cf, ps, i, m, sgn;

   if not IsRat(q) then 
      if IsSqrtFieldElement(q) and Length(q![1]) = 1 
         and IsRat(q![1][1][1]) and q![1][1][2] = [] then
         return Sqroot(q![1][1][1]);
      else
         Error("<q> must be a rational, or rational SqrtField elt"); 
      fi;
   fi;
   if q < 0 then sgn := E(4); q := -q; else sgn := 1; fi;
   if q in [0,1] then 
      return SqrtField_objectify([[sgn*q,[]]]); 
   fi;
   d  := DenominatorRat(q);
   n  := NumeratorRat(q);
   fc := FactorsInt(n);
   Append(fc,FactorsInt(d));
   fc := Collected(fc);
   cf := 1/d;
   ps := [];
   for i in [1..Length(fc)] do
      m  := fc[i][2] mod 2;
      cf := cf *(fc[i][1]^((fc[i][2]-m)/2));
      if m=1 then Add(ps,fc[i][1]); fi;
   od;
   if ps = [] then return 
      SqrtField_objectify([[sgn*cf,[]]]); 
   fi; 
   if ps[1] = 1 then ps := ps{[2..Length(ps)]}; fi;
   Sort(ps);
   MakeImmutable(ps);
   return SqrtField_objectify([[sgn*cf,ps]] );
end); 

#########################################################################

SetOne( SqrtField, Sqroot(1) );
SetZero( SqrtField, Sqroot(0) );
SetZero( SqrtFieldFam, Zero(SqrtField) );
SetOne( SqrtFieldFam, One(SqrtField) );
SetIsUFDFamily(SqrtFieldFam,true);

#########################################################################

InstallMethod( ViewObj, 
        "for SqrtField elements",
        [ IsSqrtFieldElement ],

function( sf )
local elt, l, i, printmon, printcf;
   elt := sf![1];
   l   := Length(elt);
   printmon := function(mon)
      if mon[1] = 1 then 
         Print( "Sqroot(",Product(mon[2]),")");
      elif mon[1]<0 or not (IsRat(mon[1]) or IsRat(E(4)*mon[1])) 
                    or (IsRat(E(4)*mon[1]) and E(4)*mon[1]>0) then
         Print( "(",mon[1],")*","Sqroot(",Product(mon[2]),")");
      else 
         Print( mon[1],"*","Sqroot(",Product(mon[2]),")");
      fi;
   end;
   if l=1 then
      if elt[1][2]=[] then Print(elt[1][1]); else printmon(elt[1]); fi;
   else
      if elt[1][2]=[] then 
         Print(elt[1][1]," + "); 
      else 
         printmon(elt[1]); 
         Print(" + "); 
      fi;
   fi;
   for i in [2..l-1] do printmon(elt[i]); Print(" + "); od;
   if l>1 then printmon(elt[l]); fi;  
end );

##########################################################################

InstallMethod( PrintObj, 
        "for SqrtField elements",
        [ IsSqrtFieldElement ],
function( sf )
local elt, l, i, printmon;
   elt := sf![1];
   l   := Length(elt);
   printmon := function(mon)
      if mon[1] = 1 then 
         Print( "Sqroot(",Product(mon[2]),")");
      elif mon[1]<0  or not (IsRat(mon[1]) or IsRat(E(4)*mon[1]))  
                     or (IsRat(E(4)*mon[1]) and E(4)*mon[1]>0) then
         Print( "(",mon[1],")*","Sqroot(",Product(mon[2]),")");
      else 
         Print( mon[1],"*","Sqroot(",Product(mon[2]),")");
      fi;
   end;
   if l=1 then
      if elt[1][2]=[] then Print(elt[1][1]); else printmon(elt[1]); fi;
   else
      if elt[1][2]=[] then 
         Print(elt[1][1]," + "); 
      else 
         printmon(elt[1]); 
         Print(" + "); 
      fi;
   fi;
   for i in [2..l-1] do printmon(elt[i]); Print(" + "); od;
   if l>1 then printmon(elt[l]); fi;

end );

#########################################################################


InstallGlobalFunction( SqrtFieldIsGaussRat, function(el)
   return IsSqrtFieldElement(el) and el![2];
end); 

#########################################################################


InstallGlobalFunction( SqrtFieldMakeRational, function(el)
   if IsSqrtFieldElement(el) then
      if SqrtFieldIsGaussRat(el) then
         return el![1][1][1]; 
      else 
         return false; 
      fi;
   elif IsSqrtFieldElementCollection(el) then
      if ForAll(el,SqrtFieldIsGaussRat) then 
         return List(el,x->x![1][1][1]);
      else
         return false;
      fi;
   elif IsSqrtFieldElementCollColl(el) then
      if ForAll(Flat(el),SqrtFieldIsGaussRat) then
         return List(el,x->List(x,y->y![1][1][1]));
      else
         return false;
      fi; 
   fi;
   return false;
end); 

#########################################################################

InstallMethod( Determinant,
   "for SqrtField matrices",
   [ IsSqrtFieldElementCollColl], 1000,
function(x)
local a,b;
   a := SqrtFieldMakeRational(x);
   if not a = false then 
      return DeterminantMat(a)*One(SqrtField);
   fi;
   TryNextMethod();
end);

#########################################################################

InstallMethod( \*,
   "for SqrtField matrices",
   [ IsSqrtFieldElementCollColl, IsSqrtFieldElementCollColl],
function(x,y)
local a,b;
   a := SqrtFieldMakeRational(x);
   if not a = false then 
      b := SqrtFieldMakeRational(y);      
      if not b = false then
         return One(SqrtField)*(a*b);
      fi;
   fi;
  #Print("DONT USE MAKERAT\n");
  TryNextMethod();
end);


#########################################################################

InstallMethod( \+,
   "for SqrtField matrices",
   [ IsSqrtFieldElementCollColl, IsSqrtFieldElementCollColl],
function(x,y)
local a,b;
   a := SqrtFieldMakeRational(x);
   if not a = false then 
      b := SqrtFieldMakeRational(y);      
      if not b = false then
         return One(SqrtField)*(a+b);
      fi;
   fi;
  #Print("DONT USE MAKERAT\n");
   TryNextMethod();
end);


#########################################################################

InstallMethod( \*,
   "for rational and SqrtField elements",
   [IsCyclotomic, IsSqrtFieldElement], 
function(x,y)
local b,i;
      if not IsGaussRat(x) then
         Error("cyclotomic has to be Gaussian Rational");
      fi;
      if x = 0 then return Zero(SqrtField); fi;
      if x = 1 then return y; fi;
      b := List(y![1],ShallowCopy);
      for i in b do i[1] := x*i[1]; od;
      return SqrtField_objectify( b );          
end);

#########################################################################

InstallMethod( \*,
   "for SqrtField and rational elements",
   [IsSqrtFieldElement, IsCyclotomic], 
function(x,y)
   return y*x;      
end);

#########################################################################

InstallMethod( ZeroOp,
   "for SqrtField element",
   [IsSqrtFieldElement], 
function(x)
   return Zero(SqrtField);     
end);

#########################################################################

InstallMethod( OneOp,
   "for SqrtField element",
   [IsSqrtFieldElement], 
function(x)
   return One(SqrtField);   
end);

#########################################################################

InstallMethod( \+,
   "for SqrtField and SqrtField elements",
   [IsSqrtFieldElement, IsSqrtFieldElement], 10000,
function(x,y)
local a,b,i,j,e,len,pr;

   if x = Zero(SqrtField) then return y; fi;
   if y = Zero(SqrtField) then return x; fi;
   if x![2] then
      if y![2] then
         return SqrtField_objectify([[x![1][1][1]+y![1][1][1], []]]);
      else 
         return x![1][1][1] + y;
      fi;
   fi; 
   if y![2] then return y![1][1][1] + x; fi;


   if Length(x![1])<=Length(y![1]) then
      a := List(x![1],ShallowCopy);
      b := List(y![1],ShallowCopy);
   else
      a := List(y![1],ShallowCopy);
      b := List(x![1],ShallowCopy);
   fi;
   i := 1;
   for j in [1..Length(a)] do
      e  := a[j];
      pr := Product(e[2]);
      while i<= Length(b) and Product(b[i][2])<pr do i := i+1; od;
      if i>Length(b) then 
         Append(b,a{[j..Length(a)]});
         return SqrtField_objectify( b ); 
      fi;
      if e[2]=b[i][2] then
         b[i][1] := b[i][1] + e[1];
         if b[i][1] = 0 then
            #remove position i
             len := Length(b);
             CopyListEntries( b, i + 1, 1, b, i, 1, len - i );
             Unbind( b[len] );
         fi;
      else
        #add e at position i
         CopyListEntries(b,i,1,b,i+1,1,Length(b)-i+1);
         b[i] := e;
      fi;
   od;
   if b = [] then return Zero(SqrtField); fi;
   return SqrtField_objectify( b );          
end);

#########################################################################

InstallMethod( \+,
   "for Gaussian rational and SqrtField elements",
   [IsCyclotomic, IsSqrtFieldElement], 
function(x,y)
local t;
   if not IsGaussRat(x) then
      Error("cyclotomic has to be Gaussian rational");
   fi;
   if x=0 then return y; fi;
   t := List(y![1],ShallowCopy);
   if t[1][2] = [] then
      t[1][1] := t[1][1] + x;
      if t[1][1] = 0 then t := t{[2..Length(t)]}; fi;
      if t = [] then return Zero(SqrtField); fi;
   else
     #add [y,[]] at position 1
      CopyListEntries( t, 1, 1, t, 2, 1, Length(t) );
      t[1] := [x,[]];
   fi;
   return SqrtField_objectify( t );
end);

#########################################################################

InstallMethod( \+,
   "for SqrtField and rational elements",
   [IsSqrtFieldElement, IsCyclotomic], 
function(x,y)
   return y+x;
end);

#########################################################################

InstallMethod( AdditiveInverseSameMutability,
   "for SqrtField elements",
   [IsSqrtFieldElement], 
function(x)
local a, i;
   a := List(x![1],ShallowCopy);
   for i in a do i[1] := -i[1]; od;
   return SqrtField_objectify( a );   
end);

#########################################################################

InstallMethod( \=,
   "for SqrtField elements",
   [IsSqrtFieldElement, IsSqrtFieldElement], 
function(x,y)
   return x![1] = y![1];
end);
#########################################################################

InstallMethod( \=,
   "for SqrtField and rational elements",
   [IsSqrtFieldElement, IsRat], 
function(x,y)
   if not x![2] then return false; fi;
   return SqrtFieldMakeRational(x) = y;
end);

#########################################################################

InstallMethod( \=,
   "for rational and SqrtField element",
   [IsRat, IsSqrtFieldElement], 
function(x,y)
   if not y![2] then return false; fi;
   return SqrtFieldMakeRational(y) = x;
end);

#########################################################################


InstallMethod( \*,
   "for SqrtField and SqrtField elements",
   [IsSqrtFieldElement, IsSqrtFieldElement],
   function(x,y)
   local a,b,i,j,res,int,cf,v,w, mns, cfs, isSmaller, pos, k;

      if x=Zero(SqrtField) or y=Zero(SqrtField) then 
         return Zero(SqrtField); 
      fi;
      if x=One(SqrtField) then return y; fi;
      if y=One(SqrtField) then return x; fi;

      a   := x![1];
      b   := y![1];

      if x![2] and y![2] then
         return SqrtField_objectify([[a[1][1]*b[1][1],[]]]);
      fi; 

      mns := [ ];
      cfs := [ ];

      isSmaller := function(u,v) return Product(u)< Product(v); end;

      for i in a do
         for j in b do

            int := i[2]{[1..Length(i[2])]};
            INTER_SET( int, j[2] );

            cf  := i[1]*j[1];  #*Product(int);
            for k in int do cf:= cf*k; od;

            res := i[2]{[1..Length(i[2])]};
            APPEND_LIST( res, j[2] );
            v   := Filtered(res,x->not x in int);
            SORT_LIST(v);
            pos := POSITION_SORTED_LIST_COMP( mns, v, isSmaller );
            if IsBound(mns[pos] ) and mns[pos] = v then
               cfs[pos]:= cfs[pos]+cf;
            else
               Add( mns, v, pos );
               Add( cfs, cf, pos );
            fi;

         od;
      od;
      res:= [ ];
      for i in [1..Length(mns)] do
         if cfs[i] <> 0 then
            Add( res, [ cfs[i], mns[i] ] );
         fi;
      od;
      if res = [] then return Zero(SqrtField); fi;
      return SqrtField_objectify( res );
end);


##########################################################################

InstallMethod( InverseOp,
   "for SqrtField elements",
   [IsSqrtFieldElement], 
function(el)
local b, i, j, bas, mat, elts, a, cfs, l,cf, pos, ns, inv, res, v, nrp, want;
   if el=Zero(SqrtField) then Error("elt must be nonzero"); fi;
   a    := List(el![1],ShallowCopy);
   a    := SqrtField_objectify( a );
   i    := 1;
   elts := [One(SqrtField)];
   bas  := [[]];
   cfs  := [ [1] ];
   l    := Length(bas);
   nrp  := Length(Collected(Flat(List(a![1],x->x[2]))));
   want := 2^Minimum([nrp,Length(Filtered(a![1],x->not x[2]=[]))]);
   if want > 16 then
      Info(InfoSqrtField,2,"need ",want," powers for computing the inverse");
   fi;
   repeat
      i       := i+1;
      elts[i] := elts[i-1]*a; 
      cfs[i]  := ListWithIdenticalEntries(l,0);
      for j in elts[i]![1] do
         pos := Position(bas,j[2]);
         if pos=fail then 
            Add(bas,j[2]); 
            l   := l+1;  
            pos := l; 
            for v in cfs do v[l] := 0; od;
         fi;
         cfs[i][pos] := j[1];
      od;
      if i mod 10 = 0 then
         Info(InfoSqrtField,2,"  InfoSqrtField: (inverses) computed ",i," powers");
      fi;
   until i-1 = want;
   ns  := NullspaceMat(cfs);
   res := ns[1];
   inv := Sum(List([2..Length(res)],i->res[i]*elts[i-1]))/(-res[1]);
  #Info(InfoSqrtField,3,"test inverse");
  #if not el*inv = One(SqrtField) then Error("hmpf... inv wrong!"); fi;
   return inv;       
end);


#########################################################################

InstallMethod( AdditiveInverseMutable,
   "for SqrtField elements",
   [IsSqrtFieldElement],
   function(x)
   local a, i;
      a := List(x![1],ShallowCopy);
      for i in a do i[1] := -i[1]; od;
      return SqrtField_objectify( a );
end);

#########################################################################

InstallMethod( \<,
   "for SqrtField elements",
   [IsSqrtFieldElement, IsSqrtFieldElement],
   function(x,y)
      return x![1] < y![1];
end);

#########################################################################

InstallMethod( \<,
   "for SqrtField element and Rational",
   [IsSqrtFieldElement, IsRat],
   function(x,y)
      if not x![2] then TryNextMethod(); fi;
      return SqrtFieldMakeRational(x) < y;
end);

########################################################################

InstallMethod( \<,
   "for rational and SqrtField element",
   [IsRat, IsSqrtFieldElement],
   function(x,y)
      if not y![2] then TryNextMethod(); fi;
      return SqrtFieldMakeRational(y) > x;
end);


#########################################################################


InstallGlobalFunction(SqrtFieldMinimalPolynomial, function(el)
#returns polynomial f of min deg over GaussianRats such that f(el)=0; 
local b, i, j, bas, mat, elts, a, cfs, l,cf, pos, ns, inv, res, v, nrp, 
      want, pol, x;
   if not IsSqrtFieldElement(el) then
      Error("input has to be SqrtField elt");
   fi;
   if el=Zero(SqrtField) then Error("elt must be nonzero"); fi;
   a    := List(el![1],ShallowCopy);
   a    := SqrtField_objectify( a );
   i    := 1;
   elts := [One(SqrtField)];
   bas  := [[]];
   cfs  := [ [1] ];
   l    := Length(bas);
   nrp  := Length(Collected(Flat(List(a![1],x->x[2]))));
   want := 2^Minimum([nrp,Length(Filtered(a![1],x->not x[2]=[]))]);
   if want > 16 then
      Info(InfoSqrtField,1,"need ",want," powers for computing the inverse");
   fi;
   repeat
      i       := i+1;
      elts[i] := elts[i-1]*a; 
      cfs[i]  := ListWithIdenticalEntries(l,0);
      for j in elts[i]![1] do
         pos := Position(bas,j[2]);
         if pos=fail then 
            Add(bas,j[2]); 
            l   := l+1;  
            pos := l; 
            for v in cfs do v[l] := 0; od;
         fi;
         cfs[i][pos] := j[1];
      od;
      if i mod 10 = 0 then
         Info(InfoSqrtField,2,"  - computed ",i," elements");
      fi;
   until i-1 = want;
   ns   := NullspaceMat(cfs);
   res  := ns[1]*One(SqrtField);
   x    := Indeterminate(SqrtField);
   elts := [One(SqrtField)];
   Append(elts, List([2..Length(res)],i->x^(i-1)));
   pol  := Sum(List([1..Length(res)],i->res[i]*elts[i]));
   return pol;       
end);

#########################################################################


InstallGlobalFunction(SqrtFieldEltByRationalSqrt, function(e)
local sgn;
   if not IsRat(e^2) then Error("input has to be Sqroot of a rational"); fi;
   if not e-Sqrt(e^2) = 0 then sgn := -1; else sgn := 1; fi;
   return sgn*Sqroot(e^2);
end); 

##########################################################################


InstallGlobalFunction(SqrtFieldEltRealAndComplexPart, function(e)
local real, imag, i, makecf, cf, r;
   if not IsSqrtFieldElement(e) then 
      Error("input must be SqrtFieldElement"); 
   fi;
   makecf := function(cf)
      local cj;
      cj := ComplexConjugate(cf);
      return [1/2*(cf+cj),1/2*(cf-cj)];
   end;
   real := Zero(SqrtField);
   imag := Zero(SqrtField);
   for i in e![1] do
      r  := Sqroot(Product(i[2]));
      cf := makecf(i[1]);
      real := real+cf[1]*r;
      imag := imag+cf[2]*r;
   od;
   return [real,imag];
end); 
##########################################################################


InstallGlobalFunction(IsPosSqrtFieldElt, function(e)
local real, imag, i, makecf, cf, r;
   if not IsSqrtFieldElement(e) or not Length(e![1]) = 1 or
      not ForAll(e![1],i->IsRat(i[1])) then 
      Error("input must be a real SqrtFieldElement with one monom"); 
   fi;
   return e![1][1][1]>0;
end); 

##########################################################################


InstallGlobalFunction(CoefficientsOfSqrtFieldElt, function(e)
local cfs, i, x;
   if not IsSqrtFieldElement(e) then 
      Error("input must be SqrtFieldElement"); 
   fi;
   x := List(e![1],ShallowCopy);
   return List(x,i->[i[1],Product(i[2])]);
end); 

##########################################################################


InstallGlobalFunction(SqrtFieldEltByCoefficients, function(e)
local elt, i;
   elt := Zero(SqrtField);
   for i in e do
      if not IsGaussRat(i[1]) and IsRat(i[2]) then 
          Error("wrong format of coefficients vector"); 
      fi;
      elt := elt + i[1]*Sqroot(i[2]);
   od;
   return elt;
end); 

##########################################################################


InstallGlobalFunction(SqrtFieldEltToCyclotomic, function(e)
local makemon;
   if not IsSqrtFieldElement(e) then
      Error("input has to lie in SqrtField");
   fi;
   if Length(e![1])=1 and e![1][1][2]=[] then return e![1][1][1]; fi;
  #Info(InfoSqrtField,1,"Warning: might not work for large cyclotomics.");
   makemon := function(mon)
      if mon[2]=[] then return mon[1]; fi;
      return mon[1]*Sqrt(Product(mon[2]));
   end;
   return Sum(List(e![1],makemon));
end); 

##########################################################################

InstallMethod( DefaultFieldOfMatrix,
    "method for a matrix over SqrtField",
    [ IsMatrix and IsSqrtFieldElementCollColl ],
function( mat )
   return SqrtField;
end );

##########################################################################

InstallMethod( DefaultFieldByGenerators,
    "method for a list over SqrtField",
    [ IsList and IsSqrtFieldElementCollection ],
function( mat )
   return SqrtField;
end );

##########################################################################

InstallMethod( ComplexConjugate,
    "method for SqrtField elt",
    [ IsSqrtFieldElement ],
function( el )
local x, i;
   x := List(el![1], ShallowCopy);
   for i in x do i[1] := ComplexConjugate(i[1]); od;
   return SqrtField_objectify( x );
end );

##########################################################################

InstallMethod( Sqrt,
    "method for SqrtField elt",
    [ IsSqrtFieldElement ],
function( el )
local x, i;
   x := List(el![1],ShallowCopy);
   if not (Length(x) = 1 and x[1][2] = [] and IsRat(x[1][1])) then
      Error("input has to be rational SqrtFieldElt");
   fi;
   return  Sqroot(x[1][1]);
end );

##########################################################################

InstallMethodWithRandomSource( Random,
    "for a random source and an SqrtField ",
    [ IsRandomSource, IsSqrtField ],
function( rs, F )
   return Sum(List([1..15],x-> Random(rs,Rationals)*Sqroot(Random(rs,Rationals))));
   
end );

##########################################################################

InstallMethod( AbsoluteValue,
    "for SqrtField element",
    [ IsSqrtFieldElement ], 10000,
function( el )
local x,i;
   x := List(el![1],ShallowCopy);
   if not (Length(x)=1 and IsRat(x[1][1])) then
      Error("elt must be a real SqrtFieldElt with one summand");
   fi;
   for i in x do i[1] := AbsoluteValue(i[1]); od;
   return SqrtField_objectify( x );
end );

##########################################################################

InstallMethod( String,
    "for SqrtField element",
    [ IsSqrtFieldElement ],
function( el )
local elt, l, i, printmon, printcf, str;
   elt := List(el![1],ShallowCopy);
   str := "";
   l   := Length(elt);
   printmon := function(mon)
      if mon[1] = 1 then 
         Append(str, "Sqroot("); 
         Append(str,String(Product(mon[2])));
         Append(str,")");
      elif mon[1]<0 or not (IsRat(mon[1]) or IsRat(E(4)*mon[1])) 
                    or (IsRat(E(4)*mon[1]) and E(4)*mon[1]>0) then
         Append(str,"(");
         Append(str,String(mon[1]));
         Append(str,")*Sqroot(");
         Append(str,String(Product(mon[2])));
         Append(str,")");
      else 
         Append(str, String(mon[1]));
         Append(str,"*Sqroot(");
         Append(str,String(Product(mon[2])));
         Append(str,")");
      fi;
   end;
   if l=1 then
      if elt[1][2]=[] then Append(str,String(elt[1][1])); else printmon(elt[1]); fi;
   else
      if elt[1][2]=[] then 
         Append(str,String(elt[1][1]));
         Append(str," + "); 
      else 
         printmon(elt[1]); 
         Append(str," + ");
      fi;
   fi;
   for i in [2..l-1] do printmon(elt[i]); Append(str," + "); od;
   if l>1 then printmon(elt[l]); fi;  
   return str;
end );

##########################################################################

InstallMethod( String,
    "for SqrtField",
    [ IsSqrtField ],
function( el )
   return "SqrtField";
end );

##########################################################################


InstallGlobalFunction(SqrtFieldEltByCyclotomic, function(arg)
local e, cj, re, im, solvereal, res, r1, r2, 
      sqfree, whichSqfreeRootsInCF;
   e := arg[1];
   if e in SqrtField then return e; fi;
   if not IsCyclotomic(e) then 
      Error("input must be Cyclotomic"); 
   fi;
   if IsGaussRat(e) then return e*One(SqrtField); fi;

   whichSqfreeRootsInCF := function(n)
  #returns also one for technical reasons
   local sqfree, eeven, res;
      sqfree := x -> ForAll(Collected(FactorsInt(x)),i->i[2]=1);
      eeven  := x -> Length( Filtered(List(Collected(FactorsInt(x)), i->i[1]),
                             p -> p mod 4 = 3)) mod 2 = 0;
      res    := Filtered(DivisorsInt(n),x->sqfree(x));
      if IsOddInt(n) or (IsInt(n/2) and IsOddInt(n/2)) then
         return Filtered(res,x->IsOddInt(x) and eeven(x));
      fi;
      if n mod 4 = 0 and not n mod 8 =0 then
         return Filtered(res,IsOddInt);
      fi;
      return res;
   end;

  #check if input is sqrt of rational
  #test is expensive, e.g. for e=Sqrt(1231)+Sqrt(17);;
  #if IsRat(e^2) then return SqrtFieldEltByRationalSqrt(e); fi;

  #split into real and complex part
   cj    := ComplexConjugate(e);;
   re    := 1/2*(e+cj);;
   im    := 1/2*E(4)^3*(e-cj);;

  #the main function which requires a real input   
   solvereal := function( r )
   local n, pr, prsq, bas, mat, ns;  
      if IsRat(r) then return r*One(SqrtField); fi;
      n     := Conductor(r);
      pr    := whichSqfreeRootsInCF(n);
      prsq  := List(pr,Sqrt);
      bas   := Basis(CF(n));
      mat   := List(prsq,x->Coefficients(bas,x));
      ns    := SolutionMat(mat,Coefficients(bas,r));
      if ns = fail then return fail; fi;
      return Sum(List([1..Length(pr)],i->ns[i]*Sqroot(pr[i])));
   end;

  #try to solve real and complex part
   r1  := solvereal(re); if r1 = fail then return fail; fi;
   r2  := solvereal(im); if r2 = fail then return fail; fi;
   res := r1 + E(4)*r2;
 ###TEST
  #Info(InfoSqrtField,1,"test res");
  #if not e = SqrtFieldEltToCyclotomic(res) then
  #   Error("res wrong");
  #fi;
   return res;
end); 


#########################################################################


InstallGlobalFunction( SqrtFieldPolynomialToRationalPolynomial, function(e)
local x,i,j,cf;
   x  := Indeterminate(Rationals);
   cf := List(CoefficientsOfUnivariatePolynomial(e),SqrtFieldEltToCyclotomic);
   return Sum(List([1..Length(cf)],i-> cf[i]*x^(i-1)));
end); 


#########################################################################


InstallGlobalFunction( SqrtFieldRationalPolynomialToSqrtFieldPolynomial, function(e)
local x,i,j,cf;
   x  := Indeterminate(SqrtField);
   cf := List(CoefficientsOfUnivariatePolynomial(e),SqrtFieldEltByCyclotomic);
   return Sum(List([1..Length(cf)],i-> cf[i]*x^(i-1)));
end); 



##########################################################################

InstallOtherMethod( Factors,
    "for SqrtField",
    [ IsPolynomialRing, IsUnivariatePolynomial ],
function( R, f )
local rf;
if LeftActingDomain(R) = SqrtField and 
      IsSqrtFieldElementCollection(CoefficientsOfUnivariatePolynomial(f)) then
      rf := SqrtFieldPolynomialToRationalPolynomial(f);
      rf := Factors(rf);
      return List(rf,SqrtFieldRationalPolynomialToSqrtFieldPolynomial);
   else
      TryNextMethod();
   fi;
end );


