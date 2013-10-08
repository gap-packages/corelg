############################################################################
#
#  this file contains a few preliminary functions
#
#

#################################################
#Input:  nilpotent matrix
#Output: its Jordan NF and base change
#################################################

## jb := function (n)
## local m, i;
##    m := 0*IdentityMat(n);
##    for i in [1..n-1] do m[i][i+1] := 1; od;
##    return m;
## end;
## mydbm := function ( mats )
##     local  n, M, m, d;
##     n := Sum( mats, function ( m )
##             return Length( m[1] );
##         end );
##     M := NullMat( n, n );
##     n := 0;
##     for m  in mats  do
##         d := Length( m );
##         M{[ 1 .. d ] + n}{[ 1 .. d ] + n} := m;
##         n := n + d;
##     od;
##     return M;
## end;

## rjb := function(n)
## local part, gr, k,i;
##    part := Random(Partitions(n));
##    Print("took these part ",part,"\n");
##    gr   := mydbm(List(part,x->jb(x)));
##    repeat k := RandomMat(n,n); until not Determinant(k)=0;
##    return k*gr*k^-1;
## end;

## could use that for the two adjoint matrices of nilpotent elts
## (that is, two nilp matrices in sln). If the JNF has only even 
## blocks and the determinant of base change matrices are positive
## in one case, negative in the other, then the two elts are not
## conjugate
corelg.JNFofNilpotent2 := function(m)
local F, V, kern, rws, i, new, elts, exps, preim, toExp, bas, mm;

  Print("new version\n");
  F    := DefaultFieldOfMatrix(m);
  V    := VectorSpace(F,m^0);
  kern := Subspace(V,NullspaceMat(m));
  rws := [V];
  i   := 0;
  repeat 
     i   := i+1;
     new := Subspace(V,m^i);
     if not Dimension(new)=0 then Add(rws,new); fi;
  until Dimension(new)=0;
  rws := Reversed(rws);
  elts := List(BasisVectors(Basis(rws[1])),x->[x]);;
  exps := ListWithIdenticalEntries(Length(elts),0);
  if not ForAll(elts,x->x[1]*m=0*x[1]) then Error("1"); fi;
  for i in [2..Length(rws)] do
     mm    := List(BasisVectors(Basis(rws[i])),x->x*m);
     preim := List(elts,x->SolutionMat(mm,x[Length(x)])); #need preim in rws[i]
     preim := List(preim,x-> x*BasisVectors(Basis(rws[i])));
     toExp := List(elts,x->x[1]);
     new   := BasisVectors(Basis(Intersection(kern,rws[i])));
     new   := List(BaseSteinitzVectors(new,toExp).factorspace,x->[x]);
     for i in [1..Length(elts)] do Add(elts[i],preim[i]); od;
     elts  := Concatenation(elts,new);
     exps  := Concatenation(exps+1,ListWithIdenticalEntries(Length(new),0));
     if not ForAll(elts,x->x[1]*m=0*x[1]) then Error("1"); fi;
  od;
  SortParallel(exps,elts);
  bas := Concatenation(elts); 
  return [bas*m*bas^-1,Determinant(bas)];
 #return [bas,bas*m*bas^-1,exps];
end;
	

