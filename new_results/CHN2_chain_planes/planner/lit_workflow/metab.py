"""Is the (e0,e1) 6x6 representation, restricted to W=U(e0,e1)e_1, metabelian ([L',L']=0)?  exact (sympy rationals)."""
import sys; sys.dont_write_bytecode=True
sys.path.insert(0,'.')
import sympy as sp
import phi
from fractions import Fraction as F
def span_close(vecs, ops, maxit=20):
    """closure of span(vecs) under linear ops (list of sympy Matrices)."""
    B=sp.Matrix.hstack(*vecs)
    for _ in range(maxit):
        new=[op*B for op in ops]
        B2=sp.Matrix.hstack(B,*new)
        # column space basis
        cs=B2.columnspace()
        B2=sp.Matrix.hstack(*cs) if cs else sp.zeros(B.rows,1)
        if B2.shape[1]==B.shape[1] and B2.rank()==B.rank():
            # closure reached if span unchanged
            if sp.Matrix.hstack(B,B2).rank()==B.rank(): return B2
        B=B2
    return B
def analyse(a,label):
    k=phi.kin_from_a(a); e0,e1=phi.e_mats(**k)
    E0=sp.Matrix([[sp.Rational(x.numerator,x.denominator) if isinstance(x,F) else sp.Rational(x) for x in r] for r in e0])
    E1=sp.Matrix([[sp.Rational(x.numerator,x.denominator) if isinstance(x,F) else sp.Rational(x) for x in r] for r in e1])
    v0=sp.Matrix([1,0,0,0,0,0])
    W=span_close([v0],[E0,E1])
    d=W.shape[1]
    # restrict operators to W: solve W*R = E*W
    def restrict(E):
        return (W.T*W).inv()*W.T*(E*W)
    R0=restrict(E0); R1=restrict(E1)
    assert (W*R0-E0*W).is_zero_matrix and (W*R1-E1*W).is_zero_matrix
    def com(X,Y): return X*Y-Y*X
    # Lie algebra generated: span of nested commutators
    def vec(M): return list(M)
    def basis_of(mats):
        if not mats: return []
        A=sp.Matrix([vec(m) for m in mats])
        # row space basis
        rs=A.rowspace()
        return [sp.Matrix(d,d,list(r)) for r in rs]
    g=basis_of([R0,R1])
    for _ in range(12):
        new=g+[com(x,y) for x in g for y in g]
        g2=basis_of(new)
        if len(g2)==len(g): break
        g=g2
    # derived algebra L' = [g,g]
    Lp=basis_of([com(x,y) for x in g for y in g])
    Lpp=basis_of([com(x,y) for x in Lp for y in Lp])
    print('%-28s a=%-14s dimW=%d dim g=%d dim L\'=%d dim [L\',L\']=%d'%(label,a,d,len(g),len(Lp),len(Lpp)))
    return d,len(g),len(Lp),len(Lpp)
if __name__=='__main__':
    for lab,pts in (('Pi1 (a1=a2,a3=a4)',[(1,1,2,2),(2,2,1,1),(1,1,3,3),(2,2,5,5)]),
                    ('Pi2 (a1=a3,a2=a4)',[(1,1,2,2),(1,2,1,2),(2,1,2,1),(3,5,3,5)]),
                    ('generic positive',[(1,1,1,2),(1,2,2,1),(3,1,2,5)]),
                    ('diagonal (1,1,1,1)',[(1,1,1,1)]),
                    ('Pi1 intersect Pi2',[(1,1,1,1),(2,2,2,2)])):
        for p in pts:
            if lab.startswith('Pi1 (') and p[0]!=p[1]: continue
            analyse(p,lab)
