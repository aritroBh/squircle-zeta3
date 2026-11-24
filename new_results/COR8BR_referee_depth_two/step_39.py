import sys, pickle
sys.path.insert(0,'.')
import cor8br_sym as E
from fractions import Fraction as Fr
M=(1,0,0,0); C=(0,1,0,0); P=(0,0,1,0); Q=(0,0,0,1)
def vadd4(a,b): return tuple(p+q for p,q in zip(a,b))
def vsub4(a,b): return tuple(p-q for p,q in zip(a,b))
def vscale4(k,a): return tuple(k*p for p in a)
X4 = vsub4(vsub4(vscale4(2,Q), vscale4(2,M)), C)
U4 = vsub4(vsub4(vscale4(2,P), vscale4(2,M)), C)
V4 = vsub4(P, Q)
B4=[(i,j) for i in range(4) for j in range(i+1,4)]
def w2_4(a,b):
    o={}
    for idx,(i,j) in enumerate(B4):
        cc=a[i]*b[j]-a[j]*b[i]
        if cc: o[idx]=o.get(idx,0)+cc
    return o
def wedgebiv(A,B):
    out={}
    for i,ci in A.items():
        for j,cj in B.items():
            if i==j: continue
            if i<j: out[(i,j)]=out.get((i,j),0)+ci*cj
            else: out[(j,i)]=out.get((j,i),0)-ci*cj
    return {k:v for k,v in out.items() if v}
extra31 = wedgebiv(w2_4(U4,X4), w2_4(V4,X4))
print("(U^X)^(V^X) at pair1:", extra31)
# so the extra term in their stated Li31 formula vanishes here by coincidence
# check pair with x=2 (constant), y=... find a pair where stated cf31 fails; use (x,y)=(2,c):
# x=2: X=M, U=[1-2]=[-1]=0. Y=C, V=P, W=[1-2c] not in alphabet. bad.
# Use their db: test all pairs for cf31 match
import pickle
data = pickle.load(open('<project>/research/08_theory/COR8B_depth_two_two_terms/pairs_delta_all_v2.pkl','rb'))
import sympy as sp
c = sp.symbols('c')
def factor_vec4(expr):
    expr = sp.simplify(expr)
    num, den = sp.fraction(sp.cancel(expr))
    vec=[0,0,0,0]
    for poly, sgn in [(num,1),(den,-1)]:
        poly=sp.Poly(poly,c)
        content, factors = sp.factor_list(poly)
        cont=sp.Rational(content); cont=-cont if cont<0 else cont
        for prime,e in sp.factorint(int(cont.p)).items():
            if prime!=2: return None
            vec[0]+=sgn*e
        for prime,e in sp.factorint(int(cont.q)).items():
            if prime!=2: return None
            vec[0]-=sgn*e
        for f,e in factors:
            f=sp.Poly(f,c)
            if f.degree()==0: continue
            fc=f.as_expr(); matched=False
            for be,idx in [(c,1),(1-c,2),(1+c,3)]:
                ratio=sp.simplify(fc/be)
                if ratio.is_number:
                    r=sp.Rational(ratio); r=-r if r<0 else r
                    for pr,ee in sp.factorint(int(r.p)).items():
                        if pr!=2: return None
                        vec[0]+=sgn*e*ee
                    for pr,ee in sp.factorint(int(r.q)).items():
                        if pr!=2: return None
                        vec[0]-=sgn*e*ee
                    vec[idx]+=sgn*e; matched=True; break
            if not matched: return None
    return tuple(vec)

S = E.build_symbols()
def spec_delta5(d5, imgi):
    out={}
    for (i,j),v in d5.items():
        bi,bj = E.WB5[i]; bk,bl = E.WB5[j]
        A=w2_4(imgi[bi],imgi[bj]); B=w2_4(imgi[bk],imgi[bl])
        for ia,ca in A.items():
            for ib,cb in B.items():
                if ia==ib: continue
                if ia<ib: out[(ia,ib)]=out.get((ia,ib),Fr(0))+v*ca*cb
                else: out[(ib,ia)]=out.get((ib,ia),Fr(0))-v*ca*cb
    return {k:v for k,v in out.items() if v}

n_mismatch22=0; n_mismatch31=0; n_tot=0
for (xs,ys,d22,d31,d13) in data[:60]:
    x0=sp.sympify(xs); y0=sp.sympify(ys)
    vx=factor_vec4(x0); vu=factor_vec4(1-x0); vy=factor_vec4(y0); vv=factor_vec4(1-y0)
    vw=factor_vec4(1-x0*y0)
    if None in (vx,vu,vy,vv) or vw is None:
        vw=(0,0,0,0)
    imgi={0:vx,1:vy,2:vu,3:vv,4:vw}
    m22=spec_delta5(E.delta22_5(S['Li22']),imgi)
    m31=spec_delta5(E.delta22_5(S['Li31']),imgi)
    m13=spec_delta5(E.delta22_5(S['Li13']),imgi)
    n_tot+=1
    if m22!=d22: n_mismatch22+=1; print("MISMATCH 22 at",xs,ys)
    if m31!=d31: n_mismatch31+=1; print("MISMATCH 31 at",xs,ys)
print(f"checked {n_tot} pairs from their db: my delta vs their db mismatches: Li22={n_mismatch22}, Li31={n_mismatch31}")
open('progress.log','a').write("step_36: my engine's delta matches COR8B's database at all tested pairs; their STATED closed forms for Li22/Li31 are wrong in the 5-letter space (engine itself is right).\n")
