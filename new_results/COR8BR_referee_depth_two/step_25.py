engine = r'''
# cor8br_sym.py -- COR8BR's own exact symbol engine (written from definitions)
# Two-variable letters: formal Q-vector space with basis
#   X=[x], Y=[y], U=[1-x], V=[1-y], W=[1-xy];  Z:=[xy]=X+Y.
# Symbols: dict {tuple of 5-vectors: Fraction}.
# delta_22(a b c d) = (a^b)^(c^d) with the bivector wedge ANTISYMMETRIC
#   (paper convention: vanishes on products and Li_4).
from fractions import Fraction as Fr

X=(1,0,0,0,0); Y=(0,1,0,0,0); U=(0,0,1,0,0); V=(0,0,0,1,0); W=(0,0,0,0,1)
Z=(1,1,0,0,0)
def vadd(a,b): return tuple(p+q for p,q in zip(a,b))
def vsub(a,b): return tuple(p-q for p,q in zip(a,b))
def vneg(a): return tuple(-p for p in a)
ZERO5=(0,0,0,0,0)

def sadd(s1,s2,s2scale=Fr(1)):
    out=dict(s1)
    for k,v in s2.items():
        out[k]=out.get(k,Fr(0))+s2scale*v
        if out[k]==0: del out[k]
    return out
def sscale(s,a):
    a=Fr(a); return {k:v*a for k,v in s.items() if v*a}
def tensor_right(s, letter):
    # append letter (possibly formal sum, given as dict {vec: coeff}) at the END
    out={}
    for k,c in s.items():
        for lv,lc in letter.items():
            nk=k+(lv,)
            out[nk]=out.get(nk,Fr(0))+c*lc
    return {k:v for k,v in out.items() if v}
def L(v): return {v:Fr(1)}   # single letter

# --- base symbols from S(Li_n(z)) = -[1-z]^[z]^(n-1) ---
def S_Lin(zvec, uzvec, n):
    # returns -[1-z]^[z]^...^[z]
    s={(uzvec,):Fr(-1)}
    for _ in range(n-1):
        s=tensor_right(s, L(zvec))
    return s

# dlog recursion -> symbol recursion:
# dF = sum_i F_i dlog g_i  =>  S(F) = sum_i S(F_i) ^ [g_i]
def build_symbols():
    S={}
    # Li11: d = -Li1(xy)[dlog(1-x)+dlog y-dlog(1-y)] - Li1(x) dlog(1-y)
    # S(Li1(z)) = -[1-z]
    s={}
    s=tensor_right(sadd(s, sscale(S_Lin(Z,W,1),-1)), {U:Fr(1),Y:Fr(1),V:Fr(-1)})   # -S(Li1(xy))=+W ... careful: term is S(-Li1(xy))^{letters}
    # recompute directly: S(-Li1(xy)) = +[W]; S(-Li1(x)) = +[U]
    s={}
    s=tensor_right({(W,):Fr(1)}, {U:Fr(1),Y:Fr(1),V:Fr(-1)})     # -Li1(xy) part
    s=sadd(s, tensor_right({(U,):Fr(1)}, L(V)))                    # -Li1(x) dlog(1-y)
    S['Li11']=s
    # Li12: d = -Li2(xy) dlog(1-x) + Li11 dlog y
    s=tensor_right(S_Lin(Z,W,2), L(U), -1)     # S(-Li2(xy)) ^ [1-x]
    s=sadd(s, tensor_right(S['Li11'], L(Y)))
    S['Li12']=s
    # Li21: d = Li11 dlog x - Li2(xy) dlog(y/(y-1)) - Li2(x) dlog(1-y)
    # [y/(y-1)] = [y]-[1-y] (drop [-1])
    s=tensor_right(S['Li11'], L(X))
    s=sadd(s, tensor_right(S_Lin(Z,W,2), {Y:Fr(1),V:Fr(-1)}, -1))
    s=sadd(s, tensor_right(S_Lin(X,U,2), L(V), -1))
    S['Li21']=s
    # Li22: d = Li12 dlog x + Li21 dlog y
    s=tensor_right(S['Li12'], L(X))
    s=sadd(s, tensor_right(S['Li21'], L(Y)))
    S['Li22']=s
    # Li31: d = Li21 dlog x - Li3(xy) dlog(y/(y-1)) - Li3(x) dlog(1-y)
    s=tensor_right(S['Li21'], L(X))
    s=sadd(s, tensor_right(S_Lin(Z,W,3), {Y:Fr(1),V:Fr(-1)}, -1))
    s=sadd(s, tensor_right(S_Lin(X,U,3), L(V), -1))
    S['Li31']=s
    # Li13: d = -Li3(xy) dlog(1-x) + Li12 dlog y
    s=tensor_right(S_Lin(Z,W,3), L(U), -1)
    s=sadd(s, tensor_right(S['Li12'], L(Y)))
    S['Li13']=s
    # Li4(x), Li4(xy), Li4(y) for classical comparisons
    S['Li4x']=S_Lin(X,U,4); S['Li4xy']=S_Lin(Z,W,4); S['Li4y']=S_Lin(Y,V,4)
    S['Li2x']=S_Lin(X,U,2); S['Li2y']=S_Lin(Y,V,2)
    return S

# --- shuffle product of symbols (formal) ---
def shuffle_word(w1, w2):
    # all shuffles of tuples w1,w2
    if not w1: return [w2]
    if not w2: return [w1]
    out=[]
    for t in shuffle_word(w1[1:], w2): out.append((w1[0],)+t)
    for t in shuffle_word(w1, w2[1:]): out.append((w2[0],)+t)
    return out
def shuffle_sym(s1, s2):
    out={}
    for k1,c1 in s1.items():
        for k2,c2 in s2.items():
            for t in shuffle_word(k1,k2):
                out[t]=out.get(t,Fr(0))+c1*c2
    return {k:v for k,v in out.items() if v}

# --- swap x<->y relabeling (for shuffle checks): X<->Y, U<->V, W fixed ---
def swap_xy(s):
    m={X:Y, Y:X, U:V, V:U, W:W}
    out={}
    for k,c in s.items():
        nk=tuple(m[v] for v in k)
        out[nk]=out.get(nk,Fr(0))+c
    return {k:v for k,v in out.items() if v}

# --- wedge / delta_22 in the 5-dim two-variable space ---
WB5=[(i,j) for i in range(5) for j in range(i+1,5)]  # 10 wedge basis pairs
def wedge2_5(a,b):
    out={}
    for idx,(i,j) in enumerate(WB5):
        c=a[i]*b[j]-a[j]*b[i]
        if c: out[idx]=out.get(idx,0)+c
    return out
def wedge_bivec(A,B):
    # antisymmetric wedge of two bivectors: dict over (i,j), i<j, of the 10-dim Lambda^2
    out={}
    for i,ci in A.items():
        for j,cj in B.items():
            if i==j: continue
            if i<j: out[(i,j)]=out.get((i,j),0)+ci*cj
            else: out[(j,i)]=out.get((j,i),0)-ci*cj
    return {k:v for k,v in out.items() if v}
def delta22_5(s):
    out={}
    for (a,b,c,d),cf in s.items():
        w=wedge_bivec(wedge2_5(a,b), wedge2_5(c,d))
        for k,v in w.items():
            out[k]=out.get(k,Fr(0))+cf*v
            if out[k]==0: del out[k]
    return out
'''
with open('cor8br_sym.py','w') as f:
    f.write(engine)
print("engine written")
