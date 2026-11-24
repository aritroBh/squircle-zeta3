fix = r'''
# cor8br_sym.py -- COR8BR's own exact symbol engine (written from definitions)
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
def tensor_right(s, letter, scale=Fr(1)):
    out={}
    for k,c in s.items():
        for lv,lc in letter.items():
            nk=k+(lv,)
            out[nk]=out.get(nk,Fr(0))+scale*c*lc
    return {k:v for k,v in out.items() if v}
def L(v): return {v:Fr(1)}

def S_Lin(zvec, uzvec, n):
    s={(uzvec,):Fr(-1)}
    for _ in range(n-1):
        s=tensor_right(s, L(zvec))
    return s

def build_symbols():
    S={}
    s=tensor_right({(W,):Fr(1)}, {U:Fr(1),Y:Fr(1),V:Fr(-1)})
    s=sadd(s, tensor_right({(U,):Fr(1)}, L(V)))
    S['Li11']=s
    s=tensor_right(S_Lin(Z,W,2), L(U), Fr(-1))
    s=sadd(s, tensor_right(S['Li11'], L(Y)))
    S['Li12']=s
    s=tensor_right(S['Li11'], L(X))
    s=sadd(s, tensor_right(S_Lin(Z,W,2), {Y:Fr(1),V:Fr(-1)}, Fr(-1)))
    s=sadd(s, tensor_right(S_Lin(X,U,2), L(V), Fr(-1)))
    S['Li21']=s
    s=tensor_right(S['Li12'], L(X))
    s=sadd(s, tensor_right(S['Li21'], L(Y)))
    S['Li22']=s
    s=tensor_right(S['Li21'], L(X))
    s=sadd(s, tensor_right(S_Lin(Z,W,3), {Y:Fr(1),V:Fr(-1)}, Fr(-1)))
    s=sadd(s, tensor_right(S_Lin(X,U,3), L(V), Fr(-1)))
    S['Li31']=s
    s=tensor_right(S_Lin(Z,W,3), L(U), Fr(-1))
    s=sadd(s, tensor_right(S['Li12'], L(Y)))
    S['Li13']=s
    S['Li4x']=S_Lin(X,U,4); S['Li4xy']=S_Lin(Z,W,4); S['Li4y']=S_Lin(Y,V,4)
    S['Li2x']=S_Lin(X,U,2); S['Li2y']=S_Lin(Y,V,2)
    S['Li3x']=S_Lin(X,U,3); S['Li3y']=S_Lin(Y,V,3)
    return S

def shuffle_word(w1, w2):
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

def swap_xy(s):
    m={X:Y, Y:X, U:V, V:U, W:W}
    out={}
    for k,c in s.items():
        nk=tuple(m[v] for v in k)
        out[nk]=out.get(nk,Fr(0))+c
    return {k:v for k,v in out.items() if v}

WB5=[(i,j) for i in range(5) for j in range(i+1,5)]
def wedge2_5(a,b):
    out={}
    for idx,(i,j) in enumerate(WB5):
        c=a[i]*b[j]-a[j]*b[i]
        if c: out[idx]=out.get(idx,0)+c
    return out
def wedge_bivec(A,B):
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

def show5(d):
    names=['X','Y','U','V','W']
    def bname(idx):
        i,j=WB5[idx]; return names[i]+names[j]
    parts=[]
    for (i,j),v in sorted(d.items()):
        parts.append(f"({v})*[{bname(i)}^{bname(j)}]")
    return " + ".join(parts) if parts else "0"
'''
with open('cor8br_sym.py','w') as f:
    f.write(fix)
import sys
sys.path.insert(0,'.')
import importlib, cor8br_sym as E
importlib.reload(E)
S = E.build_symbols()
for name in ['Li11','Li12','Li21','Li22','Li31','Li13']:
    print(name, "terms:", len(S[name]))
