import os, time
WORK='<project>/research/08_theory/CHN1_chain_bodies'
os.chdir(WORK)
t0=time.time()
src = r'''
"""Exact MZV toolkit: stuffle/shuffle, finite double-shuffle + Hoffman regularized
relations, reduction tables to a chosen monomial basis. All arithmetic exact (Fraction).
Convention: zeta(s1,...,sd) = sum_{n1>...>nd>=1} prod n_i^-s_i ; admissible = s1>=2."""
from fractions import Fraction
F = Fraction

def words_of_weight(W):
    res=[]
    def rec(rem, cur, first):
        if rem==0:
            res.append(tuple(cur)); return
        lo = 2 if first else 1
        for k in range(lo, rem+1):
            cur.append(k); rec(rem-k, cur, False); cur.pop()
    if W>=2: rec(W, [], True)
    return res
WORDS = {W: words_of_weight(W) for W in range(2,13)}
WT = lambda w: sum(w)

# ---------- stuffle (index notation) ----------
_ST={}
def stuffle(u,v):
    if u>v: u,v=v,u
    key=(u,v); r=_ST.get(key)
    if r is not None: return r
    if not u: r={v:1}
    elif not v: r={u:1}
    else:
        a,b=u[0],v[0]; r={}
        for w,c in stuffle(u[1:],v).items(): k=(a,)+w; r[k]=r.get(k,0)+c
        for w,c in stuffle(u,v[1:]).items(): k=(b,)+w; r[k]=r.get(k,0)+c
        for w,c in stuffle(u[1:],v[1:]).items(): k=(a+b,)+w; r[k]=r.get(k,0)+c
    _ST[key]=r; return r

# ---------- shuffle (xy alphabet) ----------
def to_xy(w): return ''.join('x'*(k-1)+'y' for k in w)
def from_xy(s):
    w=[]; c=0
    for ch in s:
        if ch=='x': c+=1
        else: w.append(c+1); c=0
    assert c==0, 'xy word must end with y'
    return tuple(w)
_SH={}
def _sh(s,t):
    key=(s,t); r=_SH.get(key)
    if r is not None: return r
    if not s: r={t:1}
    elif not t: r={s:1}
    else:
        r={}
        for w,c in _sh(s[1:],t).items(): k=s[0]+w; r[k]=r.get(k,0)+c
        for w,c in _sh(s,t[1:]).items(): k=t[0]+w; r[k]=r.get(k,0)+c
    _SH[key]=r; return r
def shuffle(u,v):
    if u>v: u,v=v,u
    return {from_xy(s):F(c) for s,c in _sh(to_xy(u),to_xy(v)).items()}

def wd_sub(a,b):
    d=dict(a)
    for w,c in b.items():
        d[w]=d.get(w,F(0))-c
        if d[w]==0: del d[w]
    return d
def wd_add_scaled(a,b,s):
    d=dict(a)
    for w,c in b.items():
        d[w]=d.get(w,F(0))+s*c
        if d[w]==0: del d[w]
    return d

def relations(W):
    """All proved relations at weight W: finite DS for admissible pairs + Hoffman w0=(1)."""
    rels=[]
    ws=WORDS
    for a in range(2, W):
        for u in ws[a]:
            for v in ws[W-a]:
                if u<=v:
                    d=wd_sub(shuffle(u,v), {w:F(c) for w,c in stuffle(u,v).items()})
                    if d: rels.append(d)
    for w in WORDS[W-1]:
        d=wd_sub(shuffle((1,),w), {ww:F(c) for ww,c in stuffle((1,),w).items()})
        assert all(ww[0]>=2 for ww in d), ('divergent word survived', W, w, d)
        if d: rels.append(d)
    return rels

# ---------- exact rref over Q ----------
def rref(rels, idx):
    rows=[]
    for r in rels:
        rows.append({idx[w]:F(c) for w,c in r.items()})
    pivots=[]  # (col, rowdict normalized)
    used=set()
    for r in rows:
        # reduce by existing pivots
        for col,pr in pivots:
            c=r.get(col)
            if c:
                for j,pc in pr.items():
                    r[j]=r.get(j,F(0))-c*pc
                    if r[j]==0: del r[j]
        if not r: continue
        col=min(r)
        inv=F(1)/r[col]
        for j in list(r): r[j]*=inv
        # eliminate this column from existing pivot rows
        for col2,pr in pivots:
            c=pr.get(col)
            if c:
                for j,rc in r.items():
                    if j==col: continue
                    pr[j]=pr.get(j,F(0))-c*rc
                    if pr[j]==0: del pr[j]
                pr.pop(col,None)
        pivots.append((col,r))
    pivots.sort()
    return pivots

# ---------- bases (task-specified monomials) ----------
BASIS = {
 2:[((2,),)],
 3:[((3,),)],
 4:[((4,),)],
 5:[((5,),),((2,),(3,))],
 6:[((6,),),((3,),(3,))],
 7:[((7,),),((2,),(5,)),((4,),(3,))],
 8:[((8,),),((3,),(5,)),((2,),(3,),(3,)),((3,5),)],
 9:[((9,),),((2,),(7,)),((4,),(5,)),((3,),(3,),(3,)),((2,),(2,),(2,),(3,))],
 10:[((10,),),((3,),(7,)),((5,),(5,)),((2,),(3,),(5,)),((4,),(3,),(3,)),((2,),(3,5)),((3,7),)],
 11:[((11,),),((2,),(9,)),((2,),(2,),(7,)),((2,),(2,),(2,),(5,)),((2,),(2,),(2,),(2,),(3,)),
     ((2,),(3,),(3,),(3,)),((3,),(3,),(5,)),((3,),(3,5)),((3,3,5),)],
}
BASIS[12]=[((12,),),((3,),(9,)),((5,),(7,)),((2,),(3,),(7,)),((2,),(5,),(5,)),((2,),(2,),(3,),(5,)),
           ((2,),(2,),(2,),(3,),(3,)),((3,),(3,),(3,),(3,)),((2,),(2,),(3,5)),((3,9),),((5,7),),((2,),(3,7))]

def wordcombo(mon):
    d={():F(1)}
    for L in mon:
        nd={}
        for w,c in d.items():
            for w2,c2 in stuffle(w,L).items():
                nd[w2]=nd.get(w2,F(0))+c*c2
        d=nd
    return d

class DSTable:
    def __init__(self, W):
        self.W=W; self.words=WORDS[W]; self.n=len(self.words)
        self.idx={w:i for i,w in enumerate(self.words)}
        rels=relations(W)
        self.pivots=rref(rels, self.idx)
        self.freecols=[j for j in range(self.n) if all(col!=j for col,_ in self.pivots)]
        self.nullity=len(self.freecols)
        self.monomials=BASIS[W]
        mv=[wordcombo(m) for m in self.monomials]
        # reduce each monomial vector onto free columns
        red=[self._elim(v) for v in mv]
        # T: matrix rows = monomial, cols = freecols -> must be square d x d
        d=self.nullity
        assert d==len(self.monomials), (W,d,len(self.monomials))
        T=[[red[i].get(fc,F(0)) for fc in self.freecols] for i in range(d)]
        self.Tinv=_inv(T)
        self._wordcache={}
    def _elim(self, v):  # v: dict word->F ; returns dict col->F supported on free cols
        r={self.idx[w]:F(c) for w,c in v.items()}
        for col,pr in self.pivots:
            c=r.get(col)
            if c:
                for j,pc in pr.items():
                    r[j]=r.get(j,F(0))-c*pc
                    if r[j]==0: del r[j]
        assert all(j in self.freecols for j in r), ('elimination left pivots', self.W, r)
        return r
    def reduce(self, v):
        r=self._elim(v)
        x=[r.get(fc,F(0)) for fc in self.freecols]
        d=self.nullity
        return [sum(x[k]*self.Tinv[k][i] for k in range(d)) for i in range(d)]
    def reduce_word(self, w):
        c=self._wordcache.get(w)
        if c is None:
            c=self.reduce({w:F(1)})
            self._wordcache[w]=c
        return c

def _inv(T):
    n=len(T)
    M=[row[:]+[F(1) if i==j else F(0) for j in range(n)] for i,row in enumerate(T)]
    for col in range(n):
        p=next(r for r in range(col,n) if M[r][col]!=0)
        M[col],M[p]=M[p],M[col]
        inv=F(1)/M[col][col]
        M[col]=[x*inv for x in M[col]]
        for r in range(n):
            if r!=col and M[r][col]!=0:
                c=M[r][col]
                M[r]=[x-c*y for x,y in zip(M[r],M[col])]
    return [row[n:] for row in M]

def dual(w):
    s=to_xy(w)
    t=s[::-1].translate(str.maketrans('xy','yx'))
    return from_xy(t)

if __name__=='__main__':
    import time
    for W in range(2,9):
        t=time.time()
        T=DSTable(W)
        print('W=%d words=%d nullity=%d (expected %d)  %.1fs'%(W,T.n,T.nullity,len(BASIS[W]),time.time()-t))
'''
open('mzv_exact.py','w').write(src)
print('written', time.time()-t0)
