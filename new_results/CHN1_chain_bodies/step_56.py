import os, time
WORK='<project>/research/08_theory/CHN1_chain_bodies'
os.chdir(WORK)
t0=time.time()
# add flint-based rref to mzv_exact
add=r'''

# ---------- flint-based rref for large weights ----------
def rref_flint(rels, idx):
    from flint import fmpq_mat, fmpq
    n=len(idx)
    rows=[]
    for r in rels:
        rows.append((r,))
    M=fmpq_mat(len(rels), n, [c for (r,) in rows for c in [r.get(w,0) for w in idx]])
    R, rank = M.rref()
    pivots=[]
    for i in range(rank):
        row={}
        lead=None
        for j in range(n):
            v=R[i,j]
            if v!=0:
                row[j]=F(v.p, v.q)
                if lead is None: lead=j
        pivots.append((lead,row))
    pivots.sort()
    return pivots

class DSTableFlint(DSTable):
    def __init__(self, W):
        self.W=W; self.words=WORDS[W]; self.n=len(self.words)
        self.idx={w:i for i,w in enumerate(self.words)}
        rels=relations(W)
        self.pivots=rref_flint(rels, self.idx)
        self.freecols=[j for j in range(self.n) if all(col!=j for col,_ in self.pivots)]
        self.nullity=len(self.freecols)
        self.monomials=BASIS[W]
        mv=[wordcombo(m) for m in self.monomials]
        red=[self._elim(v) for v in mv]
        d=self.nullity
        assert d==len(self.monomials), (W,d,len(self.monomials))
        T=[[red[i].get(fc,F(0)) for fc in self.freecols] for i in range(d)]
        self.Tinv=_inv(T)
        self._wordcache={}
'''
src=open('mzv_exact.py').read()
src=src.replace("\ndef dual(w):", add+"\ndef dual(w):")
open('mzv_exact.py','w').write(src)
import importlib, mzv_exact as M
importlib.reload(M)
# validate flint rref against pure python at W=10
import engine as E
T10=M.DSTable(10)
T10f=M.DSTableFlint(10)
print('nullities:', T10.nullity, T10f.nullity)
same=all(T10.reduce_word(w)==T10f.reduce_word(w) for w in M.WORDS[10])
print('W=10 reductions identical:', same)
print('t=%.1fs'%(time.time()-t0))
