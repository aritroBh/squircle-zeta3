"""W5C2 exact Galois-sector fold/unfold, replicating W5C1 w5_sector_fold.py exactly.
All arithmetic exact (ints for rows, Fractions for D)."""
import sys, os, pickle
sys.dont_write_bytecode=True
from fractions import Fraction as F
from math import gcd

HERE = os.path.dirname(os.path.abspath(__file__))
W5C1 = os.path.join(HERE, '..', 'W5C1_certificate_scoping')
for pth in (os.path.join(HERE,'..','A3c_claude_elementary','Q'),
            '<project>/paper/squircle_note/scripts'):
    if pth not in sys.path: sys.path.insert(0, pth)

SECTORS = [(1,1),(1,-1),(-1,1),(-1,-1)]

def tmap(w): return tuple(-1 if x==-1 else (-x)%8 for x in w)
def smap(w): return tuple(-1 if x==-1 else (3*x)%8 for x in w)
def gapply(w, g):
    if g == 0: return w
    if g == 1: return tmap(w)
    if g == 2: return smap(w)
    return tmap(smap(w))

def chi(g, sec): return (sec[0] if g in (1,3) else 1) * (sec[1] if g in (2,3) else 1)

class Ctx:
    def __init__(self):
        from relations import words
        self.W5 = words(5); N = len(self.W5)
        self.IDS = {w:i for i,w in enumerate(self.W5)}
        ti = [self.IDS[tmap(w)] for w in self.W5]
        si = [self.IDS[smap(w)] for w in self.W5]
        sti = [ti[si[i]] for i in range(N)]
        self.ti, self.si, self.sti = ti, si, sti
        oid = [-1]*N
        orbits = []
        for w in range(N):
            if oid[w] >= 0: continue
            k = len(orbits)
            ow = [(w,0),(ti[w],1),(si[w],2),(sti[w],3)]
            for x,g in ow: oid[x] = k
            orbits.append(ow)
        wm_o = [-1]*N; wm_g = [-1]*N; stab = []
        for k, ow in enumerate(orbits):
            repword = min(x[0] for x in ow)
            members = {}
            for cand, gg in ((repword,0),(ti[repword],1),(si[repword],2),(sti[repword],3)):
                members[cand] = gg
            for x in members:
                wm_o[x] = k; wm_g[x] = members[x]
            stabset = set()
            for gg, f in enumerate((lambda z:z, lambda z: ti[z], lambda z: si[z], lambda z: sti[z])):
                if f(repword) == repword: stabset.add(gg)
            stab.append(stabset)
        self.orbits, self.wm_o, self.wm_g, self.stab = orbits, wm_o, wm_g, stab
        self.N = N
        self.sec_cols = {}
        self.colmap = {}
        for sec in SECTORS:
            cols = [k for k in range(len(orbits)) if all(chi(g,sec)==1 for g in stab[k])]
            self.sec_cols[sec] = cols
            self.colmap[sec] = {k:i for i,k in enumerate(cols)}
    def check_vs_meta(self):
        meta = pickle.load(open(os.path.join(W5C1,'sector_meta.pkl'),'rb'))
        assert self.N == len(meta['wm_o'])
        assert list(meta['wm_o']) == list(self.wm_o)
        assert list(meta['wm_g']) == list(self.wm_g)
        assert [sorted(s) for s in meta['stab']] == [sorted(s) for s in self.stab]
        for sec in SECTORS:
            assert meta['sec_cols'][sec] == self.sec_cols[sec]
        return True

def sectorize_raw(row, sec, ctx):
    """row: iterable of (word_id, int coeff). Returns dict sector_col -> int (raw fold)."""
    cm = ctx.colmap[sec]; wm_o = ctx.wm_o; wm_g = ctx.wm_g
    d = {}
    for c, a in row:
        k = wm_o[c]
        if k not in cm: continue
        kk = cm[k]
        d[kk] = d.get(kk, 0) + chi(wm_g[c], sec)*a
    return {k:v for k,v in d.items() if v}

def norm(d):
    """-> (tuple_sorted, lam) with raw = lam * normalized."""
    if not d: return (), 1
    g = 0
    for v in d.values(): g = gcd(g, abs(v))
    items = sorted(d.items())
    if items[0][1] < 0: g = -g
    return tuple((k, v//g) for k, v in items), g

def fold_D(dct, sec, ctx):
    """dct: word_id -> Fraction. Returns dict sector_col -> Fraction."""
    cm = ctx.colmap[sec]; wm_o = ctx.wm_o; wm_g = ctx.wm_g
    d = {}
    for c, a in dct.items():
        k = wm_o[c]
        if k not in cm: continue
        kk = cm[k]
        d[kk] = d.get(kk, F(0)) + chi(wm_g[c], sec)*a
    return {k:v for k,v in d.items() if v}
