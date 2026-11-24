import sys
sys.path.insert(0,'.')
import point_engine as PE, mzv_exact as M
from fractions import Fraction as F

def f_coord(a, W, idx):
    """single irreducible coordinate (regular part) at Fraction point a."""
    a=[F(x) for x in a]
    vw = PE.volreg_words_point(a, W); vwr = PE.volreg_words_point(a[::-1], W)
    T = PE.table(W); rw_acc=F(0)
    for src in (vw[W], vwr[W]):
        for w,pp in src.items():
            c=T.reduce_word(w)[idx]
            if c: rw_acc+=pp*c
    return rw_acc

def f11(a): return f_coord(a, 11, 8)

# Test divisibility: f11 should vanish when a2=0, a3=0, a4=0, or Sum a = 0
import itertools
tests = {
 'a2=0':  [(1,0,2,3,1),(2,0,1,1,3),(1,0,1,1,1),(3,0,2,1,2)],
 'a3=0':  [(1,2,0,3,1),(2,1,0,1,3),(1,1,0,1,1),(2,3,0,2,1)],
 'a4=0':  [(1,2,3,0,1),(2,1,1,0,3),(1,1,1,0,1),(3,1,2,0,2)],
 'Sum=0': [(1,2,3,4,-10),(1,1,1,1,-4),(2,1,3,1,-7),(1,3,1,2,-7)],
}
for name, pts in tests.items():
    vals=[f11(p) for p in pts]
    print(f"{name}: {vals}  all zero: {all(v==0 for v in vals)}")
