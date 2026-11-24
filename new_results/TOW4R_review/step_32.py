
import os, sys, pickle, time
import sympy as sp
from fractions import Fraction as Fr
base = "<project>"
work = base + "/research/08_theory/TOW4R_review"
sys.path.insert(0, base + "/research/08_theory/TOW2_tower_all_depths")  # towlib, fastred, genintegrand
sys.path.insert(0, work + "/tow4_copy")
import towlib
import tow4lib as T4
d_, m_ = T4.d, T4.m
t0 = time.time()
# keys to test: include off-diagonal, adjacent, diagonal, fresh high degrees
keys = [(9,0),(7,3),(10,5),(12,1),(26,26),(20,20),(15,8),(31,7),(40,12),(52,0),(33,33),(25,25),(6,6),(7,7),(8,8)]
a7 = T4.arho(7); a8 = T4.arho(8)
bad = 0; n = 0
for (p,q) in keys:
    dd = p+q
    if 2*dd-4 < 8: continue
    for delta, aa in [(7,a7),(8,a8)]:
        r0 = towlib.rho((p,q), 0); rd = towlib.rho((p,q), delta)
        # QSI -> tuple (p0,p1,p2,p3): a + b sqrt2 + i(c + d sqrt2)?
        t0q = towlib.q4(r0); t1q = towlib.q4(rd)
        # divide rd/r0 in Q(sqrt2,i): representation (a,b,c,d) = a + b*r2 + i*(c + d*r2)
        def gm(A,B): return (A[0]*B[0]-A[1]*B[1], A[0]*B[1]+A[1]*B[0])
        def mul(x,y):  # (a + i b)(c + i d) in Q(r2) parts
            A=(x[0],x[2]); B=(x[1],x[3]); C=(y[0],y[2]); D=(y[1],y[3])
            return None
        # careful: tuple layout unknown; print once
        if n == 0:
            print("rho((9,0),0) tuple:", t0q)
        n += 1
print("n cells attempted:", n, " time", round(time.time()-t0,1))
