# Planner (W4): independent checks of CHN4 with CHN3's ORIGINAL point_engine (not Kimi's copy).
# (1) f11 == -(1/10) a2 a3 a4 (Sum a) P7 at fresh random rational points (exact).
# (2) On Gamma the W12 coordinates looked proportional: z2z37/z2^2z35 = 75/14? z2118/z2^2z35 = -25/3?
#     Test whether these are polynomial IDENTITIES (at generic rational points) or only hold on V10.
import sys, os, random, time
from fractions import Fraction as F
CHN3 = '<project>/research/08_theory/CHN3_p5_locus'
HERE = '<project>/research/08_theory/CHN4_p5_w11_locus'
os.chdir(CHN3); sys.path.insert(0, CHN3)
import point_engine as PE
import sympy as sp
a = sp.symbols('a1:6')
P7 = sp.Poly(sp.sympify(open(HERE + '/P7.txt').read().replace('a1','a1')), *a)
def coords(pt, W):
    vw = PE.volreg_words_point(pt, W); vwr = PE.volreg_words_point(pt[::-1], W)
    T = PE.table(W); acc = [F(0)] * T.nullity
    for src in (vw[W], vwr[W]):
        for w, pp in src.items():
            rw = T.reduce_word(w)
            for b in range(T.nullity):
                if rw[b]: acc[b] += pp * rw[b]
    return acc
random.seed(29)
for k in range(3):
    pt = tuple(F(random.randint(1, 40), random.randint(1, 13)) for _ in range(5))
    t0 = time.time(); c11 = coords(pt, 11)
    s = pt[0] + pt[1] + pt[2] + pt[3] + pt[4]
    pred = -F(1, 10) * pt[1] * pt[2] * pt[3] * s * F(P7.eval(dict(zip(a, pt))))
    c12 = coords(pt, 12)
    A, B, Cz9, Cz57, Z = c12[7], c12[8], c12[9], c12[10], c12[11]
    print(f'pt {k}: f11 match {c11[8] == pred}; z2z37/z2^2z35 = {float(B/A):.6f}; z2118/z2^2z35 = {float(Z/A):.6f}; '
          f'CLASS/z2^2z35 = {float((Cz9 + F(14, 9) * Cz57)/A):.6f} ({time.time()-t0:.0f}s)', flush=True)
