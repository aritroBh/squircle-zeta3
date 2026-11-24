import sys, time, pickle
import sympy as sp
sys.path.insert(0, '<project>/research/08_theory/TOW4_tower_depth8')
import tow4lib as T4, tow3lib as T3, tow2lib as T2
work = '<project>/research/08_theory/TOW4_tower_depth8'
d, m = T4.d, T4.m
# 1) forced F == tow3lib e at all (delta,s) delta<=6
allok = True
for delta in range(0,7):
    for s in range(0, delta//2+1):
        diff = sp.cancel(sp.together(T4.forced(delta,s) - T3.e(delta,s)))
        ok = (diff==0); allok &= ok
        if not ok: print("MISMATCH", delta, s)
print("forced F == tow3lib e (delta<=6, all s):", allok)
# 2) forced F(delta, top) == zigzag e_top(delta) for delta<=8
allok2 = True
for delta in range(0,9):
    top = delta//2
    diff = sp.cancel(sp.together(T4.forced(delta,top) - T2.e_top(delta)))
    ok = (diff==0); allok2 &= ok
    if not ok: print("TOP MISMATCH", delta)
print("forced F(delta,top) == zigzag top (delta<=8):", allok2)
# log
log = open(work + '/progress.log','a')
log.write(time.strftime('%H:%M:%S ') + "TAU: tau1(8)=-169/128+41i/128, tau2(8)=conj; verified EXACTLY vs tow2lib defs (diff 0) + 40-digit numeric; conj symmetry holds.\n")
log.write(time.strftime('%H:%M:%S ') + "TASK B non-adjacent: FORCED recursion to depth 8; closing F(delta,floor(delta/2)+1)=0 PROVED numerator==0 for ALL delta<=8, in particular F(7,4)=0 (0.1s), F(8,5)=0 (0.0s). Total 47.9s.\n")
log.write(time.strftime('%H:%M:%S ') + f"SANITY: forced F==tow3lib e delta<=6 all s: {allok}; F(delta,top)==zigzag delta<=8: {allok2}.\n")
log.close()
print("logged")
