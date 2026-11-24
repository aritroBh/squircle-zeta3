import os, sys, time
os.chdir("<project>/research/08_theory/OCT5D2_octagon_depth_test")
sys.path.insert(0, ".")
import gengine, mpmath as mp
G1, letters = gengine.gtab(170, 480)
L0,L1,LP,LQ = letters['L0'],letters['L1'],letters['LP'],letters['LQ']

print("--- polylogs: G(0^{k-1},a;1) = -Li_k(1/a)")
for a, an in [(LP,'p'), (LQ,'q')]:
    for k in range(2, 6):
        w = (L0,)*(k-1) + (a,)
        print(f"  k={k} {an}:", mp.nstr(G1[w]+mp.polylog(k, 1/a), 4))
print("--- direct quadrature checks")
tests = [((LP,LQ), mp.quad(lambda t: mp.log(1-t/LQ)/(t-LP), [0,1])),
         ((LQ,LP), mp.quad(lambda t: mp.log(1-t/LP)/(t-LQ), [0,1])),
         ((LP,L1), mp.quad(lambda t: mp.log(1-t)/(t-LP), [0,1])),
         ((LQ,L1), mp.quad(lambda t: mp.log(1-t)/(t-LQ), [0,1])),
         ((L0,LP), mp.quad(lambda t: mp.log(1-t/LP)/t, [0,1]))]
for w, ref in tests:
    nm = ''.join('0' if x==0 else ('1' if x==1 else ('p' if x==LP else 'q')) for x in w)
    print("  G(%s;1) residual:" % nm, mp.nstr(G1[w]-ref, 4))

print("--- nested-quad checks of weight-3 mixed words")
# G(p,0,q;1) = int_0^1 dt/(t-p) int_0^t ds/(s-0) G(q;s)... careful: G(0,q;t)=int_0^t ds/s log(1-s/q)
v = mp.quad(lambda t: mp.quad(lambda s: mp.log(1-s/LQ)/s, [0,t])/(t-LP), [0,1])
print("  G(p,0,q;1) residual:", mp.nstr(G1[(LP,L0,LQ)]-v, 4))
v = mp.quad(lambda t: mp.quad(lambda s: mp.log(1-s/LP)/s, [0,t])/(t-LQ), [0,1])
print("  G(q,0,p;1) residual:", mp.nstr(G1[(LQ,L0,LP)]-v, 4))

print("--- shuffle identities (|u|+|v|<=5)")
import random
random.seed(7)
def shuf_words(u, v):
    if not u: return {v:1}
    if not v: return {u:1}
    out = {}
    for w,c in shuf_words(u[1:], v).items(): out[(u[0],)+w] = out.get((u[0],)+w,0)+c
    for w,c in shuf_words(u, v[1:]).items(): out[(v[0],)+w] = out.get((v[0],)+w,0)+c
    return out
pool = [w for w in G1 if 1 <= len(w) <= 4 and w[0] != L1]
worst = mp.mpf(0)
for _ in range(60):
    u = random.choice(pool); v = random.choice(pool)
    if len(u)+len(v) > 5: continue
    lhs = G1[u]*G1[v]
    rhs = sum(c*G1[w] for w,c in shuf_words(u,v).items())
    worst = max(worst, abs(lhs-rhs))
print("  worst shuffle residual:", mp.nstr(worst, 4))
print("--- regularized trailing-zero values: G_*(p,0;1) should = -G(0,p;1) via shuffle G(p)G(0)")
print("  G(p,0;1) + G(0,p;1):", mp.nstr(G1[(LP,L0)]+G1[(L0,LP)], 4), " vs G(p;1)G(0;1)=reg 0 check:",
      mp.nstr(G1[(LP,)]*G1[(L0,)], 4))
