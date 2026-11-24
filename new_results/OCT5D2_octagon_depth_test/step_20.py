# Chen concatenation: G(w;1) = sum_{uv=w} G(u;1/2) (-1)^{|v|} G(mrev(v);1/2)
import mpmath as mp
mp.mp.dps = 170
t0=time.time()
def mrev(v):
    return tuple(mp.mpf(1)-x for x in reversed(v))
G1 = {}
for w in Ghalf:
    tot = mp.mpf(0)
    for i in range(len(w)+1):
        u, v = w[:i], w[i:]
        tot += Ghalf[u] * ((-1)**len(v)) * Ghalf[mrev(v)]
    G1[w] = tot
print("Chen done in %.1fs, %d words" % (time.time()-t0, len(G1)))

# --- validation 1: polylogs  G(0^{k-1}, a; 1) = -Li_k(1/a)
for a, an in [(LP,'p'), (LQ,'q')]:
    for k in range(2, 6):
        w = (L0,)*(k-1) + (a,)
        val = G1[w]
        ref = -mp.polylog(k, 1/a)
        print(f"G(0^{k-1},{an};1) + Li_{k}(1/{an}) =", mp.nstr(val-ref, 4))
# --- validation 2: direct quad of convergent words (closed-form inner integrals)
tests = []
tests.append(((LP,LQ), mp.quad(lambda t: mp.log(1-t/LQ)/(t-LP), [0,1])))
tests.append(((LQ,LP), mp.quad(lambda t: mp.log(1-t/LP)/(t-LQ), [0,1])))
tests.append(((LP,L1), mp.quad(lambda t: mp.log(1-t)/(t-LP), [0,1])))
tests.append(((L0,LP), -mp.polylog(2, 1/LP)))
tests.append(((L0,LP,LQ), mp.quad(lambda t: mp.quad(lambda s: 1/(s-LQ), [0,t])/(t-LP), [0,1])))
for w, ref in tests:
    print("G%s;1 vs direct:" % (w,), mp.nstr(G1[w]-ref, 4))
