import sys, pickle, time
sys.path.insert(0, "."); sys.path.insert(0, "<project>/research/08_theory/E5c_nonlocal_z5")
import importlib, mpmath as mp
mp.mp.dps = 50
import dod5d_exact as DX
importlib.reload(DX)
r3 = mp.sqrt(3)
KAPm = [mp.mpf(str(k[0]))+mp.mpf(str(k[1]))*r3 for k in DX.KAP]
C0m  = [mp.mpf(str(c[0]))+mp.mpf(str(c[1]))*r3 for c in DX.C0S]
def eval_d(d, t):
    LT=mp.log(t); LU=mp.log(1-t)
    W=[mp.log(1+KAPm[i]*t) for i in range(4)]
    CA=mp.log((2-r3)/2); CC0=mp.log(2-r3); CC1=mp.log((2*r3-3)/2)
    logs=[LT,LU,W[0],W[1],W[2],W[3],CA,CC0,CC1]
    tot=mp.mpf(0)
    for (n,le),c in d.items():
        v=t**n
        for i,e in enumerate(le):
            if e: v*=logs[i]**e
        tot+=(mp.mpf(str(c[0]))+mp.mpf(str(c[1]))*r3)*v
    return tot
t = mp.mpf('0.37')
ap=[C0m[i]*KAPm[i] for i in range(4)]
a=[C0m[i]*(1+KAPm[i]*t) for i in range(4)]
CCv = {7: mp.log(2-r3), 8: mp.log((2*r3-3)/2), 6: mp.log((2-r3)/2)}
la=[CCv[DX.CCslot[i]]+mp.log(1+KAPm[i]*t) for i in range(4)]
for (kind,k,e) in [('gam',1,1),('gam',1,2),('gam',4,1),('gam',2,2),('tau',1,3),('tau',2,2),('tau',5,1)]:
    d = DX.taugam(kind,k,e)
    v = eval_d(d, t)
    ref = (sum(a[i]*la[i]**k for i in range(4)) if kind=='tau' else sum(ap[i]*la[i]**k for i in range(4)))**e
    print(kind,k,e, mp.nstr(abs(v-ref),3))
# full re-expand
terms = pickle.load(open("E5_terms.pkl","rb"))
t0 = time.time()
E5X = DX.expand_E5(terms)
print("expand_E5: %.1fs, monomials: %d" % (time.time()-t0, len(E5X)))
pickle.dump(E5X, open("E5_dod_expanded.pkl","wb"))
