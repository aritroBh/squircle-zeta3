import pickle, mpmath as mp, time
from q1s_exact import femb, feq
t0=time.time()
mp.mp.dps=130
V3=pickle.load(open('q3_accs.pkl','rb'))['V3p']
# direct: weight-2 words via mp.quad, products via logs
def G2(r,s):
    # G(r,s;1) = int_0^1 log(1-t/s)/(t-r) dt
    if s==0:
        return mp.quad(lambda t: mp.log(t)/(t-r), [0,1])
    if r==0:
        return mp.quad(lambda t: mp.log(1-t/s)/t, [0,1])
    return mp.quad(lambda t: mp.log(1-t/s)/(t-r), [0,1])
W2=mp.mpc(0)
for (w,lvm),c in V3.W.items():
    if len(w)+len(lvm)==3: continue
    wn=tuple(femb(s,1) for s in w)
    if len(w)==1:
        v=mp.log(1-1/wn[0])
    else:
        v=G2(wn[0],wn[1])
    for f in lvm: v*=mp.log(femb(f[1],1))
    W2+=femb(c,1)*v
P=mp.mpc(0)
for lvm,c in V3.P.items():
    v=mp.mpc(1)
    for f in lvm: v*=mp.log(femb(f[1],1))
    P+=femb(c,1)*v
W2P_direct=mp.re(W2)+mp.re(P)
print("W2+P direct (130d):", mp.nstr(W2P_direct,25))
Kold=pickle.load(open('K2_V3p_400.pkl','rb')); Rold=pickle.load(open('R_V3p_400.pkl','rb'))
W2P = mp.mpf(Kold['re'])-mp.mpf(Rold['re'])
print("W2+P from K-R   :", mp.nstr(W2P,25))
print("difference:", mp.nstr(abs(W2P_direct-W2P),6))
print("time %.0fs"%(time.time()-t0))
