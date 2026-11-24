import time, sys, json, os
sys.path.insert(0,'<project>/research/08_theory/DCYA1_corner_decay_proof/')
import mpmath as mp
import corner_series as cs
# Ground truth M_k via corner_coeffs at high dps
mp.mp.dps=100
t=time.time()
cc=cs.corner_coeffs(mp.pi/2, 55)
print("corner_coeffs N=55 dps=100 time:",round(time.time()-t,1),"s")
# cache as strings
json.dump([mp.nstr(x,80) for x in cc], open('mk_pi2_k55.json','w'))
# Build series S(u) and compare to M_split at hard points
def gamma_m(c, M):
    Q=[mp.mpf(0)]*(M+1); Q[0]=1
    if M>=1: Q[1]=-2*c
    if M>=2: Q[2]=1
    b=[mp.sqrt(Q[0])]
    for n in range(1,M+1):
        s=Q[n]-mp.fsum(b[k]*b[n-k] for k in range(1,n))
        b.append(s/(2*b[0]))
    return b
def beta_j(u, J):
    out=[mp.mpf(1)]
    for j in range(1,J+1): out.append(-(u+j)/j*out[-1])
    return out
def M_split2(u, a, Jm=30):
    c=mp.cos(a); g0=mp.mpf(1)/64
    gm=gamma_m(c,Jm); bj=beta_j(u,Jm)
    ana=mp.mpf(0)
    for j in range(Jm+1):
        for m in range(Jm+1):
            if (j,m)==(0,0): continue
            nu=j+(m-1)*(1-u); ana+=bj[j]*gm[m]*g0**(nu+1)/(nu+1)
    def F(g):
        return ((1+g)**(-1-u)*mp.sqrt(1+g**(2-2*u)-2*c*g**(1-u))-1)*g**(u-1)
    return ana+mp.quad(F, [g0, mp.mpf(1)/8, mp.mpf(1)/2, 1])
def S_series(u):
    return mp.fsum(cc[k]*u**k for k in range(len(cc)))
for u in [mp.mpf('-0.7'),mp.mpf('-0.8'),mp.mpf('-0.85'),mp.mpf('-0.9')]:
    ms=M_split2(u,mp.pi/2); ss=S_series(u)
    # tail estimate: next term cc[55]*u^55 roughly
    print(f"u={u}: M_split-S = {mp.nstr(abs(ms-ss),3)}  (tail~{mp.nstr(abs(cc[55]*u**55),2)})")
