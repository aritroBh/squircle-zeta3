# B4: mutation tests. Baseline must give 0 bad; each mutation must be DETECTED (bad>0) by the chain comparison.
import sys, time, pickle
sys.dont_write_bytecode=True
from core import *
import core, cand
from math import comb
H = pickle.load(open('<project>/research/08_theory/TOWN1_tower_large_n/data/towers_D8_K52.pkl','rb')); H.pop('_meta')
def toG(q):
    return G(fr(q[1]), fr(q[3]))
KS = list(range(12, 21))    # levels
def T1_param(tau1, tau2, binom_fn=None, conj=True, foldsign=1):
    def T1(i,j,delta,U):
        d=i+j; tot=ZERO
        for dp in range(delta+1):
            den=2*d-dp-1
            tot = tot + tau1[delta-dp]*U((i+1,j),dp)*fmpq(i+1,den)
            if i>j: tot = tot + tau2[delta-dp]*U((i,j+1),dp)*fmpq(j+1,den)
        if j==i-1:
            m=4*i-4-delta; ft=ZERO
            for mu in range(m,4*i-3):
                t=4*i-4-mu; g=ZERO
                for dp in range(t+1):
                    g = g - U((i,i),dp)*tau1[t-dp]*fmpq(1,4*i-3-dp)
                gc = g.conj() if conj else g
                bf = binom_fn or comb
                ft = ft + gc*ipow(-mu)*ipow(m-2)*bf(mu-2,mu-m)*foldsign
            tot = tot + ft*(-i)
        return tot
    return T1

def run(name, tau1=None, tau2=None, binom_fn=None, conj=True, foldsign=1, rho_fn=None, cand_kw=None, cand_fn=None, e3_mut=None):
    cand.reset()
    tau1 = tau1 or TAU1; tau2 = tau2 or TAU2
    T1 = T1_param(tau1, tau2, binom_fn, conj, foldsign)
    rf = rho_fn or rho_tilde
    bad_t1=0; bad_c=0; n=0; bad_c_by_delta={}
    badt1_by_delta={}
    for k in KS:
        T = tower_levelk(k, max(6,k-8), rho_fn=rf, T1=T1)
        for key, vals in T.items():
            p,q=key
            if p+q==k: continue
            for dl,u in enumerate(vals):
                n+=1
                pk=toG(H[key][dl][k])
                if not (pk==u):
                    bad_t1+=1; badt1_by_delta[dl]=badt1_by_delta.get(dl,0)+1
                kw = cand_kw or {}
                try:
                    c = (cand_fn or cand.cand_value)(p,q,dl,k,**kw)
                except ValueError:
                    c = None
                if c is None or not (c==pk):
                    bad_c+=1; bad_c_by_delta[dl]=bad_c_by_delta.get(dl,0)+1
    cand.reset()
    print(f'{name:55s} cells={n} T1-vs-data bad={bad_t1} {dict(sorted(badt1_by_delta.items()))} | cand-vs-data bad={bad_c} {dict(sorted(bad_c_by_delta.items()))}', flush=True)
    return bad_t1, bad_c

def shifted(t, n, delta_val):
    t2=list(t); t2[n] = t2[n] + delta_val; return t2

print('--- baseline')
run('baseline')
print('--- tau mutations')
run('tau1_8 += 1/128 (T1 and cand)', tau1=shifted(TAU1,8,G(fmpq(1,128))), cand_kw=dict(tau1=shifted(TAU1,8,G(fmpq(1,128)))))
run('tau2_8 conj-sign flipped (T1 and cand)', tau2=shifted(TAU2,8,G(0,fmpq(1,64))), cand_kw=dict(tau2=shifted(TAU2,8,G(0,fmpq(1,64)))))
run('tau1_7 += i/64 (T1 and cand)', tau1=shifted(TAU1,7,G(0,fmpq(1,64))), cand_kw=dict(tau1=shifted(TAU1,7,G(0,fmpq(1,64)))))
run('tau1_6 += 1/32 (T1 and cand)', tau1=shifted(TAU1,6,G(fmpq(1,32))), cand_kw=dict(tau1=shifted(TAU1,6,G(fmpq(1,32)))))
print('--- fold mutations')
run('fold binom(mu-2,mu-m) -> binom(mu-1,mu-m)', binom_fn=lambda n_,k_: comb(n_+1,k_))
run('fold binom -> binom(mu-2,mu-m-1)... shifted lower index', binom_fn=lambda n_,k_: comb(n_,k_+1) if k_+1<=n_ else 0)
run('fold: no complex conjugation', conj=False)
run('fold: sign flipped', foldsign=-1)
print('--- initial-data mutations (T1 run from mutated rho; cand base mutated)')
def rho_mut_a8(p,q,N=8):
    r = rho_tilde(p,q,N); d=p+q
    r = list(r); r[8] = r[8] + rho0_tilde(p,q)*G(fmpq(1,1000))   # a_8 -> a_8 + 1/1000
    return r
run('rho(8) += rho(0)/1000 (initial data, T1+cand base)', rho_fn=rho_mut_a8)   # cand uses true a_8 -> cand mismatch at delta=8 only via data cmp
def rho_mut_a7(p,q,N=8):
    r = list(rho_tilde(p,q,N)); r[7] = r[7] + rho0_tilde(p,q)*G(0,fmpq(1,1000)); return r
run('rho(7) += i rho(0)/1000', rho_fn=rho_mut_a7)
# candidate with mutated base a_8 / a_7 (closing then fails): base returns list
def base_mut(delta_m, eps):
    def base(D,M):
        a=list(a_of(D,M)); a[delta_m] = a[delta_m] + eps*fmpq(1,1)/ (2*D-3)   # perturb by small rational function
        return a
    return base
run('cand base a_8 perturbed by eps/(2d-3)', cand_kw=dict(base=base_mut(8,G(1))))
run('cand base a_7 perturbed by eps/(2d-3)', cand_kw=dict(base=base_mut(7,G(0,1))))
def base_swap(D,M): return a_of(D, D-M)
run('cand base a_delta(d, d-m) (m <-> q swap)', cand_kw=dict(base=base_swap))
print('--- basis / range mutations')
def cand_basis_off(p,q,dl,k,**kw): return cand.cand_value(p,q,dl,k-1) if True else None
run('cand evaluated at k-1 (off-by-one level)', cand_fn=lambda p,q,dl,k,**kw: cand.cand_value(p,q,dl,k-1))
run('cand evaluated at k+1', cand_fn=lambda p,q,dl,k,**kw: cand.cand_value(p,q,dl,k+1))
