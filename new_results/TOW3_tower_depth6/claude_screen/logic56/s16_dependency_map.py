# s16: which non-adjacent Newton identities N(delta,s) (symbolic in d,m) does each fitted form ENTER?  Perturb one form by a generic
# multiple of a generic polynomial/Pi and list the identities whose symbolic residual becomes nonzero.
import sys, time
sys.dont_write_bytecode = True
import sympy as sp
import tow2lib as T2, tow3lib as T3
d, m = T3.d, T3.m; I = sp.I
tau1 = T2.tau1f; tau2 = T2.tau2f
base = {(dl, s): T3.e(dl, s) for dl in range(0, 7) for s in range(0, dl//2+1)}
def make_resid(F):
    def Fz(dl, s): return F[(dl, s)] if (dl, s) in F else sp.Integer(0)
    def lam(delta, s):
        tot = 0
        for dp in range(delta+1):
            fA = Fz(dp, s).subs({d: d+1, m: m+1}, simultaneous=True); fB = Fz(dp, s).subs({d: d+1, m: m}, simultaneous=True)
            tot += I*(2*d-1)/(2*d-dp-1)*(tau1(delta-dp)*fA - tau2(delta-dp)*fB)
        return tot
    def resid(delta, s):
        return sp.cancel(sp.together(Fz(delta, s) + Fz(delta, s+1) - lam(delta, s)))
    return resid
t0 = time.time()
r0 = make_resid(base)
bad0 = [(dl, s) for dl in range(3, 7) for s in range(0, dl//2+1) if r0(dl, s) != 0]
print('unperturbed: nonzero non-adjacent identities among delta=3..6:', bad0, f'({round(time.time()-t0,1)}s)', flush=True)
pert = {(4,1): (d**2+m+1)/T3.Pi(2), (5,1): (d**2+m+1)/T3.Pi(2), (6,1): (d**2+m+1)/T3.Pi(3), (6,2): (d**2+m+1)/T3.Pi(3),
        (5,2): (d-1)/T3.Pi(2), (6,3): (d+7)/T3.Pi(3)}
for key, x in pert.items():
    F = dict(base); F[key] = F[key] + x
    r = make_resid(F)
    hit = [(dl, s) for dl in range(2, 7) for s in range(0, dl//2+1) if r(dl, s) != 0]
    print(f'  perturb e{key}: non-adjacent identities that FAIL (= identities it enters): {hit}   ({round(time.time()-t0,1)}s)', flush=True)
