code = r'''# step_07.py: INDEPENDENT high-precision P^(k)(1), k=2,3,4,5, for the sheared parallelogram
# vertices +-(2,0), +-(1,2), power convention.  No use of the symbolic E-series:
# for each t, g_t(d) = |x_{1+d}'(t)| is analytic in d on |d|<1 (F-zeros have Re p = 0, p=1+d),
# so Cauchy on |d|=0.5 with N nodes gives g^{(k)}(0); then tanh-sinh over t.
import mpmath as mp, time, json, sys
mp.mp.dps = int(sys.argv[1]) if len(sys.argv) > 1 else 70

v0 = (mp.mpf(0), mp.mpf('0.5')); v1 = (mp.mpf('-0.5'), mp.mpf('0.25'))
def edge_geom(V0, V1):
    L = mp.sqrt((V1[0]-V0[0])**2 + (V1[1]-V0[1])**2)
    tau = ((V1[0]-V0[0])/L, (V1[1]-V0[1])/L)
    h = abs(V0[0]*V1[1] - V0[1]*V1[0])/L
    s0 = V0[0]*tau[0] + V0[1]*tau[1]
    return V0, V1, L, tau, h, s0

N = 220; R = mp.mpf('0.5')
om = [mp.exp(2*mp.pi*mp.j*k/N) for k in range(N)]
def gk_of_t(V0, V1, L, tau, h, s0, t, kmax=5):
    """return [g^{(k)}(0)/k! for k=0..kmax] of |x_{1+d}'(t)| at this t (E_k(t) values)"""
    y0 = (1-t)*V0[0] + t*V1[0]; y1 = (1-t)*V0[1] + t*V1[1]
    sig = s0 + L*t
    def g(d):
        p = 1 + d
        a0 = t; a1 = 1-t
        l0 = mp.log(a0) if a0 > 0 else mp.mpf('-inf')
        l1 = mp.log(a1) if a1 > 0 else mp.mpf('-inf')
        e0 = mp.exp(d*l0) if a0 > 0 else mp.mpf(0)
        e1 = mp.exp(d*l1) if a1 > 0 else mp.mpf(0)
        F = a0*e0 + a1*e1
        Phi = (a0*l0*e0 - a1*l1*e1)/F   # a_j' = +-w|<n_j,dir>| = +1*t-deriv... a0'=+1 in a-units? see below
        # exact: dF/dt = p*(a0' e0 + a1' e1), a0' = +c0, a1' = -c1 with c_j = |<v_j, V1-V0>| = 1 here
        Phi = (e0*l0*0 + 0)  # placeholder replaced below
        return None
    return None
'''
# The Phi line needs care: a_j(t) are affine with slopes a0' = |<v0,V1-V0>| = +c, a1' = -c; but also
# a_j*l_j' = a_j'. Rather than hand-derive, evaluate Phi = (1/p) d/dt log F numerically via the
# closed form Phi = (a0' e^{d l0} + a1' e^{d l1})/F with a0' = +1, a1' = -1 (since |<v_j,y>| = t or 1-t).
# Rewrite the script cleanly:
code = r'''# step_07.py: INDEPENDENT high-precision P^(k)(1), k=2..5, sheared parallelogram +-(2,0), +-(1,2).
# g_t(d) = |x_{1+d}'(t)| analytic on |d|<1; Cauchy radius 0.5, N nodes; tanh-sinh in t.
import mpmath as mp, time, json, sys
mp.mp.dps = int(sys.argv[1]) if len(sys.argv) > 1 else 70

def edge_geom(V0, V1):
    L = mp.sqrt((V1[0]-V0[0])**2 + (V1[1]-V0[1])**2)
    tau = ((V1[0]-V0[0])/L, (V1[1]-V0[1])/L)
    h = abs(V0[0]*V1[1] - V0[1]*V1[0])/L
    s0 = V0[0]*tau[0] + V0[1]*tau[1]
    return L, tau, h, s0

N = 220; R = mp.mpf('0.5')
om = [mp.exp(2*mp.pi*mp.j*k/N) for k in range(N)]
def Ek_of_t(L, tau, h, s0, t, kmax=5):
    """[E_k(t)] = [g^{(k)}(0)/k!] for the edge with |<v0,y>|=t, |<v1,y>|=1-t (power convention)"""
    l0 = mp.log(t); l1 = mp.log(1-t)
    sig = s0 + L*t
    vals = []
    for k in range(N):
        d = R*om[k]; p = 1 + d
        e0 = mp.exp(d*l0); e1 = mp.exp(d*l1)
        F = t*e0 + (1-t)*e1
        Phi = (e0 - e1)/F           # a0' = +1, a1' = -1 (in t-derivative of a_j)
        w2 = (L - sig*Phi)**2 + h**2*Phi**2
        vals.append(mp.sqrt(w2)*mp.exp(-mp.log(F)/p))
    out = []
    for k in range(kmax+1):
        s = mp.fsum(vals[j]*om[j]**(-k) for j in range(N))
        out.append(mp.re(s/(N*R**k)))
    return out

def edge_sum(V0, V1, kmax=5):
    L, tau, h, s0 = edge_geom(V0, V1)
    def integrand(t, k):
        return Ek_of_t(L, tau, h, s0, t, kmax)[k]
    res = []
    for k in range(kmax+1):
        f = lambda t: integrand(t, k)
        val = mp.quad(f, [0, mp.mpf('0.5')]) + mp.quad(f, [mp.mpf('0.5'), 1])
        res.append(val)
    return res

t0 = time.time()
VA0, VA1 = (mp.mpf(2), mp.mpf(0)), (mp.mpf(1), mp.mpf(2))
VB0, VB1 = (mp.mpf(1), mp.mpf(2)), (mp.mpf(-2), mp.mpf(0))
SA = edge_sum(VA0, VA1); print("edge A done %.1f s" % (time.time()-t0))
SB = edge_sum(VB0, VB1); print("edge B done %.1f s" % (time.time()-t0))
Pk_num = [mp.factorial(k)*(2*SA[k] + 2*SB[k]) for k in range(6)]
z2 = mp.zeta(2); z3 = mp.zeta(3); z5 = mp.zeta(5)
r5 = mp.sqrt(5); r13 = mp.sqrt(13)
exact = {
 2: z2*(64*r13/169 + 64*r5/25),
 3: z2*(-384*r5/25 - 384*r13/169) + z3*(192*r13/169 + 192*r5/25),
 4: z2**2*(-29952*r5/625 - 315648*r13/142805) + z2*(2304*r13/169 + 2304*r5/25) + z3*(-2304*r5/25 - 2304*r13/169),
 5: z2**2*(1262592*r13/28561 + 119808*r5/125) + z2*z3*(-10752*r5/25 - 821760*r13/28561)
    + z2*(-3072*r5/5 - 15360*r13/169) + z3*(23040*r13/169 + 4608*r5/5) + z5*(-6912*r5/25 - 172800*r13/28561)}
print("dps", mp.mp.dps, "wall %.1f s" % (time.time()-t0))
for k in [2, 3, 4, 5]:
    print(f"P^({k})(1) numeric = {mp.nstr(Pk_num[k], 30)}")
    print(f"   exact          = {mp.nstr(exact[k], 30)}")
    print(f"   |diff|         = {mp.nstr(abs(Pk_num[k]-exact[k]), 3)}")
json.dump({str(k): mp.nstr(Pk_num[k], 60) for k in range(2, 6)}, open('steps/step_07_num.json', 'w'))
'''
open('steps/step_07.py','w').write(code)
import subprocess, time
t0=time.time()
r=subprocess.run(['<project>/.venv/bin/python','steps/step_07.py','70'],capture_output=True,text=True,timeout=3000)
print("wall: %.1f s"%(time.time()-t0))
print(r.stdout); print("ERR:", r.stderr[-2000:])
