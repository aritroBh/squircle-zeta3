# step_09.py: validate the FULL E5 (tau,gam symbols) on a NON-fully-local edge:
# regular decagon (m=5) edge, where the 3 non-local kinks have nonconstant logs.
# Compares symbolic E5(t) (tau_k = sum a_i l_i^k, gam_k = sum a_i' l_i^k evaluated exactly)
# against the per-t Cauchy d-series of |x_p'(t)| at several t.  Also evaluates the local
# approximation (l_i -> mu_i, frozen) to exhibit the non-local defect.
import mpmath as mp, sympy as sp, pickle
mp.mp.dps = 60
N = 200
om = [mp.exp(2*mp.pi*mp.j*k/N) for k in range(N)]

m = 5
# decagon K_1: normals n_j = (cos th_j, sin th_j), th_j = pi/2 - j pi/m, weights 1, fixed-weight convention (lam=0)
th = [mp.pi/2 - j*mp.pi/m for j in range(m)]
h_ = mp.pi/(2*m); R_ = mp.tan(h_)
V0 = (R_, mp.mpf(0)); V1 = (R_*mp.cos(2*h_), R_*mp.sin(2*h_))
L = mp.sqrt((V1[0]-V0[0])**2+(V1[1]-V0[1])**2)
tau_ = ((V1[0]-V0[0])/L, (V1[1]-V0[1])/L)
he = abs(V0[0]*V1[1]-V0[1]*V1[0])/L
s0 = V0[0]*tau_[0]+V0[1]*tau_[1]
# kinks: j0 = 0 (n=(0,1) through V0), j1 = m-1 = 4 (through V1); non-local j = 1,2,3
def aj(j, t):
    y0 = (1-t)*V0[0]+t*V1[0]; y1 = (1-t)*V0[1]+t*V1[1]
    return abs(mp.cos(th[j])*y0 + mp.sin(th[j])*y1)
def cauchy_E5(t):
    sig = s0 + L*t
    lj = [mp.log(aj(j, t)) for j in range(m)]
    ajv = [aj(j, t) for j in range(m)]
    # slopes a_j' = w |<n_j, V1-V0>| * sign adjustment: a_j(t) affine; a_j l_j' = a_j'
    # compute a_j' numerically: a_j' = d/dt a_j(t)
    ap = []
    for j in range(m):
        c0 = mp.cos(th[j])*V0[0]+mp.sin(th[j])*V0[1]
        c1 = mp.cos(th[j])*V1[0]+mp.sin(th[j])*V1[1]
        sg = 1 if c0+1e-30 > 0 else -1
        ap.append(sg*(c1-c0))
    vals = []
    aF = aw = None
    for jj in range(N):
        d = 0.35*om[jj]; p = 1+d
        F = mp.fsum(ajv[j]*mp.exp(d*lj[j]) for j in range(m))
        Phi = mp.fsum(ap[j]*mp.exp(d*lj[j]) for j in range(m))/F
        w2 = (L - sig*Phi)**2 + he**2*Phi**2
        a1 = mp.arg(F); a2 = mp.arg(w2)
        if aF is not None:
            a1 += 2*mp.pi*int(mp.nint((aF - a1)/(2*mp.pi)))
            a2 += 2*mp.pi*int(mp.nint((aw - a2)/(2*mp.pi)))
        aF, aw = a1, a2
        logF = mp.log(abs(F)) + 1j*a1; logw = mp.log(abs(w2)) + 1j*a2
        vals.append(mp.exp(logw/2 - logF/p))
    return [mp.re(mp.fsum(vals[jj]*om[jj]**(-k) for jj in range(N))/(N*mp.mpf('0.35')**k)) for k in range(6)]

# symbolic full E5
t_, al, be, Ls, hs, s0s = sp.symbols('t alpha beta L h s0', real=True)
l0s, l1s = sp.symbols('l0 l1')
taus = sp.symbols('tau1:6'); gams = sp.symbols('gam1:6')
D = pickle.load(open('steps/step_02_full.pkl','rb'))
loc = {'t': t_, 'alpha': al, 'beta': be, 'L': Ls, 'h': hs, 's0': s0s, 'l0': l0s, 'l1': l1s}
for i, ts_ in enumerate(taus): loc[f'tau{i+1}'] = ts_
for i, gs_ in enumerate(gams): loc[f'gam{i+1}'] = gs_
E5full = sp.sympify(D['E5'], locals=loc)
# E5 reduced (local, tau=gam=0)
Rr = pickle.load(open('steps/step_01_reduced.pkl','rb'))
E5red = sp.sympify(Rr['E'][5], locals=loc)

alpha_ = aj(0, 1); beta_ = aj(m-1, 0)
def symE5(t, mode):
    lj = [mp.log(aj(j, t)) for j in range(m)]
    ap = []
    for j in range(m):
        c0 = mp.cos(th[j])*V0[0]+mp.sin(th[j])*V0[1]
        c1 = mp.cos(th[j])*V1[0]+mp.sin(th[j])*V1[1]
        sg = 1 if c0+1e-30 > 0 else -1
        ap.append(sg*(c1-c0))
    sub = {al: alpha_, be: beta_, Ls: L, hs: he, s0s: s0, t_: t,
           l0s: lj[0], l1s: lj[m-1]}
    if mode == 'full':
        for k in range(1, 6):
            sub[taus[k-1]] = mp.fsum(aj(j, t)*lj[j]**k for j in range(1, m-1))
            sub[gams[k-1]] = mp.fsum(ap[j]*lj[j]**k for j in range(1, m-1))
    expr = E5full if mode == 'full' else E5red
    return mp.mpf(sp.N(expr.subs({k_: mp.mpf(str(v_)) for k_, v_ in sub.items()}), 60))

for tv in ['0.15', '0.35', '0.5', '0.75', '0.9']:
    t = mp.mpf(tv)
    Ec = cauchy_E5(t)[5]
    Ef = symE5(t, 'full'); El = symE5(t, 'local')
    print(f"t={tv}: Cauchy E5={mp.nstr(Ec,16)}")
    print(f"   full sym E5={mp.nstr(Ef,16)}  |diff|={mp.nstr(abs(Ec-Ef),3)}   local-reduced E5={mp.nstr(El,12)}  |local defect|={mp.nstr(abs(Ef-El),3)}")
