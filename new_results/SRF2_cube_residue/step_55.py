# verify_srf2.py -- SRF2: pole of the vertex term V at eps=-1 and its residue
# R = -1 - 2*pi^{3/2}/Gamma(1/4)^2.  Runs in ~2 min.  Work dir: research/08_theory/SRF2_cube_residue.
import sys, json, time
sys.path.insert(0, '<project>/research/07_push/S3_superellipsoid/verify')
sys.path.insert(0, '<project>/research/08_theory/SRF2_cube_residue')
import mpmath as mp
import flint
from flint import arb, ctx
from vlib import ts_nodes, V_series, tompf
from srf2_lib import expm1_s, expm1m, u_minus_log1p, Phi, tomp   # step_24, step_26 (stable integrands)

t0 = time.time(); ctx.dps = 60; mp.mp.dps = 60
G4 = mp.gamma(mp.mpf('0.25'))
R  = -1 - 2*mp.pi**mp.mpf('1.5')/G4**2
I  = 1 - 2*mp.pi**mp.mpf('1.5')/G4**2

# PART 1 (step_38): closed form of I and R, every proof step numerically
I_dir = mp.quad(lambda a: a*a/(mp.sqrt(1+a**4)+1), [0,1])
K     = mp.quad(lambda a: a*a/mp.sqrt(1+a**4), [0,1])
Lp    = mp.quad(lambda t: t*t/mp.sqrt(1-t**4), [0,1])
assert abs(I_dir - I) < mp.mpf('1e-50')
assert abs(1-mp.sqrt(2)+2*K - I) < mp.mpf('1e-50')                      # integration by parts
assert abs(Lp - mp.sqrt(2)*mp.pi**mp.mpf('1.5')/G4**2) < mp.mpf('1e-50') # (1/4)B(3/4,1/2)
assert abs(Lp - 2*mp.sqrt(2)*mp.quad(lambda a: a*a/(1+a**4)**mp.mpf('1.5'), [0,1])) < mp.mpf('1e-50')  # substitution
assert abs(mp.quad(lambda a: a*a/(1+a**4)**mp.mpf('1.5'), [0,1]) - (1/(2*mp.sqrt(2))-K/2)) < mp.mpf('1e-50')
print("PART 1 ok: I =", mp.nstr(I_dir, 30), " R =", mp.nstr(R, 30))

# PART 2 (step_36/37): residue via boundary-layer-aware split integration of J
def J_small(eps, a0=mp.mpf('0.5'), Lam=mp.mpf('160')):
    s = -2-2*eps; Lmax = 200
    g1 = [mp.mpf(1)]
    for l in range(1, Lmax+1): g1.append(g1[-1]*(s-l+1)/l)
    bh = [mp.mpf(1)]
    for j in range(1, 60): bh.append(bh[-1]*(mp.mpf('0.5')-j+1)/j)
    two2e = 2-2*eps; tot = mp.mpf(0)
    for j in range(0, 60):
        for l in range(0, Lmax+1):
            if l == 0 and j == 0: continue
            pw = l + two2e*j
            if pw > Lam: break
            c = g1[l]*bh[j]
            if c == 0: continue
            mu = 2*eps + pw
            tot += c * a0**(mu+1)/(mu+1)
    return tot
def J_big(eps, h=arb(1)/64, T=5.0):
    a0 = arb('0.5'); s = arb(0)
    for (x, xm, w, par) in ts_nodes(h, T):
        a = a0 + (1-a0)*x; la = a.log(); l1a = (1+a).log()
        A2 = a*a*(arb(-2*eps)*la).exp(); S0 = (1+A2).sqrt()
        f = S0*expm1_s(arb(-(2+2*eps))*l1a) + A2/(S0+1)
        s += (arb(2*eps)*la).exp()*f*w*(1-a0)
    return s
def J_rel(e):
    return 1/(2*e+1) + J_small(e) + tomp(J_big(arb(mp.nstr(e, 30))))
assert abs(J_rel(mp.mpf(-1)) - (-1 + I)) < mp.mpf('1e-24')              # face value = -1 + I
pts = []
for eta in ['0.05', '0.02', '0.01', '0.005']:
    e = -1 + mp.mpf(eta); q = 1+2*e
    pts.append((mp.mpf(eta), -q*J_rel(e)))
A = mp.matrix([[x**c for c in range(3)] for x, y in pts]); b = mp.matrix([y for x, y in pts])
R_est = mp.lu_solve(A, b)[0]
print("PART 2 ok: -qJ extrapolated R =", mp.nstr(R_est, 15), " residual:", mp.nstr(abs(R_est-R), 3))
assert abs(R_est - R) < mp.mpf('2e-6')

# PART 3 (step_47): rho_{-2} = 3(-2 + I2')
I2 = mp.quad(lambda a: (1+a)*a**3/(mp.sqrt(1+a**6)+1), [0, mp.mpf('0.3'), 1])
rho2 = 3*(-2+I2)
print("PART 3 ok: rho_{-2} =", mp.nstr(rho2, 20))

# PART 4 (step_49/50): V_k to order 40 (h=1/32) vs stored high-precision values + delta_k behaviour
flint.ctx.cap = 64; ctx.dps = 100
Vf, _ = V_series(0, 41, h=arb(1)/32, T=4.6)
V = [tompf(x, 110) for x in Vf]
stored = [mp.mpf(x) for x in json.load(open('EV_coeffs_K40.json'))['V']]
worst = max(abs(V[k]-stored[k])/abs(stored[k]) for k in range(1, 41))
print("PART 4 ok: worst rel diff vs K40 store =", mp.nstr(worst, 3))
d40 = (-1)**40*V[40] - R
d10 = ((-1)**10*V[10] - R)*2**11
print("  (-1)^40 V_40 - R =", mp.nstr(d40, 5), " (accumulation);  2^11 delta_10 =", mp.nstr(d10, 8), " vs rho_-2 =", mp.nstr(rho2, 8))
assert abs(d40) < mp.mpf('1e-7') and abs(d10 - rho2)/abs(rho2) < mp.mpf('0.06')

# PART 5 (step_45/46): E regular at -1
def E_int(eps, h=arb(1)/32, T=4.6):
    q = 1+2*eps; s = arb(0)
    for (a, am, wa, pa) in ts_nodes(h, T):
        la = a.log(); l1a = (1+a).log()
        A2 = a*a*(arb(-2*eps)*la).exp(); S0 = (1+A2).sqrt()
        f = Phi(a, q)*S0 + (1-q*a)*A2/(S0+1)
        s += (arb(2*eps)*la).exp()/a*f*wa
    return -arb(1+2*eps)/(2*eps+1) + s
E_m1 = tomp(E_int(arb(-1)))
print("PART 5 ok: E(-1) (subtracted-integral value) =", mp.nstr(E_m1, 20), " (finite: no pole)")

# PART 6 (step_50): s_k = 48(E_{k-1}+V_{k-2}) vs SRF1
E3 = [mp.mpf(x) for x in json.load(open('EV_coeffs_K40.json'))['E3']]
s = [mp.mpf(24)] + [48*(E3[k-1] + (stored[k-2] if k >= 2 else 0)) for k in range(1, 42)]
print("PART 6 ok: s5 =", mp.nstr(s[5], 12), " s6 =", mp.nstr(s[6], 12),
      " 48R(-1)^5 =", mp.nstr(48*R*(-1)**5, 12))
assert abs(s[5] - mp.mpf('90.7414583339')) < mp.mpf('1e-8') and abs(s[6] - mp.mpf('-93.2703787987')) < mp.mpf('1e-8')
print(f"ALL CHECKS PASSED in {time.time()-t0:.0f}s")
