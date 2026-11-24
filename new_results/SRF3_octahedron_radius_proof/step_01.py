
"""
SRF3 step_01: exact representation of S(p) at the octahedron end, numerical check.

CLAIM (to be proved in proof.md):
  S(p) = 8 (1-tau)^2 T(tau),   tau = 1 - 1/p,
  T(tau) = ∫_Δ (t1 t2 t3)^(-tau) sqrt(t1^(2tau)+t2^(2tau)+t3^(2tau)) dt1 dt2,
  Δ = {t1,t2 >= 0, t1+t2 <= 1}, t3 = 1-t1-t2.

Equivalent forms:
 (square)  t1=(1-u)s, t2=(1-u)(1-s), t3=u:
  T(tau) = ∫_0^1∫_0^1 [u(1-u)^2 s(1-s)]^(-tau)
           sqrt((1-u)^(2tau)(s^(2tau)+(1-s)^(2tau)) + u^(2tau)) (1-u) ds du.
 (Laplace) u = e^{-V}:
  T(tau) = ∫_0^∞ e^{-V} Φ(tau,V) dV,
  Φ(tau,V) = e^{tau V}(1-e^{-V})^{1-2tau} ∫_0^1 [s(1-s)]^(-tau)
             sqrt((1-e^{-V})^(2tau)(s^(2tau)+(1-s)^(2tau)) + e^{-2 tau V}) ds.

Check: S(1+d) against H1's independent direct surface quadrature (40 dps),
targeting >= 30 digits of agreement.
"""
import mpmath as mp, json, os
mp.mp.dps = 60
base = "<project>"
direct = json.load(open(os.path.join(base,"research/10_frontier/H1_octahedron_end/direct_S_dps40.json")))

def T_square(tau):
    tau = mp.mpf(tau)
    f = lambda s, u: ((u*(1-u)**2*s*(1-s))**(-tau)
                      * mp.sqrt((1-u)**(2*tau)*(s**(2*tau)+(1-s)**(2*tau)) + u**(2*tau))
                      * (1-u))
    return mp.quad(f, [0,1],[0,1])

def S_from_repr(p):
    tau = 1 - 1/mp.mpf(p)
    return 8*(1-tau)**2 * T_square(tau)

def T_laplace(tau):
    tau = mp.mpf(tau)
    def Phi(V):
        ev = mp.exp(-V)
        inner = lambda s: ((s*(1-s))**(-tau)
                           * mp.sqrt((1-ev)**(2*tau)*(s**(2*tau)+(1-s)**(2*tau)) + mp.exp(-2*tau*V)))
        return mp.exp(tau*V)*(1-ev)**(1-2*tau)*mp.quad(inner,[0,1])
    return mp.quad(lambda V: mp.exp(-V)*Phi(V), [0, mp.inf])

rows = []
for dstr in ["-0.1","0.1","-0.03","0.03","-0.01","0.01"]:
    d = mp.mpf(dstr)
    p = 1 + d
    ref = mp.mpf(direct[dstr]["N1"])
    val = S_from_repr(p)
    err = abs(val-ref)
    digits = -mp.log10(err/abs(ref))
    rows.append((dstr, val, ref, err, digits))
    print(f"d={dstr:>6}: S={mp.nstr(val,40)}")
    print(f"          ref={mp.nstr(ref,40)}  |err|={mp.nstr(err,3)}  digits={mp.nstr(digits,4)}")

# Laplace form cross-check at two d values
for dstr in ["-0.03","0.03"]:
    d = mp.mpf(dstr); p = 1+d; tau = 1-1/p
    v1 = 8*(1-tau)**2*T_laplace(tau)
    v2 = S_from_repr(p)
    print(f"Laplace-form check d={dstr}: agree to {mp.nstr(-mp.log10(abs(v1-v2)/abs(v2)),4)} digits")

worst = min(float(r[4]) for r in rows)
print("WORST digits vs direct quadrature:", worst)
assert worst > 30
print("step_01 OK")
