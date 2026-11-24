import mpmath as mp
from flint import fmpz_mat, fmpz
import time

def find_relation(vals, dps, scale_exp_margin=25, delta=0.99):
    """vals: list of mpf at current mp.mp.dps >= dps. Look for integer relation sum a_i vals[i] = 0.
    Returns (coeffs, residual) with residual computed at full current precision."""
    mp.mp.dps = dps
    C = mp.mpf(10) ** (dps - scale_exp_margin)
    n = len(vals)
    M = fmpz_mat(n, n + 1)
    for i, v in enumerate(vals):
        M[i, i] = fmpz(1)
        cv = mp.nint(C * v)
        M[i, n] = fmpz(int(cv))
    t0 = time.time()
    R = M.lll(delta=delta, eta=0.51)
    t1 = time.time()
    # shortest row = candidate
    rows = []
    for r in range(n):
        row = [R[r, c] for c in range(n + 1)]
        rows.append(row)
    return rows, t1 - t0

# ---- test on W3 with the G1 11-element weight-3 basis ----
mp.mp.dps = 1050
W3 = mp.mpf(open("hex_W3.txt").read())
l2, l3 = mp.log(2), mp.log(3)
z2 = mp.pi**2/6
def Li(s, x): return mp.polylog(s, x)
basis3_names = ["zeta3", "Li3(1/3)", "Li3(-1/2)", "pi2*log2", "pi2*log3",
                "Li2(1/3)*log2", "Li2(1/3)*log3", "log2^3", "log2^2*log3", "log2*log3^2", "log3^3"]
basis3 = [mp.zeta(3), Li(3, mp.mpf(1)/3), Li(3, -mp.mpf(1)/2), mp.pi**2*l2, mp.pi**2*l3,
          Li(2, mp.mpf(1)/3)*l2, Li(2, mp.mpf(1)/3)*l3, l2**3, l2**2*l3, l2*l3**2, l3**3]
rows, dt = find_relation([W3] + basis3, 600)
print(f"LLL time: {dt:.1f}s")
for r in rows[:3]:
    coeffs = r[:-1]
    resid = sum(int(c)*v for c, v in zip(coeffs, [W3] + basis3))
    height = max(abs(int(c)) for c in coeffs)
    print(f"height {height}, residual {mp.nstr(abs(resid), 3)}")
    if abs(resid) < mp.mpf(10)**(-500):
        a0 = coeffs[0]
        terms = [(basis3_names[i-1], int(-coeffs[i])/int(a0)) for i in range(1, len(coeffs)) if coeffs[i]]
        print("W3 =", " + ".join(f"({mp.nstr(mp.mpf(t[1]), 6)})*{t[0]}" for t in terms))
