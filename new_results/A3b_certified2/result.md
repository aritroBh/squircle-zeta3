(no FINAL answer; see transcript.json)

Last reply:
FINAL:

# Task A3b-2 — Squircle W4 certified-numerics run (flag fix applied; prec-400 run complete; prec-1100 run **not executed** — session ended)

**Bottom line: PARTIAL / NUMERICAL for the 10⁻³⁰⁰ goal.** The rigor fix (forwarding `analytic` to every branch-cut call in all six integrands) was implemented, audited, and verified behaviorally; two additional latent crashes in the draft were found and fixed; the fixed script `w4_certified.py` ran end-to-end at PREC=400 / EPS=2⁻⁵⁰⁰ and **printed certified bounds `|c0| < 10^-112 : True` and `|c1| < 10^-111 : True`** (certified ball-comparison semantics), with midpoints identical to the previous session's unflagged run. The production PREC=1100 / EPS=2⁻¹³⁰⁰ run was **not executed**: this turn was declared the last before I could launch it. Per the rules I therefore do **not** claim the `10^-300` theorem. Everything below was printed by code that ran here, except where marked as quoted baseline from `../A3b_certified/result.md`.

## 1. Rigor fix applied (required item)

python-flint's own docstring (printed here) states the flag "*must* be checked for functions … with branch cuts such as √x and log(x), even if the integration path does not touch any non-analytic points", and I verified the mechanism in-session:

```
zball = acb(arb(-1, 0.5))            # straddles the negative real axis
zball.log()                -> [+/- 0.694] + [3.1415926535897932384626433832795029 +/- 1.68e-35]j
zball.log(analytic=True)   -> nan + nanj        (non-finite => acb.integral subdivides)
```

**Branch-cut audit of the six integrands:** the only branch-cut functions are `log` — one call in each of `fU_mid, fU_left, fU_right`, two in each of `fV_mid, fV_left, fV_right` (8 total). All are now called as `.log(analytic=analytic)`. There is **no `sqrt` and no `pow`** anywhere in the integrands (all powers are products `l*l*l`, `l*l*m`); `exp` is entire; `+, -, *, /` are meromorphic, and per the printed python-flint documentation poles automatically yield non-finite balls, so division needs no flag.

## 2. Two further draft bugs found and fixed (the draft could not run at all)

1. **`P(C,k)` crashed on its `j=0` term**: `comb(k,0)*mul(*[w+C]*0)*factorial(k)` calls `mul()` with zero arguments → reproduced here: `IndexError: tuple index out of range`. Replaced by an explicit loop using only `*` (no `arb_pow`).
2. **`certified_decimal_exp` used `int(arb)`** → reproduced here: `TypeError: int() argument must be ... not 'flint.types.arb.arb'`. Replaced by a pure certified-comparison loop: increment `k` while `bool(mag(ci) < arb(10)**(-(k+1)))` (arb `<` returns True only when certified), so the printed `k` is backed by a printed `True` comparison.

Also fixed (cosmetic/robustness, not affecting values): integral radii now printed from `IU.real.rad()`/`IU.imag.rad()` (the draft's `mag(IU-U+U-IU)` printed nonsense); tail balls added via `rU.abs_upper()`; the mpmath cross-check's `num()` made `Fraction`-safe.

## 3. Evidence — fixed `w4_certified.py`, run as `python w4_certified.py 400 500` (returncode 0, wall 2.2 s; full stdout saved to `out_400.txt`)

**`acb.integral time: 0.4s`** (all six pieces, with the analytic flag active).

| quantity | printed value (midpoint) | printed radius | certified digits |
|---|---|---|---|
| `IU_left` | `-1.012713010069575783392292461494559209433321481941954588935871354802448878286946547937038738642345425319316007906134353` | `3.65e-118` | 117.4 |
| `IU_mid` | `-0.003446309720200776041845442335850591249926660865127092298813263024988615340010336037455911665322176438204535244263511771` | `5.24e-121` | — |
| `IU_right` | `-0.0002707836536969512028882312213522812820586996282385570594162752930346215263679675802968661421986253172197284707720592559` | `8.48e-122` | — |
| `IU` (= U integral) | `-1.016430103443473510637026135051762081965306842435320238294100893120472115153324851554791516449866227074740271621169924` | re rad `1.594e-118`, im 0 | 117.8 |
| `IV` (= V integral) | `1.281429630921344022065641268801917196814522298172349739197686393033896817414004840408672212708216615058232091817961969` | re rad `2.012e-118`, im 0 | 117.7 |
| tail `T_U_b` | `3.6666195543427906013746548775665725369731873731824264678193626901512174651923116835804506377210482206070173710753089522e-144` | `7.00e-263` | — |
| tail `T_V_b` | `4.5372102140872220585855570925461599860827380106898171338226029004525593707778911579530581844761187492076743726399121062e-144` | `8.05e-263` | — |
| tail `T_U_a` | `0` | `±1.93e-631` | — |
| tail `T_V_a` | `0` | `±5.54e-390` | — |
| tail rad `rU` / `rV` | — | `3.6666195…e-144` / `4.5372102…e-144` | — |
| `U` | `-1.016430103443473510637026135051762081965306842435320238294100893120472115153324851554791516449866227074740271621169924` | `±3.41e-118` | 117.5 |
| `V` | `1.281429630921344022065641268801917196814522298172349739197686393033896817414004840408672212708216615058232091817961969` | `±2.88e-118` | 117.5 |
| `R0` | `-1.01643010344347351063702613505176208196530684243532023829410089312047211515332485155479151644986622707474027162` | re rad `2.660e-111` | 110.6 |
| `R1` | `1.2814296309213440220656412688019171968145222981723497391976863930338968174140048404086722127082166150582320918` | re rad `1.695e-110` | 109.8 |
| **`c0` (residual)** | `0` | **`±2.74e-113`** | — |
| **`c1` (residual)** | `0` | **`±1.75e-112`** | — |

(Digit counts are `−log10(rad/|mid|)` computed by hand from the printed radii — derived arithmetic on printed balls, not itself printed.)

**CERTIFIED lines, copied verbatim from the run output:**

```
c0 = [+/- 2.74e-113]
  |c0| bound ball: [+/- 2.74e-113]
  CERTIFIED |c0| < 10^-112 : True
  CERTIFIED |c0| < 1e-300 : False
c1 = [+/- 1.75e-112]
  |c1| bound ball: [+/- 1.75e-112]
  CERTIFIED |c1| < 10^-111 : True
  CERTIFIED |c1| < 1e-300 : False
NUMERICAL mp  U = -1.0164301034434735106370261350517620819653068424353
NUMERICAL mp  V = 1.2814296309213440220656412688019171968145222981723
```

**Exactly which bound was certified (k as printed): `|c0| < 10^-112` and `|c1| < 10^-111`, at PREC=400 / EPS=2⁻⁵⁰⁰** — consistent with the previous session's hand-assembled estimates `3.9e-113`, `3.6e-112`. The `1e-300` comparisons printed `False`, as they must at 400 bits (the script only asserts the 1e-300 target when `PREC ≥ 1100`).

**Flag fix did not change midpoints (prec-400 comparison):** the fixed script's `IU, IV, R0, R1` midpoints printed above are character-for-character the same strings as the previous session's unflagged prec-400 midpoints quoted in `../A3b_certified/result.md` §4 (e.g. `IU = -1.016430103443473510637026135051762081965306842435320238294100893120472115153324851554791516449866227074740271621169924`; `R1 = 1.2814296309213440220656412688019171968145222981723497391976863930338968174140048404086722127082166150582320918`; tails `T_U_b = 3.6666195e-144`, `T_V_b = 4.5372102e-144`, factors `L1 = ±1.42e-120`, `H_U_a = ±2.22e-121` likewise). Radii moved only at display level (`±3.41e-118` vs quoted `3.39e-118`), consistent with slightly different adaptive subdivision under the analytic flag. (My programmatic string-comparison cell had a dict-key typo and printed `common prefix 0`; the comparison stated here is by direct inspection of the two printed outputs — the run output above is in `out_400.txt`, the baseline in `result.md`.)

## 4. What remains (not done — session ended)

Launch `python w4_certified.py` (defaults `PREC=1100`, `EPS_EXP=1300`, exactly the draft's production settings). At prec 400 the six integrals took 0.4 s, so the 1100-bit run is expected to be short — but this is an expectation, **not** a printed result, and the `10^-300` bound is **not claimed**. The script will print `CERTIFIED |c0| < 10^-k` / `|c1| < 10^-k` with the k it certifies and will assert `k ≥ 300`.

## 5. Final `w4_certified.py` (as run above)

```python
"""Certified bounds for |c0|, |c1| (squircle W4 residual constants).
Method: rigorous ball arithmetic (python-flint arb/acb).
  U, V: acb.integral on three holomorphic pieces (subs v=b+e^s, v=a-e^s),
        plus endpoint tails bounded by exact formulas for int(log(1/s)+C)^k ds.
  R0, R1: all atoms of conjecture_terms() as balls (z.polylog(s) = Li_s(z)).
RIGOR (analytic flag): acb.integral calls f(x, analytic); with analytic=True
  the integrand must be non-finite whenever it is not holomorphic on the ball,
  else the Gauss-Legendre error bound is not rigorous. Every branch-cut
  function in the six integrands is log; each is called as .log(analytic=analytic).
  (exp is entire; division is meromorphic - poles give non-finite automatically.)
Run: python w4_certified.py [PREC=1100] [EPS_EXP=PREC+100]
"""
import sys, time
sys.dont_write_bytecode = True
from math import comb, factorial
from flint import arb, acb, ctx
from w4_symbol import conjecture_terms

PREC = int(sys.argv[1]) if len(sys.argv) > 1 else 1100
EPS_EXP = int(sys.argv[2]) if len(sys.argv) > 2 else PREC + 100
ctx.prec = PREC
TARGET = arb(10) ** (-300)
print('PREC =', PREC, ' EPS_EXP =', EPS_EXP, flush=True)

def mag(z):                          # arb ball containing |z'| for z' in z
    return abs(z.real) + abs(z.imag)

eps   = arb(2) ** (-EPS_EXP)
sqrt2 = acb(2).sqrt()
a, b  = 1 + sqrt2, sqrt2 - 1
log2  = acb(2).log()
pi    = acb.pi()
c1_6, c1_3   = acb(arb(1)/arb(6)), acb(arb(1)/arb(3))
c7_48, c1_16 = acb(arb(7)/arb(48)), acb(arb(1)/arb(16))
c7_48r, c1_16r = arb(7)/arb(48), arb(1)/arb(16)

# ---- integrands: analytic flag forwarded to EVERY branch-cut call (log) ----
def fU_mid(v, analytic=False):
    l = ((v - b)*(v + a)/(4*v)).log(analytic=analytic)
    return l*l*l * (c1_6/v - c1_3/(v + a))
def fV_mid(v, analytic=False):
    l = ((v - b)*(v + a)/(4*v)).log(analytic=analytic)
    m = ((a - v)*(v + b)/(4*v)).log(analytic=analytic)
    return -c7_48*l*l*l/v + c1_16*l*l*m/v
def fU_left(s, analytic=False):      # v = b + e^s in [b+eps, 1]
    e = s.exp(); v = b + e
    l = s + ((v + a)/(4*v)).log(analytic=analytic)
    return l*l*l * (c1_6/v - c1_3/(v + a)) * e
def fV_left(s, analytic=False):
    e = s.exp(); v = b + e
    l = s + ((v + a)/(4*v)).log(analytic=analytic)
    m = ((a - v)*(v + b)/(4*v)).log(analytic=analytic)
    return (-c7_48*l*l*l/v + c1_16*l*l*m/v) * e
def fU_right(s, analytic=False):     # v = a - e^s in [a-1, a-eps]
    e = s.exp(); v = a - e
    l = ((v - b)*(v + a)/(4*v)).log(analytic=analytic)
    return l*l*l * (c1_6/v - c1_3/(v + a)) * e
def fV_right(s, analytic=False):
    e = s.exp(); v = a - e
    l = ((v - b)*(v + a)/(4*v)).log(analytic=analytic)
    m = s + ((v + b)/(4*v)).log(analytic=analytic)
    return (-c7_48*l*l*l/v + c1_16*l*l*m/v) * e
# branch-cut audit: integrands contain only log (flagged), exp (entire),
# +,-,*,/ (meromorphic). No sqrt, no pow anywhere (powers are products).

s_lo, s_ml, s_mr = -EPS_EXP*log2, (1 - b).log(), acb(0)
t0 = time.time()
IU_left  = acb.integral(fU_left,  s_lo, s_ml)
IU_mid   = acb.integral(fU_mid,   acb(1), a - 1)
IU_right = acb.integral(fU_right, s_lo, s_mr)
IV_left  = acb.integral(fV_left,  s_lo, s_ml)
IV_mid   = acb.integral(fV_mid,   acb(1), a - 1)
IV_right = acb.integral(fV_right, s_lo, s_mr)
IU = IU_left + IU_mid + IU_right
IV = IV_left + IV_mid + IV_right
print('acb.integral time: %.1fs' % (time.time() - t0), flush=True)
print('IU_left  =', IU_left,  '\nIU_mid   =', IU_mid,  '\nIU_right =', IU_right)
print('IV_left  =', IV_left,  '\nIV_mid   =', IV_mid,  '\nIV_right =', IV_right)
for name, I in [('IU', IU), ('IV', IV)]:
    print('%s = %s\n  integral rad: re %.3e im %.3e' % (name, I, I.real.rad(), I.imag.rad()))
assert IU.is_finite() and IV.is_finite(), 'integral not finite'

# ---- endpoint tails (rigorous sup bounds over covering balls) ----
w = EPS_EXP * arb(2).log()                       # w = log(1/eps)
def P(C, k):                                     # (1/eps) * int_0^eps (log(1/s)+C)^k ds
    out = arb(0)                                 #   = sum_j C(k,j) (w+C)^j (k-j)!
    for j in range(k + 1):
        t = arb(comb(k, j)) * arb(factorial(k - j))
        for _ in range(j):                       # x*x*x, not x**j (arb_pow NaN on 0-balls)
            t = t * (w + C)
        out = out + t
    return out
vb = acb(b) + acb(arb(0, eps))                   # covers [b, b+eps]
C1, invb = mag(((vb + a)/(4*vb)).log()), mag(1/vb)
H_U_b = mag(c1_6/vb - c1_3/(vb + a))
M1    = mag(((a - vb)*(vb + b)/(4*vb)).log())
T_U_b = H_U_b * eps * P(C1, 3)
T_V_b = c7_48r*invb*eps*P(C1, 3) + c1_16r*invb*M1*eps*P(C1, 2)
va = acb(a) + acb(arb(0, eps))                   # covers [a-eps, a]
C2, inva = mag(((va + b)/(4*va)).log()), mag(1/va)
L1    = mag(((va - b)*(va + a)/(4*va)).log())
H_U_a = mag(c1_6/va - c1_3/(va + a))
S1    = (w + C2) + 1
T_U_a = H_U_a * eps * L1*L1*L1
T_V_a = c7_48r*inva*eps*L1*L1*L1 + c1_16r*inva*eps*L1*L1*S1
print('tail factors: C1=%s invb=%s H_U_b=%s M1=%s' % (C1, invb, H_U_b, M1))
print('tail factors: C2=%s inva=%s L1=%s H_U_a=%s S1=%s' % (C2, inva, L1, H_U_a, S1))
print('T_U_b =', T_U_b, '\nT_V_b =', T_V_b, '\nT_U_a =', T_U_a, '\nT_V_a =', T_V_a)
rU, rV = T_U_b + T_U_a, T_V_b + T_V_a
print('tail rad rU =', rU, '\ntail rad rV =', rV)
U = IU + acb(arb(0, rU.abs_upper()))
V = IV + acb(arb(0, rV.abs_upper()))
print('U =', U, '\n  tail rad %.3e' % rU.abs_upper())
print('V =', V, '\n  tail rad %.3e' % rV.abs_upper())

# ---- R0, R1 from the literal conjecture terms ----
def qs2(p, q):                                   # p + q*sqrt(2), Fractions
    v = acb(arb(p.numerator)/arb(p.denominator))
    return v + acb(arb(q.numerator)/arb(q.denominator))*sqrt2 if q else v
def atom_value(at):
    if at[0] == 'log':
        z = qs2(at[1].p, at[1].q)
        assert z.real > 0
        return acb(abs(z.real)).log()
    if at[0] == 'li':
        return qs2(at[2].p, at[2].q).polylog(acb(at[1]))
    if at == ('pi2',):
        return pi*pi
    assert at == ('z3',)
    return acb(3).zeta()
R0, R1 = acb(0), acb(0)
for c, atoms in conjecture_terms():
    val = acb(1)
    for at in atoms:
        val = val * atom_value(at)
    R0 += acb(arb(c.p.numerator)/arb(c.p.denominator)) * val
    R1 += acb(arb(c.q.numerator)/arb(c.q.denominator)) * val
print('R0 =', R0, '\n  rad: re %.3e im %.3e' % (R0.real.rad(), R0.imag.rad()))
print('R1 =', R1, '\n  rad: re %.3e im %.3e' % (R1.real.rad(), R1.imag.rad()))

# ---- residuals and the certified test ----
pi4 = pi*pi*pi*pi
c0, c1 = (U - R0)/pi4, (V - R1)/pi4
def certified_decimal_exp(ci):
    """Maximal integer k with mag(ci) < 10^-k certifiably (arb < semantics)."""
    m = mag(ci)
    k = 0
    while bool(m < arb(10) ** (-(k + 1))):
        k += 1
    return k
for name, ci in [('c0', c0), ('c1', c1)]:
    print(name, '=', ci)
    print('  |%s| bound ball:' % name, mag(ci))
    k = certified_decimal_exp(ci)
    okk = bool(mag(ci) < arb(10) ** (-k))
    ok3 = bool(mag(ci) < TARGET)                # certified ball comparison
    print('  CERTIFIED |%s| < 10^-%d :' % (name, k), okk)
    print('  CERTIFIED |%s| < 1e-300 :' % name, ok3)
    assert okk
    if PREC >= 1100:
        assert ok3, 'production run must beat 1e-300'

# ---- NUMERICAL cross-check (mpmath, not part of the certificate) ----
import mpmath as mp
mp.mp.dps = 100
def fmp(fr): return mp.mpf(fr.numerator) / mp.mpf(fr.denominator)
def num(x): return fmp(x.p) + mp.sqrt(2)*fmp(x.q)
am, bm = num(__import__('w4_symbol').A), num(__import__('w4_symbol').B)
fU = lambda v: mp.log((v-bm)*(v+am)/(4*v))**3*(1/(6*v)-1/(3*(v+am)))
fV = lambda v: (-7*mp.log((v-bm)*(v+am)/(4*v))**3/(48*v)
                + mp.log((v-bm)*(v+am)/(4*v))**2*mp.log((am-v)*(v+bm)/(4*v))/(16*v))
print('NUMERICAL mp  U =', mp.nstr(mp.quad(fU, [bm, 1, am-1, am]), 50))
print('NUMERICAL mp  V =', mp.nstr(mp.quad(fV, [bm, 1, am-1, am]), 50))
```

## 6. Session files (working directory `…/new_results/A3b_certified2`; the tool auto-numbers each code cell as `step_NN.py`)

- Exploration/reading: directory listing (`step_01.py`); reading `w4_symbol.py` (`conjecture_terms`, `w4_symbol`, `q` definitions) and `../A3b_certified/result.md` §§1–7; flint version/path checks (python-flint 0.9.0, `acb.log(s, analytic)` confirmed).
- API/pitfall checks (the cell printing `mul() -> IndexError`, `log(analytic=True) -> nan + nanj`, the `acb.integral` docstring, and `arb` `.rad()/.upper()/.abs_upper()` checks).
- `w4_certified.py` — written by one cell, then patched in place by a later cell (the `certified_decimal_exp` fix shown in §5); both cells are in the auto-numbered `step_NN.py` series (the patch cell is the one beginning `old = '''def certified_decimal_exp...`).
- `out_400.txt` — full stdout of `python w4_certified.py 400 500` (returncode 0, 2.2 s), the run quoted in §3; written by the midpoint-comparison cell (which itself had the harmless dict-key typo noted above).

**Label: PARTIAL / NUMERICAL.** Certified and printed here: `|c0| < 10^-112`, `|c1| < 10^-111` (PREC=400). Not done: the PREC=1100/EPS=2⁻¹³⁰⁰ run that would print the `k ≈ 300`-level balls — the pipeline is fixed, audited, and ready; one command remains.