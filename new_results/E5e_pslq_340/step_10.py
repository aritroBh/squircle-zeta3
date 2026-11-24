import json
from fractions import Fraction
import mpmath as mp

C = json.load(open("consts360.json"))
rel = [7680000, 0, 1175040000, 49152000, 2574432000, -501811200, -111744000, 927640800,
       6229440, 397036800, 45175200, 287114250, 381542400, 77518080, 52715520, -1240896,
       -1751040, -10005440, -2977497, -24330240, -56401920]
labels = ["X"] + C["labels"]
strs = [C["vals"]["X"]] + [C["vals"][l] for l in C["labels"]]

# ---- exact rational residual on 350-digit decimal strings (truncated data) ----
def to_frac(s):
    return Fraction(s)  # exact decimal parse
fr = [to_frac(s) for s in strs]
resF = sum(Fraction(c)*x for c, x in zip(rel, fr))
print("exact rational residual (350-digit truncated data):")
print("  numerator bits:", resF.numerator.bit_length())
af = abs(resF)
print("  |residual| =", mp.nstr(mp.mpf(af.numerator)/mp.mpf(af.denominator), 6))
print("  bound maxcoeff*1e-330 = 2.574432e-321 -> PASS:",
      mp.mpf(af.numerator)/mp.mpf(af.denominator) < mp.mpf('2.574432e-321'))

# ---- PSLQ at 340 digits (4th precision) ----
mp.mp.dps = 340
v = [mp.mpf(s) for s in strs]
import time; t0=time.time()
rel340 = mp.pslq(v, tol=mp.mpf(10)**(-315), maxcoeff=int(1e14), maxsteps=200000)
print("\nPSLQ @340 digits:", "FOUND" if rel340 else "none", f"({time.time()-t0:.0f}s)")
print("identical to 250/280/310 vector:", [int(c) for c in rel340] == rel)
