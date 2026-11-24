# Planner check: the tetrakis T2 formula (with the c_A mixing) against the model-A grids h192hp and h256hp.
import json, mpmath as mp
mp.mp.dps = 120
S0 = 2*mp.sqrt(5)/3
la = mp.log(mp.mpf(2)**(-mp.mpf(6)/5)*mp.mpf(3)**(mp.mpf(13)/5))
L = mp.log(mp.mpf(4)/3)
T2 = mp.sqrt(5)*(mp.mpf(254)/675 - 79*mp.pi**2/8100 + mp.mpf(118)/75*mp.polylog(2, -mp.mpf(1)/3) + mp.mpf(4)/15*L - mp.mpf(76)/75*L**2)
pred1 = S0*(mp.mpf(13)/15 + la)
pred2 = T2 + S0*mp.mpf(13)/15*la + S0*(la**2/2 - la)
for tag in ("h192hp", "h256hp"):
    S = json.load(open(f"modelA_{tag}.json"))["Sdk"]
    print(tag, "d1 resid", mp.nstr(mp.mpf(S[1]) - pred1, 3), " d2 resid", mp.nstr(mp.mpf(S[2]) - pred2, 3))
