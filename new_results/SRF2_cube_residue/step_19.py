import json, mpmath as mp
mp.mp.dps = 60
base = "<project>"
d = json.load(open(base + "/research/08_theory/SRF1_radii_3d/EV_coeffs_K30.json"))
V = [mp.mpf(x) for x in d["V"]]
E3 = [mp.mpf(x) for x in d["E3"]]
E2 = [mp.mpf(x) for x in d["E2"]]
R_planner = -1 - 2*mp.pi**mp.mpf('1.5')/mp.gamma(mp.mpf('0.25'))**2
print("R_planner =", mp.nstr(R_planner, 40))
# alternating limit of V_k
for k in range(3, 31):
    print(k, mp.nstr((-1)**k*V[k], 25))
