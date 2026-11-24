from mpmath import mp, mpf, pi, sin, cos, log, sqrt, factorial
import json
mp.dps = 60
z0 = (-log(2)+pi*1j)/2; R = abs(z0); th = mp.arg(z0)
alpha = pi - th
bound = sin(alpha/2)
print("alpha = pi - theta =", mp.nstr(alpha,14), " sin(alpha/2) =", mp.nstr(bound,14))
# check: among any 3 consecutive k, one has |cos((k-1/2)th)| >= sin(alpha/2)
bad = 0; minmax = 1.0
for k in range(1, 200001):
    m = max(abs(cos((kk-mpf('0.5'))*th)) for kk in (k,k+1,k+2))
    minmax = min(minmax, float(m))
    if m < bound: bad += 1
print("windows k..k+2 (k<=200000) failing:", bad, "; weakest window max|cos| =", minmax, ">= sin(alpha/2)?", minmax >= float(bound))

# limsup behavior of actual coefficients
j120 = json.load(open("../SRF3_octahedron_radius_proof/SRF3_hi_coeffs_K120.json"))
bt = [mpf(j120['eta'][k])/factorial(k) for k in range(121)]
print("\n k : |a_k/k!|^{1/k} * R   (should approach 1; dips when cos small)")
for k in range(90,121):
    v = abs(bt[k])**(mpf(1)/k)*R
    print(f"{k:4d}: {mp.nstr(v,10)}")
