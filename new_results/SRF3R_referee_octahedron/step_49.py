import json
from mpmath import mp, mpf, log, pi, sqrt, factorial, cos, sin, atan2
mp.dps = 60
z0 = (-log(2)+pi*1j)/2; R = abs(z0); th = mp.arg(z0)
j120 = json.load(open("../SRF3_octahedron_radius_proof/SRF3_hi_coeffs_K120.json"))
A_tau = [mpf(x) for x in j120['eta']]
bt = [A_tau[k]/factorial(k) for k in range(121)]

k = 100
u = bt[k]*R**k*k**mpf('1.5')
print("data u_100 =", mp.nstr(u,12))
C = sqrt(3*pi*R/8)
for phi_name, phi in [("+2.2476162 (= pi-th/2)", pi-th/2), ("-2.2476162", -(pi-th/2))]:
    print(f"model C cos(k th + phi), phi={phi_name}: {mp.nstr(C*cos(k*th+phi),12)}")
print("model -C cos((k-1/2) th):", mp.nstr(-C*cos((k-mpf('0.5'))*th),12))
# so the true phase phi_tau (in cos(k th + phi)) is:
A = -C*cos(th/2); B = -C*sin(th/2)
print("expand: A cos(kth)+B sin(kth): A =", mp.nstr(A,10), " B =", mp.nstr(B,10))
print("atan2(B,A) =", mp.nstr(atan2(B,A),12), "  2pi + atan2(B,A) =", mp.nstr(2*pi+atan2(B,A),12))
print("check cos(kth+atan2):", mp.nstr(C*cos(k*th+atan2(B,A)),12))
