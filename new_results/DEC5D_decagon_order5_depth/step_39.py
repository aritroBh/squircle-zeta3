
import sys, mpmath as mp
mp.mp.dps = 50
base = "<project>/research/08_theory/"
sys.path.insert(0, base + "E5c_nonlocal_z5")
import symeng as SE

sqrt5 = mp.sqrt(5); phi = (1+mp.sqrt(5))/2
h = 1/(2*phi)
R = 1/(2*phi**2*mp.sin(mp.pi/5))
Ry = R*mp.sin(mp.pi/10)   # = R/(2 phi)
V0 = [-h, Ry]; V1 = [-h, -Ry]
M = [[mp.cos(j*mp.pi/5), mp.sin(j*mp.pi/5)] for j in range(5)]
P = SE.edge_params(M, V0, V1)
for k in ['L','s0','h','alpha','beta','j0','j1','nonloc']:
    print(k, "=", mp.nstr(P[k], 30))
print("c0 =", [mp.nstr(c, 20) for c in P['c0']])
print("c1 =", [mp.nstr(c, 20) for c in P['c1']])
# kappas for nonlocal terms: a_j(t) = c0_j(1-t)+c1_j t = B_j(1+kap_j t)
for j in P['nonloc']:
    c0, c1 = P['c0'][j], P['c1'][j]
    print(f"j={j}: c0={mp.nstr(c0,20)} c1={mp.nstr(c1,20)} B={mp.nstr(abs(c0),20)} kap={mp.nstr(c1/c0-1,20)}")
print("\nexact target values:")
print("alpha=beta=1/(2 phi^2) =", mp.nstr(1/(2*phi**2), 30))
print("L = sqrt((25-11 sqrt5)/10) =", mp.nstr(mp.sqrt((25-11*sqrt5)/10), 30))
print("h = 1/(2 phi) =", mp.nstr(1/(2*phi), 30), " s0/L = -1/2? s0=", mp.nstr(P['s0'],20))
print("kappa_1 = 1/phi =", mp.nstr(1/phi, 30), " kappa_4 = -1/phi^2 =", mp.nstr(-1/phi**2, 30))
