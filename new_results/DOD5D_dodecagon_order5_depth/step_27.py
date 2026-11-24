import sys
sys.path.insert(0, "<project>/research/08_theory/E5c_nonlocal_z5")
import mpmath as mp
mp.mp.dps = 50
import symeng
r3 = mp.sqrt(3)
# signed normals for the edge in wedge (0, pi/6): signs (+,+,+,+,-,-)
M = []
for j in range(6):
    c, s = mp.cos(j*mp.pi/6), mp.sin(j*mp.pi/6)
    sg = 1 if j <= 3 else -1
    M.append((sg*c, sg*s))
y0 = (2-r3, mp.mpf(0))
y1 = ((2*r3-3)/2, (2-r3)/2)
P = symeng.edge_params(M, y0, y1)
for k in ['L','s0','h','alpha','beta','j0','j1','nonloc']:
    print(k, "=", mp.nstr(P[k], 30) if not isinstance(P[k],(int,list)) else P[k])
print("c0:", [mp.nstr(c, 20) for c in P['c0']])
print("c1:", [mp.nstr(c, 20) for c in P['c1']])
kappas = [(P['c1'][j]-P['c0'][j])/P['c0'][j] for j in P['nonloc']]
print("nonloc:", P['nonloc'])
print("kappas:", [mp.nstr(k, 30) for k in kappas])
print("letters -1/kappa:", [mp.nstr(-1/k, 30) for k in kappas])
# exact forms
print("exact kappa check:")
print(mp.nstr(abs(kappas[0]-(r3-2)/2),4), mp.nstr(abs(kappas[1]-(2*r3-3)/3),4),
      mp.nstr(abs(kappas[2]-(r3-1)),4), mp.nstr(abs(kappas[3]-(r3-3)/3),4))
print("exact letter check:")
print(mp.nstr(abs(-1/kappas[0]-(4+2*r3)),4), mp.nstr(abs(-1/kappas[1]+(3+2*r3)),4),
      mp.nstr(abs(-1/kappas[2]+(1+r3)/2),4), mp.nstr(abs(-1/kappas[3]-(3+r3)/2),4))
# L, h, s0 exact
u = mp.sqrt(2-r3)
print("L = u*(2-r3)?", mp.nstr(abs(P['L']-u*(2-r3)),4), " h = u/2?", mp.nstr(abs(P['h']-u/2),4))
print("s0 =", mp.nstr(P['s0'], 30), " s0/L =", mp.nstr(P['s0']/P['L'],30))
