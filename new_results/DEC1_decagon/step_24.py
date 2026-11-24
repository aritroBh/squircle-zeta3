import os, mpmath as mp
print("CWD:", os.getcwd())
mp.mp.dps = 50
# normals n_i = (cos(i pi/5), sin(i pi/5)), i=1..5
ths = [mp.pi*i/5 for i in range(1,6)]
ns = [(mp.cos(t), mp.sin(t)) for t in ths]
def N(x):
    return sum(abs(n[0]*x[0]+n[1]*x[1]) for n in ns)
# vertex on kink line j: x = c * perp(n_j), N(x)=1
def vertex(j):
    nj = ns[j]
    perp = (-nj[1], nj[0])
    s = sum(abs(n[0]*perp[0]+n[1]*perp[1]) for n in ns)
    return (perp[0]/s, perp[1]/s)
import numpy as np
Vs = [vertex(j) for j in range(5)]
for j,V in enumerate(Vs):
    ang = mp.atan2(V[1],V[0])*180/mp.pi
    print(f"kink {j+1}: V=({mp.nstr(V[0],30)},{mp.nstr(V[1],30)})  ang={mp.nstr(ang,12)}  |V|={mp.nstr(mp.sqrt(V[0]**2+V[1]**2),20)}")
# edge from vertex on kink 2 (index1) to kink 3 (index2)
V0, V1 = Vs[1], Vs[2]
print("V0:", [mp.nstr(v,25) for v in V0])
print("V1:", [mp.nstr(v,25) for v in V1])
def avals(t):
    y = (V0[0]+t*(V1[0]-V0[0]), V0[1]+t*(V1[1]-V0[1]))
    return [abs(n[0]*y[0]+n[1]*y[1]) for n in ns]
a0 = avals(0); a1 = avals(1); am = avals(mp.mpf('0.37'))
print("a(V0):", [mp.nstr(v,20) for v in a0])
print("a(V1):", [mp.nstr(v,20) for v in a1])
print("ratios r_j = a(1)/a(0):", [mp.nstr(a1[j]/a0[j],20) if a0[j]!=0 else 'nan' for j in range(5)])
# check affine-ness + constancy of term 5 (index 4, normal (cos180,sin180)=(-1,0))
for j in range(5):
    lin = am[j] - (a0[j] + mp.mpf('0.37')*(a1[j]-a0[j]))
    print(j+1, "affine residual:", mp.nstr(lin,6))
phi = (1+mp.sqrt(5))/2
print("phi:", mp.nstr(phi,30))
print("r1 vs phi:", mp.nstr(a1[0]/a0[0]-phi,6), " r4 vs 1/phi:", mp.nstr(a1[3]/a0[3]-1/phi,6))
print("a5 const check a5(0),a5(1),a5(0.37):", mp.nstr(a0[4],25), mp.nstr(a1[4],25), mp.nstr(am[4],25))
print("a5 vs 1/(2phi):", mp.nstr(a0[4]-1/(2*phi),8))
print("B1 vs 1/(2phi^2):", mp.nstr(a0[0]-1/(2*phi**2),8))
print("B4 vs 1/(2phi):", mp.nstr(a0[3]-1/(2*phi),8))
print("A0(vanish@0 term idx1) a(V1) vs 1/(2phi^2):", mp.nstr(a1[1]-1/(2*phi**2),8))
print("A1(vanish@1 term idx2) a(V0) vs 1/(2phi^2):", mp.nstr(a0[2]-1/(2*phi**2),8))