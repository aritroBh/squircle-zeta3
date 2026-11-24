import mpmath as mp, sympy as sp, pickle
mp.mp.dps = 60
N = 200
om = [mp.exp(2*mp.pi*mp.j*k/N) for k in range(N)]
m = 5
th = [mp.pi/2 - j*mp.pi/m for j in range(m)]
h_ = mp.pi/(2*m); R_ = mp.tan(h_)
V0 = (R_, mp.mpf(0)); V1 = (R_*mp.cos(2*h_), R_*mp.sin(2*h_))
L = mp.sqrt((V1[0]-V0[0])**2+(V1[1]-V0[1])**2)
tau_ = ((V1[0]-V0[0])/L, (V1[1]-V0[1])/L)
he = abs(V0[0]*V1[1]-V0[1]*V1[0])/L
s0 = V0[0]*tau_[0]+V0[1]*tau_[1]
def aj(j, t):
    y0 = (1-t)*V0[0]+t*V1[0]; y1 = (1-t)*V0[1]+t*V1[1]
    return abs(mp.cos(th[j])*y0 + mp.sin(th[j])*y1)
def apj(j):
    c0 = mp.cos(th[j])*V0[0]+mp.sin(th[j])*V0[1]
    c1 = mp.cos(th[j])*V1[0]+mp.sin(th[j])*V1[1]
    return (1 if c0 > 0 else -1)*(c1-c0)
def cauchy_all(t, R):
    sig = s0 + L*t
    lj = [mp.log(aj(j, t)) for j in range(m)]
    ajv = [aj(j, t) for j in range(m)]
    ap = [apj(j) for j in range(m)]
    vals = []; aF = aw = None
    for jj in range(N):
        d = R*om[jj]; p = 1+d
        F = mp.fsum(ajv[j]*mp.exp(d*lj[j]) for j in range(m))
        Phi = mp.fsum(ap[j]*mp.exp(d*lj[j]) for j in range(m))/F
        w2 = (L - sig*Phi)**2 + he**2*Phi**2
        a1 = mp.arg(F); a2 = mp.arg(w2)
        if aF is not None:
            a1 += 2*mp.pi*int(mp.nint((aF - a1)/(2*mp.pi)))
            a2 += 2*mp.pi*int(mp.nint((aw - a2)/(2*mp.pi)))
        aF, aw = a1, a2
        logF = mp.log(abs(F)) + 1j*a1; logw = mp.log(abs(w2)) + 1j*a2
        vals.append(mp.exp(logw/2 - logF/p))
    return [mp.re(mp.fsum(vals[jj]*om[jj]**(-k) for jj in range(N))/(N*R**k)) for k in range(6)]
t = mp.mpf('0.5')
for R in ['0.35', '0.175', '0.0875']:
    Ec = cauchy_all(t, mp.mpf(R))
    print(f"R={R}: E0={mp.nstr(Ec[0],10)} (L={mp.nstr(L,10)}) E3={mp.nstr(Ec[3],10)} E5={mp.nstr(Ec[5],12)}")
# check F at d=0 and Phi at d=0 directly
print("sum a_j =", mp.fsum(aj(j,t) for j in range(m)), " (should be 1)")
print("sum a_j' =", mp.fsum(apj(j) for j in range(m)), " (should be 0)")
print("L =", mp.nstr(L,15), "he =", mp.nstr(he,15), "s0 =", mp.nstr(s0,15), "alpha=", mp.nstr(aj(0,1),15), "beta=", mp.nstr(aj(1,0),15))
