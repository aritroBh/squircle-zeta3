
import mpmath as mp
import json

P = json.load(open("pade9090.json"))

def de_quad(phi, h, Tlo, Thi):
    tot = mp.mpf(0)
    klo = int(Tlo/h)+1; khi = int(Thi/h)+1
    for k in range(-klo, khi+1):
        t = k*h
        sh = mp.sinh(t); u = mp.exp(mp.pi/2*sh)
        up = mp.pi/2*mp.cosh(t)*u
        tot += phi(u)*up
    return tot*h

refs = None
mp.mp.dps = 110
refs = {k: mp.mpf(v) for k, v in json.load(open("iref_de160.json")).items()}

# (a) conditioning: evaluate PB at zeta=-23,-45 at dps 110 vs 260
for z in ['-23','-45','30','45']:
    vals = []
    for dps in [110, 260]:
        mp.mp.dps = dps
        num = [mp.mpf(v) for v in P['num']]; den = [mp.mpf(v) for v in P['den']]
        zz = mp.mpf(z)
        nv = sum(num[k]*zz**k for k in range(len(num)))
        dv = sum(den[k]*zz**k for k in range(len(den)))
        vals.append(nv/dv)
    mp.mp.dps = 110
    print(f"PB({z}): dps110 vs dps260 diff = {mp.nstr(abs(vals[0]-mp.mpf(str(vals[1]))),3)}")

# (b) build PB[70/70] at dps 400 and compare the eta=-0.2,-0.3 sums
mp.mp.dps = 400
d = json.load(open("<project>/research/08_theory/D11b2_stokes_jump/step2_num450_N180.json"))
f = [mp.mpf(v) for v in d['f450']]
b = [f[n]/mp.factorial(n) for n in range(len(f))]
def pade_coeffs(a, L, M):
    A = mp.matrix(M, M); rhs = mp.matrix(M, 1)
    for i in range(M):
        for j in range(M):
            k = L + i - j
            A[i, j] = a[k] if k >= 0 else mp.mpf(0)
        rhs[i] = -a[L + 1 + i]
    q = (A**-1)*rhs
    den = [mp.mpf(1)] + [q[j] for j in range(M)]
    num = [sum(den[j]*a[n-j] for j in range(0, min(M, n)+1)) for n in range(L+1)]
    return num, den
num70, den70 = pade_coeffs(b, 70, 70)
json.dump({"num":[mp.nstr(v,380) for v in num70], "den":[mp.nstr(v,380) for v in den70]}, open("pade7070.json","w"))
mp.mp.dps = 110
num70 = [mp.mpf(v) for v in json.load(open("pade7070.json"))['num']]
den70 = [mp.mpf(v) for v in json.load(open("pade7070.json"))['den']]
def peval70(t):
    nv = num70[0]; dv = den70[0]
    for k in range(1, len(num70)): nv += num70[k]*t**k
    for k in range(1, len(den70)): dv += den70[k]*t**k
    return nv/dv
num = [mp.mpf(v) for v in P['num']]; den = [mp.mpf(v) for v in P['den']]
def peval90(t):
    nv = num[0]; dv = den[0]
    for k in range(1, len(num)): nv += num[k]*t**k
    for k in range(1, len(den)): dv += den[k]*t**k
    return nv/dv
Tlo = mp.asinh(2/mp.pi*mp.log(mp.mpf('1e60')))
Thi = mp.asinh(2/mp.pi*mp.log(160))
for eta in ['-0.2','-0.3','0.3']:
    e = mp.mpf(eta)
    v90 = de_quad(lambda u: mp.exp(-u)*peval90(e*u), mp.mpf(1)/64, Tlo, Thi)
    v70 = de_quad(lambda u: mp.exp(-u)*peval70(e*u), mp.mpf(1)/64, Tlo, Thi)
    print(f"eta={eta}: S90-S70={mp.nstr(abs(v90-v70),3)}  |S90-I|={mp.nstr(abs(v90-refs[eta]),3)}  |S70-I|={mp.nstr(abs(v70-refs[eta]),3)}")
