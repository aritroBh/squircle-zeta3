
import mpmath as mp
mp.mp.dps = 400
import json, time

d = json.load(open("<project>/research/08_theory/D11b2_stokes_jump/step2_num450_N180.json"))
f = [mp.mpf(v) for v in d['f450']]
b = [f[n]/mp.factorial(n) for n in range(len(f))]  # Borel coefficients

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

t0 = time.time()
num, den = pade_coeffs(b, 90, 90)
print("pade build s:", round(time.time()-t0, 1))

def peval(t):
    nv = num[0]; dv = den[0]
    for k in range(1, len(num)): nv += num[k]*t**k
    for k in range(1, len(den)): dv += den[k]*t**k
    return nv/dv

# certification: reproduce b_n from the Pade approximant's series
worst = mp.mpf(0)
# series expansion of num/den to 181 terms
ser = [mp.mpf(0)]*181
ser[0] = num[0]
for n in range(1, 181):
    ser[n] = (num[n] if n < len(num) else mp.mpf(0)) - sum(den[j]*ser[n-j] for j in range(1, min(n, 90)+1))
for n in range(181):
    worst = max(worst, abs(ser[n]-b[n]))
print("Pade reproduces b_0..b_181 with max err:", mp.nstr(worst, 4))

json.dump({"num":[mp.nstr(v,380) for v in num], "den":[mp.nstr(v,380) for v in den]}, open("pade9090.json","w"))

# reference I(eta) by direct quadrature (the Laplace form)
def I_ref(eta):
    f2 = lambda u: mp.exp(-u)*(2*mp.exp(-eta*u/2)*mp.sqrt(mp.cosh(eta*u)))*(1+mp.exp(-u))**(-2-eta)
    return mp.quad(f2, [0, 1, 5, 20, 100])
for eta in ['0.1','0.2','0.3','-0.1','-0.2','-0.3']:
    v = I_ref(mp.mpf(eta))
    print(f"I_ref({eta}) = {mp.nstr(v, 60)}")
json.dump({eta: mp.nstr(I_ref(mp.mpf(eta)), 90) for eta in ['0.1','0.2','0.3','-0.1','-0.2','-0.3']}, open("iref.json","w"))
