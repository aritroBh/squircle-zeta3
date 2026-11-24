import mpmath as mp, pickle, sys
sys.path.insert(0,'.')
import ftfield as F
mp.mp.dps = 80
m = 7
g = pickle.load(open('geom14.pkl','rb'))
LL = mp.sqrt(F.feval(g['L'],80))

def fsum(th, p):
    return sum(mp.exp(p*mp.log(abs(mp.cos(th - i*mp.pi/m)))) for i in range(m))
def r_and_rp(th, p):
    f = fsum(th,p)
    fp = sum(-mp.sin(th - i*mp.pi/m)*mp.sign(mp.cos(th-i*mp.pi/m))*p*mp.exp((p-1)*mp.log(abs(mp.cos(th-i*mp.pi/m)))) for i in range(m))
    r = f**(-1/p)
    rp = -(1/p)*f**(-1/p-1)*fp
    return r, rp
def integrand(th, p):
    r, rp = r_and_rp(th,p)
    return mp.sqrt(r*r+rp*rp)
def perim(p):
    return 2*m*mp.quad(lambda th: integrand(th,p), [-mp.pi/14, mp.pi/14], method='tanh-sinh', maxdegree=8)

# 9-point central difference 3rd derivative stencil, error O(h^8)
c3 = { -4: -mp.mpf(7)/240, -3: mp.mpf(3)/10, -2: -mp.mpf(169)/120, -1: mp.mpf(61)/30,
       1: -mp.mpf(61)/30, 2: mp.mpf(169)/120, 3: -mp.mpf(3)/10, 4: mp.mpf(7)/240 }
results = {}
for hs in ['0.02','0.01','0.005','0.0025']:
    h = mp.mpf(hs)
    vals = {k: perim(1+k*h) for k in range(-4,5)}
    P3 = sum(c3[k]*vals[k] for k in c3)/h**3
    results[hs] = P3
    print(f"h={hs}:  P3(1) = {mp.nstr(P3,40)}", flush=True)
# Richardson: assume P3(h) = P3 + A h^2 + B h^4 + C h^6 (central stencil even powers)
hs_list = [mp.mpf(x) for x in ['0.02','0.01','0.005','0.0025']]
Ps = [results[x] for x in ['0.02','0.01','0.005','0.0025']]
# fit polynomial in h^2
import itertools
X = [[h**0, h**2, h**4, h**6] for h in hs_list]
M = mp.matrix(X); b = mp.matrix(Ps)
sol = mp.lu_solve(M, b)
print("Richardson extrapolated P3(1) =", mp.nstr(sol[0],40))
print("chain value                     = -5.5080836980534783")
print("gold value                      = -5.5080732424212877")
open('progress.log','a').write('PG6 s_A2: Richardson P3(1) = %s\n'%mp.nstr(sol[0],30))
