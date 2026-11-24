
import mpmath as mp
mp.mp.dps = 70
import pickle
from fractions import Fraction as Fr

def Dtilde_GF(tau, m):
    tau = mp.mpf(tau)
    def f(y):
        t = (1-y)/(1+y)
        prof = mp.sqrt((1+y)/2) * (1 + t**(1/(1-tau)))**((1-tau)/2) - 1
        return prof * mp.asin(y)**(2*m) / mp.sqrt(1-y**2)
    return 2*mp.sqrt(2)*mp.quad(f, [0, mp.mpf('0.5'), 1])

K = 9
tau0 = mp.mpf('0.004')
taus = [tau0*(j+1) for j in range(K)]
def extract(m):
    vals = [Dtilde_GF(t, m) for t in taus]
    # solve Vandermonde: sum_{n=2}^{10} c_n tau_j^n = val_j  -> use n=2..10 (9 unknowns)
    nn = list(range(2, 11))
    A = mp.matrix(K, len(nn))
    for i, t in enumerate(taus):
        for j, n in enumerate(nn):
            A[i,j] = t**n
    b = mp.matrix(vals)
    sol = mp.lu_solve(A, b)
    return {nn[j]: sol[j] for j in range(len(nn))}

Dtarget = pickle.load(open('Dtilde_relations.pkl','rb'))
def zval_target(n, m):
    coeffs = [Fr(c) for c in Dtarget[f"{n},{m}"]]
    if n==2: B=[mp.zeta(2)]
    if n==3: B=[mp.zeta(3)]
    if n==4: B=[mp.zeta(4)]
    if n==5: B=[mp.zeta(5), mp.zeta(2)*mp.zeta(3)]
    if n==6: B=[mp.zeta(6), mp.zeta(3)**2]
    if n==7: B=[mp.zeta(7), mp.zeta(2)*mp.zeta(5), mp.zeta(2)**2*mp.zeta(3)]
    return sum(mp.mpf(c.numerator)/mp.mpf(c.denominator)*b for c,b in zip(coeffs,B))

for m in [0, 1]:
    cs = extract(m)
    for n in range(2, 8):
        if n < 2*m+2: continue
        try: tgt = zval_target(n, m)
        except Exception: continue
        print(f"m={m} n={n}: extracted={mp.nstr(cs[n],18)} target={mp.nstr(tgt,18)} diff={mp.nstr(abs(cs[n]-tgt),4)}")
