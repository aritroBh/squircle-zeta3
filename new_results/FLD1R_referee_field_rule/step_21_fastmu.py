"""step21_fastmu.py -- FAST exact mu_{k,n} in Q[z2,...,zk], via the Stirling formula (no symbolic-n series).

[u^k]E_n = sum_{i+j=k} A^i/i! * ( sum_{m<=j} binom(1-2n,m) y^m m! S2(j,m) ) * l^j/j!,
A = -(n-1)a + n b,  l = b - a.  Each monomial y^j a^i b^l (j<=k, i+l<=k) integrates to
(n)_j/(2n)_j * D_{i,l}(n+j,n),  D = moment polynomials of the cumulants
  kap_{r,0} = psi^(r-1)(p)-psi^(r-1)(p+q), kap_{0,m} = psi^(m-1)(q)-psi^(m-1)(p+q),
  kap_{r,m} = -psi^(r+m-1)(p+q)  (r,m>=1);   p=n+j, q=n.
At integers: psi^(m)(N) = (-1)^(m+1) m! (zeta(m+1) - H^(m+1)_{N-1}) for m>=1;  psi(N)+gamma = H_{N-1}.
Only MIXED cumulants carry a zeta (single zeta(r+m)); pure ones are rational.
mu_{k,n} = -(1/(n(2n-1))) * sum ...
Output: dict {(e2,...,ek): Fraction} per (k,n); verified against flat_exact k<=3 and G1 table.
"""
import sys, json, time
from fractions import Fraction as Fr
from math import comb, factorial

def S2(j, m):  # Stirling 2nd kind
    if m > j or m < 0: return 0
    if j == 0: return 1 if m == 0 else 0
    return sum(((-1)**(m-t))*comb(m, t)*t**j for t in range(m+1))//factorial(m)

_H = {}
def H(r, m):
    key=(r,m)
    if key not in _H:
        _H[key] = sum((Fr(1, i**r) for i in range(1, m+1)), Fr(0))
    return _H[key]

# psi values at integer N as dict over zeta-monomials.  Represent element as dict {zeta_exp_tuple: Fraction}.
def padd(P,Q,s=1):
    R=dict(P)
    for k,v in Q.items():
        R[k]=R.get(k,0)+s*v
        if R[k]==0: del R[k]
    return R
def pmul(P,Q):
    R={}
    for k1,v1 in P.items():
        for k2,v2 in Q.items():
            kk=tuple(a+b for a,b in zip(k1,k2)); R[kk]=R.get(kk,0)+v1*v2
            if R[kk]==0: del R[kk]
    return R
def pscale(P,c):
    if c==0: return {}
    return {k:v*c for k,v in P.items()}
def const(x, K):
    return {tuple([0]*(K-1)): Fr(x)} if x else {}
def zeta(idx, K):  # zeta(idx), idx>=2 ; exponent tuple length K-1 (z2..zK)
    e=[0]*(K-1); e[idx-2]=1; return {tuple(e): Fr(1)}

def psi(m, N, K):
    """psi^{(m)}(N); for m=0 returns H_{N-1} - gamma (gamma as 0 -- it cancels in differences)."""
    if m == 0:
        return const(H(1, N-1), K)
    c = (-1)**(m+1)*factorial(m)
    return padd(pscale(zeta(m+1, K), c), const(-c*H(m+1, N-1), K))

def cumulant(r, m, p, q, K):
    out = pscale(psi(r+m-1, p+q, K), -1)
    if m == 0: out = padd(out, psi(r+m-1, p, K))
    if r == 0: out = padd(out, psi(r+m-1, q, K))
    return out

def moments(K):
    """all D_{i,l} for i+l<=K as functions of (p,q); returns dict (i,l)->(p,q)->val; we memo per (p,q)."""
    memopq = {}
    def D(i, l, p, q):
        key=(i,l,p,q)
        if key in memopq: return memopq[key]
        if i==0 and l==0:
            res = const(1, K)
        else:
            res = {}
            if i>0:
                for r in range(1, i+1):
                    for m in range(0, l+1):
                        t = pmul(cumulant(r,m,p,q,K), D(i-r, l-m, p, q))
                        res = padd(res, pscale(t, comb(i-1,r-1)*comb(l,m)))
            else:
                for m in range(1, l+1):
                    t = pmul(cumulant(0,m,p,q,K), D(0,l-m,p,q))
                    res = padd(res, pscale(t, comb(l-1,m-1)))
        memopq[key]=res
        return res
    return D

def Ek_coeffs(k, n):
    """[u^k]E_n as dict {(j,ia,lb): Fraction}; j=y-exp, ia=a-exp, lb=b-exp."""
    out = {}
    def addterm(j,ia,lb,c):
        out[(j,ia,lb)] = out.get((j,ia,lb),Fr(0))+c
        if out[(j,ia,lb)]==0: del out[(j,ia,lb)]
    for i in range(k+1):
        j = k-i
        Ai = {}  # A^i as dict {(ia,lb): coef(n)*} with numeric n: A = -(n-1)a+nb
        # expand (-(n-1)a+nb)^i
        for ta in range(i+1):
            tb = i-ta
            Ai[(ta,tb)] = Fr(comb(i,ta))*Fr(-(n-1))**ta * Fr(n)**tb
        # binom(1-2n, m) general: product (1-2n)(1-2n-1)...(1-2n-m+1)/m!
        inner=[]
        for m in range(j+1):
            bc=Fr(1)
            for t in range(m): bc *= Fr(1-2*n-t, t+1)
            inner.append((m, bc*factorial(m)*S2(j,m)))
        lj = [(ta, Fr(comb(j,ta))*Fr(-1)**ta) for ta in range(j+1)]  # (b-a)^j = sum C(j,ta) b^{j-ta}(-a)^ta
        fi = Fr(1, factorial(i)); fj = Fr(1, factorial(j))
        for (ia,ib),ca in Ai.items():
            for (m, cm) in inner:
                for (ta, cl) in lj:
                    # a-exp = ia + ta, b-exp = ib + (j-ta)
                    addterm(m, ia+ta, ib+(j-ta), ca*cm*cl*fi*fj)
    return out

def mu(k, n, K, D):
    tot = {}
    for (j, ia, lb), co in Ek_coeffs(k, n).items():
        pr = Fr(1)
        for t in range(j): pr *= Fr(n+t, 2*n+t)
        term = pscale(D(ia, lb, n+j, n), co*pr)
        tot = padd(tot, term)
    return pscale(tot, Fr(-1, n*(2*n-1)))

if __name__ == '__main__':
    KMAX = int(sys.argv[1]) if len(sys.argv)>1 else 5
    NMAX = int(sys.argv[2]) if len(sys.argv)>2 else 100
    t0=time.time()
    D = moments(KMAX)
    table = {}
    for k in range(KMAX+1):
        table[k]={}
        for n in range(1, NMAX+1):
            P = mu(k, n, KMAX, D)
            table[k][n] = {','.join(map(str,e)): str(v) for e,v in P.items()}
    print(f'computed mu_{{k,n}} for k<={KMAX}, n<={NMAX} in {time.time()-t0:.1f}s')
    # verify k<=3 against G1 table flat_exact_N320_K3.json (n<=100); normalize arity
    G1='<project>/research/10_frontier/G1_zeta3_any_angle/'
    ref = json.load(open(G1+'flat_exact_N320_K3.json'))
    def normkey(kk):
        comps = [int(x) for x in kk.strip('()').split(',')]
        comps = (comps + [0,0,0,0])[:4]
        return ','.join(map(str,comps))
    ok=True
    for k in (0,1,2,3):
        for n in range(1,101):
            mine = table[k][n]
            theirs = {normkey(kk): vv for kk,vv in ref[str(k)][str(n)].items()}
            for kk in set(list(mine)+list(theirs)):
                if Fr(mine.get(kk,'0')) != Fr(theirs.get(kk,'0')):
                    ok=False; print('MISMATCH',k,n,kk,mine.get(kk),theirs.get(kk))
    print('k<=3 matches G1 flat_exact table (n<=100):', ok)
    json.dump({str(k):{str(n):table[k][n] for n in table[k]} for k in table}, open(f'mu_exact_K{KMAX}_N{NMAX}.json','w'))
    print('saved', f'mu_exact_K{KMAX}_N{NMAX}.json')
    # print a few
    for k in (2,3,4,5):
        print(f'mu_{{{k},2}} =', table[k][2])
