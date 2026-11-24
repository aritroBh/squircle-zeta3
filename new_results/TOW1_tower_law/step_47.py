verify = '''

def U_e(key, k, delta):
    """tower-law candidate from e-coefficients; diagonal at k=2p uses rho (initial)."""
    i, j = key; d = i+j
    if delta > 2*d-4: return None
    if i == j:
        p = i
        if k == 2*p:
            return rho_c(key, delta)
        return sum( (ehat_pp(p, delta, s).scale(comb(k-2*p-1, s)) for s in range(0, delta//2+1)), ZERO )
    return sum( (e_off(key, delta, s).scale(comb(k-d, s)) for s in range(0, delta//2+1)), ZERO )

def check_tower_law(kmax=14, dmax=10):
    bad = []
    for k in range(2, kmax+1):
        for d in range(2, min(dmax, k)+1):
            for key in keys_of_degree(d):
                i, j = key
                for delta in range(0, min(2, 2*d-4)+1):
                    ue = U_e(key, k, delta)
                    actual = u(i,j,k,delta)
                    onlaw = (i > j) or (k >= 2*i+1)
                    if onlaw and ue != actual:
                        bad.append((key,k,delta,ue,actual))
    return bad

def check_newton(dmax=9):
    """k-free Newton identities: e(delta,sigma)+e(delta,sigma+1) == LambdaE_ij(i,j,delta,sigma),
       for off-diagonal keys; diagonal keys: ehat == parent formula (by construction, but check vs data)."""
    bad = []
    for d in range(2, dmax+1):
        for key in keys_of_degree(d):
            i, j = key
            if i == j: continue
            for delta in (0,1,2):
                if delta > 2*d-4: continue
                for sigma in range(0, delta//2+1):
                    lhs = e_off(key, delta, sigma) + e_off(key, delta, sigma+1)
                    rhs = LambdaE_ij(i, j, delta, sigma)
                    if lhs != rhs:
                        bad.append((key,delta,sigma,lhs,rhs))
    return bad
'''
open('towlib.py','a').write(verify)
import importlib, towlib
importlib.reload(towlib)
from towlib import *
bad1 = check_tower_law()
print('tower law delta<=2, all keys k<=14: bad =', len(bad1))
for b in bad1[:6]: print(b)
