rec = '''

def fold_feed(i, k, delta):
    """the +i order-(2d-2-delta) coefficient of -i*F(Gm_ii DL1) at level k (d=2i-1), from tower data.
       = -i * sum_{mu=m}^{4i-4} f_mu (-i)^mu binom(mu-2, mu-m) i^{m-2}, f_mu = conj(g_mu),
       g_mu = -sum_{dp=0}^{4i-4-mu} u_ii(k,dp) tau1[4i-4-mu-dp]/(4i-3-dp)."""
    d = 2*i-1
    m = 4*i-4-delta
    tot = ZERO
    from math import comb
    for mu in range(m, 4*i-3):   # mu <= 4i-4
        g = ZERO
        for dp in range(0, 4*i-4-mu+1):
            g = g - u(i,i,k,dp) * tau1(4*i-4-mu-dp).scale(Fr(1, 4*i-3-dp))
        f = g.conj_i()
        tot = tot + f * (IQ**(-mu)) * Fr(comb(mu-2, mu-m)) * (IQ**(m-2))
        # careful: (-i)^mu = IQ^{-mu}? (-i) = i^{-1}? i^{-1} = -i. yes (-i)^mu = i^{-mu}.
    return tot.scale(-1) * IQ   # times -i

def Lambda_ij(i, j, k, delta, U):
    """the recursion feed (k>d part): value of (Lambda_d U_{d+1})(i,j) at depth delta,
       U(key, dp) = degree-(d+1) depth-dp data at level k."""
    d = i+j
    tot = ZERO
    for dp in range(0, delta+1):
        c1 = tau1(delta-dp).scale(Fr(i+1, 2*d-dp-1))
        tot = tot + U((i+1,j), dp) * c1
        if i > j:
            c2 = tau2(delta-dp).scale(Fr(j+1, 2*d-dp-1))
            tot = tot + U((i,j+1), dp) * c2
    if j == i-1:
        # fold term: need g_mu built from U((i,i), dp)
        from math import comb
        m = 4*i-4-delta
        ft = ZERO
        for mu in range(m, 4*i-3):
            g = ZERO
            for dp in range(0, 4*i-4-mu+1):
                g = g - U((i,i), dp) * tau1(4*i-4-mu-dp).scale(Fr(1, 4*i-3-dp))
            ft = ft + g.conj_i() * (IQ**(-mu)) * Fr(comb(mu-2, mu-m)) * (IQ**(m-2))
        tot = tot + (ft.scale(-i) * IQ)
    return tot
'''
open('towlib.py','a').write(rec)
import importlib, towlib
importlib.reload(towlib)
from towlib import *

# verify recursion: for all keys, levels k>d, depths delta: u_ij(k,delta) == Lambda_ij with U = u(k)
bad = []
tested = 0
for k in range(3, 13):
    for d in range(2, k):
        for key in keys_of_degree(d):
            i,j = key
            for delta in range(0, 2*d-3):
                lhs = u(i,j,k,delta)
                rhs = Lambda_ij(i,j,k,delta, lambda kk,dp: u(kk[0],kk[1],k,dp))
                tested += 1
                if lhs != rhs:
                    bad.append((k,key,delta,lhs,rhs))
print('recursion T1 tested:', tested, 'bad:', len(bad))
for b in bad[:5]: print(b)
