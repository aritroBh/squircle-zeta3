import sys, os
work = '<project>/research/08_theory/TOW2_tower_all_depths'
sys.path.insert(0, work); os.chdir(work)
from towlib import *
from math import comb

def e_from_data(key, delta, s):
    i,j = key; dd = i+j
    base = (2*i+1) if i==j else dd
    return sum( u(i,j,base+t, delta).scale(Fr((-1)**(s-t)*comb(s,t))) for t in range(0,s+1) )

def Ediag_data(key, delta, sigma):
    i,j = key
    if i != j: return e_from_data(key, delta, sigma)
    out = ZERO; t = 0
    while sigma+t <= delta//2:
        out = out + e_from_data(key, delta, sigma+t).scale(Fr((-1)**t))
        t += 1
    return out

def LambdaE_data(i, j, delta, sigma):
    d = i+j
    tot = ZERO
    for dp in range(0, delta+1):
        c1 = tau1(delta-dp).scale(Fr(i+1, 2*d-dp-1))
        tot = tot + e_from_data((i+1,j), dp, sigma) * c1
        if i > j:
            c2 = tau2(delta-dp).scale(Fr(j+1, 2*d-dp-1))
            tot = tot + Ediag_data((i,j+1), dp, sigma) * c2
    if j == i-1:
        mm = 4*i-4-delta
        ft = ZERO
        for mu in range(mm, 4*i-3):
            g = ZERO
            for dp in range(0, 4*i-4-mu+1):
                g = g - Ediag_data((i,i), dp, sigma) * tau1(4*i-4-mu-dp).scale(Fr(1, 4*i-3-dp))
            ft = ft + g.conj_i() * (IQ**(-mu)) * Fr(comb(mu-2, mu-mm)) * (IQ**(mm-2))
        tot = tot + ft.scale(-i)
    return tot

for (delta, s) in [(4,0),(4,1),(4,2)]:
    bad = 0; n = 0
    for pp in range(2, 6):
        key = (pp, pp-1); dd = 2*pp-1
        if delta > 2*dd-4 or dd+s > 14: continue
        n += 1
        lhs = e_from_data(key, delta, s) + e_from_data(key, delta, s+1)
        rhs = LambdaE_data(pp, pp-1, delta, s)
        if lhs != rhs: bad += 1
    print(f'adjacent ({delta},{s}) vs DATA: {bad} bad of {n}')
