from mpmath import (mp, mpf, pi, gamma, sin, cos, exp, log, nsum, inf, nstr, fac, zeta)
mp.dps = 100
def N(x,d=30): return nstr(x,d)
def logline(msg):
    with open("progress.log","a") as f: f.write(msg+"\n")
import sympy as sp
z3, z5, z7, p_ = sp.symbols('z3 z5 z7 pi')
def sym_zeta(k):
    if k==3: return z3
    if k==5: return z5
    if k==7: return z7
    return sp.zeta(k)
def sym_eta(k): return (1-2**(1-k))*sym_zeta(k)
pi_s = p_
def terms_a(s):
    if s==1: return [(2/pi_s,'a',2),(-2/pi_s,'c',2)]
    if s==2: return [(4,'a',2)]
    if s==3: return [(6*pi_s,'a',2),(12/pi_s,'c',4),(-12/pi_s,'a',4)]
    if s==4: return [(8*pi_s**2,'a',2),(-48,'a',4)]
def terms_b(s):
    if s==1: return [(-2,'a',1)]
    if s==2: return [(-2*pi_s,'a',1),(-4/pi_s,'c',3),(4/pi_s,'a',3)]
    if s==3: return [(-2*pi_s**2,'a',1),(12,'a',3)]
    if s==4: return [(-2*pi_s**3,'a',1),(24*pi_s,'a',3),(48/pi_s,'c',5),(-48/pi_s,'a',5)]
def Ssum(w,k): return sym_zeta(k) if w=='c' else -sym_eta(k)
def B_sym(s,sp):
    tot = 0
    for (c1,w1,k1) in terms_a(s):
        for (c2,w2,k2) in terms_b(sp):
            w = 'c' if w1==w2 else 'a'
            tot += c1*c2*Ssum(w,k1+k2)
    return sp.simplify(tot)
pairs = [(1,1),(1,2),(2,1),(2,2),(1,3),(2,3),(3,3),(1,4),(2,4),(3,4),(4,4)]
for (s,sp_) in pairs:
    expr = sp.simplify(B_sym(s,sp_))
    print(f"B({s},{sp_}) = {expr}")
logline("step_35: B(s,s') closed forms via sympy")
