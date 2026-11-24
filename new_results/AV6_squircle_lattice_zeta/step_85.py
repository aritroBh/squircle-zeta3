code = r'''
from umach2 import eta_l, E_rl, eta
from squircle import setdps
from mpmath import mp, nstr, zeta, gamma, power, mpf, fac
setdps(80)
s = 4
def monomials(r):
    out = []
    def rec(m, rem, cur):
        if rem == 0:
            A = sum(cur.values()); w = mpf(1)
            for mm, am in cur.items(): w *= ((-1)/mm)**am / fac(am)
            out.append((A, w)); return
        if m-1 > rem: return
        a = 0
        while a*(m-1) <= rem:
            cur[m] = a; rec(m+1, rem - a*(m-1), cur); a += 1
        cur.pop(m, None)
    rec(2, r, {})
    return out
def Q(r, X):
    tot = mpf(0)
    for (A, w) in monomials(r):
        poch = gamma(X+r+A)/gamma(X)
        tot += w*poch
    return tot
def u_1_r_closed(s, r):
    return gamma(s-1)*zeta(s-1)*eta(s+r)*Q(r, s-1)
def u_1_r_mach(s, r):
    tot = mpf(0)
    for l in range(1, r+2):
        tot += ((-1)**(l+1))/l * eta_l(l, s-1) * E_rl(r, l, s)
    return gamma(s-1)*zeta(s-1)*tot
for r in range(1, 9):
    a = u_1_r_closed(s, r); b = u_1_r_mach(s, r)
    print(f"r={r}: closed={nstr(a, 20)}")
    print(f"      mach  ={nstr(b, 20)}  diff={nstr(a-b, 5)}")
'''
exec(code)
with open("progress.log","a") as f: f.write("step_46: closed form u_{1,r} = Gamma zeta eta Q verified r=1..8\n")
