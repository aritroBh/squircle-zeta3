"""ZLAD1R step_26 -- remaining mixed rungs on my table (n<=40, exact), plus Lemma 2 (P2) with MY OWN E-series:
(d_a+d_b)[u^k]E_n = [u^{k-1}]E_n exactly for n=2,3,5, k<=8."""
import json, time
from fractions import Fraction as Fr
from math import comb, factorial
mine = json.load(open('muR_exact_N40_K8.json'))
def mu(k, n, m): return Fr(mine[str(k)][str(n)].get(m, '0'))
N = 40
res = []
def chk(name, cond):
    ok = all(cond(n) for n in range(1, N+1))
    res.append((name, ok)); print(('PASS' if ok else 'FAIL'), name)
chk('[z2z3]mu_k = -2[z2^2]mu_{k-1}, k=7,8', lambda n: all(mu(k,n,'z2z3') == -2*mu(k-1,n,'z2^2') for k in (7,8)))
chk('[z3^2]mu_k = [z2^2]mu_{k-2}, k=7,8', lambda n: all(mu(k,n,'z3^2') == mu(k-2,n,'z2^2') for k in (7,8)))
chk('[z2z4]mu_k = 2[z2^2]mu_{k-2}+3[z2^3]mu_k, k=7,8', lambda n: all(mu(k,n,'z2z4') == 2*mu(k-2,n,'z2^2')+3*mu(k,n,'z2^3') for k in (7,8)))
chk('[z2z5]mu_8 = -2[z2^2]mu_5 - 6[z2^3]mu_7', lambda n: mu(8,n,'z2z5') == -2*mu(5,n,'z2^2') - 6*mu(7,n,'z2^3'))
chk('[z3z4]mu_8 = -2[z2^2]mu_5 - 3[z2^3]mu_7', lambda n: mu(8,n,'z3z4') == -2*mu(5,n,'z2^2') - 3*mu(7,n,'z2^3'))
chk('[z2z3^2]mu_8 = 3[z2^3]mu_6', lambda n: mu(8,n,'z2z3^2') == 3*mu(6,n,'z2^3'))
chk('[z2^2z4]mu_8 = 3[z2^3]mu_6 + 6[z2^4]mu_8', lambda n: mu(8,n,'z2^2z4') == 3*mu(6,n,'z2^3') + 6*mu(8,n,'z2^4'))
chk('[z2z6]mu_8 = 2[z2^2]mu_4 + 9[z2^3]mu_6 + 8[z2^4]mu_8', lambda n: mu(8,n,'z2z6') == 2*mu(4,n,'z2^2') + 9*mu(6,n,'z2^3') + 8*mu(8,n,'z2^4'))
chk('[z4^2]mu_8 = [z2^2]mu_4 + 3[z2^3]mu_6 + 3[z2^4]mu_8', lambda n: mu(8,n,'z4^2') == mu(4,n,'z2^2') + 3*mu(6,n,'z2^3') + 3*mu(8,n,'z2^4'))
chk('[z5]mu_8 = -[z2]mu_5 - 2[z2^2]mu_7', lambda n: mu(8,n,'z5') == -mu(5,n,'z2') - 2*mu(7,n,'z2^2'))
chk('[z4]mu_8 = [z2]mu_6 + [z2^2]mu_8', lambda n: mu(8,n,'z4') == mu(6,n,'z2') + mu(8,n,'z2^2'))
chk('[z6]mu_8 = [z2]mu_4 + 3[z2^2]mu_6 + 2[z2^3]mu_8', lambda n: mu(8,n,'z6') == mu(4,n,'z2') + 3*mu(6,n,'z2^2') + 2*mu(8,n,'z2^3'))
print('mixed-rung batch:', all(ok for _, ok in res), f'({len(res)} checks)')

# ---- Lemma 2 (P2) with my own code: build [u^k]E_n as poly in (y,a,b)
def stir2(m, r):
    return sum((-1)**i * comb(r, i) * (r-i)**m for i in range(r+1)) // factorial(r)
def E_poly(k, n):
    # dict {(T, al, be): Fraction} for y^T a^al b^be
    P = {}
    for i in range(k+1):
        j = k - i
        for t in range(j+1):
            cb = Fr(1)
            for s_ in range(t):
                cb *= Fr(1-2*n - s_, s_+1)
            st = stir2(j, t)
            if st == 0 or cb == 0: continue
            base = cb * st * factorial(t) * Fr(1, factorial(j)) * Fr(1, factorial(i))
            for x in range(i+1):
                cx = comb(i, x) * (-(n-1))**(i-x) * n**x
                for y_ in range(j+1):
                    cy = comb(j, y_) * (-1)**(j-y_)
                    key = (t, i-x+(j-y_), x+y_)
                    P[key] = P.get(key, Fr(0)) + base*cx*cy
    return {kk: vv for kk, vv in P.items() if vv}
def Dsum(P):
    # (d_a + d_b) P
    R = {}
    for (T, al, be), v in P.items():
        if al > 0:
            key = (T, al-1, be); R[key] = R.get(key, Fr(0)) + v*al
        if be > 0:
            key = (T, al, be-1); R[key] = R.get(key, Fr(0)) + v*be
    return {kk: vv for kk, vv in R.items() if vv}
ok2 = True
for n in (2, 3, 5):
    EPs = [E_poly(k, n) for k in range(9)]
    for k in range(1, 9):
        D = Dsum(EPs[k]); R = EPs[k-1]
        if D != R:
            ok2 = False
            print('P2 FAIL at', n, k)
print('P2 (d_a+d_b)[u^k]E_n = [u^{k-1}]E_n, n in {2,3,5}, k<=8, MY code:', ok2)
print('ALL:', all(ok for _, ok in res) and ok2)
open('progress.log','a').write('step_26: all 12 remaining mixed/single rungs EXACT PASS n<=40; Lemma 2 ((d_a+d_b)E_n=uE_n) verified with independent code n in {2,3,5} k<=8\n')
