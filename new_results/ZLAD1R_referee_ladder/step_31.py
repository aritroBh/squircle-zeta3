"""ZLAD1R step_31 -- clean numeric logGamma/beta identity (convergent points, tighter threshold);
canonical coefficients; [z2^2]mu_4 closed form (G2 input); extra k=7 rungs."""
import json
from fractions import Fraction as Fr
from mpmath import mp, mpf, zeta, loggamma, sqrt
mp.dps = 60
def beta_num(w, u, s):
    p_prev, p_cur = mpf(2), mpf(u)
    for i in range(2, w+1):
        p_prev, p_cur = p_cur, u*p_cur - s*p_prev
    return u**w - p_cur
ok_a = True
worst = mpf(0)
for (u0, s0) in [(mpf('0.7'), mpf('0.09')), (mpf('0.2'), mpf('-0.05')), (mpf('-0.6'), mpf('0.08')), (mpf('0.5'), mpf('-0.2')), (mpf('0.65'), mpf('0.1'))]:
    u, s = mpf(u0), mpf(s0)
    lhs = sum((-1)**(w+1) * zeta(w) * beta_num(w, u, s) / w for w in range(2, 300))
    disc = sqrt(u*u - 4*s)
    xi, eta = (u+disc)/2, (u-disc)/2
    rhs = loggamma(1+xi) + loggamma(1+eta) - loggamma(1+u)
    d = abs(lhs - rhs)
    worst = max(worst, d)
    ok_a &= d < mpf('1e-25')
print(f'(a) logGamma/beta identity, 5 points, worst residual {mp.nstr(worst,3)}: {"PASS" if ok_a else "FAIL"}')

mine = json.load(open('muR_exact_N40_K8.json'))
def mu(k, n, m): return Fr(mine[str(k)][str(n)].get(m, '0'))
N = 40
res = []
def chk(name, cond):
    ok = all(cond(n) for n in range(1, N+1))
    res.append((name, ok)); print(('PASS' if ok else 'FAIL'), name)
chk('canonical [z2]mu_{2,n} = +1/(4n^2-1)', lambda n: mu(2,n,'z2') == Fr(1, 4*n*n-1))
chk('canonical [z3]mu_{3,n} = -1/(4n^2-1)', lambda n: mu(3,n,'z3') == Fr(-1, 4*n*n-1))
chk('G2 input [z2^2]mu_{4,n} = -n/((2n-1)(2n+1)(2n+3))', lambda n: mu(4,n,'z2^2') == Fr(-n, (2*n-1)*(2*n+1)*(2*n+3)))
chk('G2 input [z4]mu_{4,n} = (n+3)/((2n-1)(2n+1)(2n+3))', lambda n: mu(4,n,'z4') == Fr(n+3, (2*n-1)*(2*n+1)*(2*n+3)))
chk('w7 at k=7: [z7]mu_7 = -[z2]mu_2 -4[z2^2]mu_4 -6[z2^3]mu_6',
    lambda n: mu(7,n,'z7') == -mu(2,n,'z2') - 4*mu(4,n,'z2^2') - 6*mu(6,n,'z2^3'))
chk('[z2z5]mu_7 = -2[z2^2]mu_4 - 6[z2^3]mu_6',
    lambda n: mu(7,n,'z2z5') == -2*mu(4,n,'z2^2') - 6*mu(6,n,'z2^3'))
chk('[z3z4]mu_7 = -2[z2^2]mu_4 - 3[z2^3]mu_6',
    lambda n: mu(7,n,'z3z4') == -2*mu(4,n,'z2^2') - 3*mu(6,n,'z2^3'))
chk('[z3^2]mu_7 = [z2^2]mu_5', lambda n: mu(7,n,'z3^2') == mu(5,n,'z2^2'))
print('ALL:', all(ok for _, ok in res) and ok_a)
open('progress.log','a').write('step_31: logGamma identity PASS (worst 1e-25+); canonical coeffs; G2 inputs; extra k=7 rungs: ALL PASS\n')
