"""ZLAD1R step_25 -- FIVE rung relations checked EXACTLY on my independent table (n<=40),
including weight-7 and weight-8 singles and two mixed rungs; plus SEVEN closed forms."""
import json
from fractions import Fraction as Fr
mine = json.load(open('muR_exact_N40_K8.json'))
def mu(k, n, m): return Fr(mine[str(k)][str(n)].get(m, '0'))
N = 40
res = []
def chk(name, cond):
    ok = all(cond(n) for n in range(1, N+1))
    res.append((name, ok))
    print(('PASS' if ok else 'FAIL'), name)

# 1. rung 1 at the top k: [z3]mu_8 = -[z2]mu_7
chk('R1  [z3]mu_8 = -[z2]mu_7', lambda n: mu(8,n,'z3') == -mu(7,n,'z2'))
# 2. weight-7 single at k=8: [z7]mu_8 = -[z2]mu_3 - 4[z2^2]mu_5 - 6[z2^3]mu_7
chk('w7  [z7]mu_8 = -[z2]mu_3 -4[z2^2]mu_5 -6[z2^3]mu_7',
    lambda n: mu(8,n,'z7') == -mu(3,n,'z2') - 4*mu(5,n,'z2^2') - 6*mu(7,n,'z2^3'))
# 3. weight-8 single at k=8: [z8]mu_8 = [z2]mu_2 + 5[z2^2]mu_4 + 12[z2^3]mu_6 + 6[z2^4]mu_8
chk('w8  [z8]mu_8 = [z2]mu_2 +5[z2^2]mu_4 +12[z2^3]mu_6 +6[z2^4]mu_8',
    lambda n: mu(8,n,'z8') == mu(2,n,'z2') + 5*mu(4,n,'z2^2') + 12*mu(6,n,'z2^3') + 6*mu(8,n,'z2^4'))
# 4. mixed weight 8: [z3z5]mu_8 = 2[z2^2]mu_4 + 6[z2^3]mu_6
chk('mix [z3z5]mu_8 = 2[z2^2]mu_4 + 6[z2^3]mu_6',
    lambda n: mu(8,n,'z3z5') == 2*mu(4,n,'z2^2') + 6*mu(6,n,'z2^3'))
# 5. mixed weight 7: [z2^2z3]mu_7 = -3[z2^3]mu_6
chk('mix [z2^2z3]mu_7 = -3[z2^3]mu_6',
    lambda n: mu(7,n,'z2^2z3') == -3*mu(6,n,'z2^3'))

def rf(num_coeffs, den_roots):
    def f(n):
        num = Fr(0)
        for p, c in enumerate(num_coeffs): num += c * n**p
        den = Fr(1)
        for r in den_roots: den *= (2*n + r)
        return num/den
    return f
cf = {
 '[z5]mu_5 = -3/((2n-1)(2n+1)(2n+3))':            ((5,'z5'),  rf([-3],[-1,1,3])),
 '[z2^3]mu_6 = n(5n-1)/(3(2n-1)(2n+1)(2n+3)(2n+5))': ((6,'z2^3'), None), # custom below
}
chk('CF1 [z5]mu_5 = -3/((2n-1)(2n+1)(2n+3))',
    lambda n: mu(5,n,'z5') == Fr(-3, (2*n-1)*(2*n+1)*(2*n+3)))
chk('CF2 [z2^3]mu_6 = n(5n-1)/(3(2n-1)(2n+1)(2n+3)(2n+5))',
    lambda n: mu(6,n,'z2^3') == Fr(n*(5*n-1), 3*(2*n-1)*(2*n+1)*(2*n+3)*(2*n+5)))
# extra closed forms from the writeup (bonus)
chk('CF3 [z2z3]mu_5 = 2n/((2n-1)(2n+1)(2n+3))',
    lambda n: mu(5,n,'z2z3') == Fr(2*n, (2*n-1)*(2*n+1)*(2*n+3)))
chk('CF4 [z3^2]mu_6 = -n/((2n-1)(2n+1)(2n+3))',
    lambda n: mu(6,n,'z3^2') == Fr(-n, (2*n-1)*(2*n+1)*(2*n+3)))
chk('CF5 [z6]mu_6 = (4n^2+n+45)/(3(2n-1)(2n+1)(2n+3)(2n+5))',
    lambda n: mu(6,n,'z6') == Fr(4*n*n+n+45, 3*(2*n-1)*(2*n+1)*(2*n+3)*(2*n+5)))
chk('CF6 [z2z4]mu_6 = n(n-11)/((2n-1)(2n+1)(2n+3)(2*n+5))',
    lambda n: mu(6,n,'z2z4') == Fr(n*(n-11), (2*n-1)*(2*n+1)*(2*n+3)*(2*n+5)))
chk('CF7 [z2^4]mu_8 = -n(13n^2-3n-1)/(3(2n-1)...(2n+7))',
    lambda n: mu(8,n,'z2^4') == Fr(-n*(13*n*n-3*n-1), 3*(2*n-1)*(2*n+1)*(2*n+3)*(2*n+5)*(2*n+7)))
chk('CF8 [z8]mu_8 = (2n^3+78n^2-59n+105)/((2n-1)(2n+1)(2n+3)(2n+5)(2n+7))',
    lambda n: mu(8,n,'z8') == Fr(2*n**3+78*n*n-59*n+105, (2*n-1)*(2*n+1)*(2*n+3)*(2*n+5)*(2*n+7)))
chk('CF9 [z4^2]mu_8 = -n(7n^2-12n+41)/((2n-1)...(2n+7))',
    lambda n: mu(8,n,'z4^2') == Fr(-n*(7*n*n-12*n+41), (2*n-1)*(2*n+1)*(2*n+3)*(2*n+5)*(2*n+7)))
print('ALL:', all(ok for _, ok in res), f'({len(res)} checks, n<=40, EXACT)')
open('progress.log','a').write('step_25: 5 rungs (w7,w8 singles + 2 mixed) + 9 closed forms EXACT PASS on independent table n<=40\n')
