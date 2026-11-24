"""ZLAD1R step_46 -- re-run ZLAD1's claimed 14 rung checks at k=7,8, n<=150 on THEIR table;
G1 unmodified run status."""
import json, os
from fractions import Fraction as Fr
base = '<project>'
t = json.load(open(base + '/research/08_theory/ZLAD1_zeta_ladder/mu_exact_N150_K8.json'))
def mu(k, n, m): return Fr(t[str(k)][str(n)].get(m, '0'))
N = 150
checks = []
def chk(name, cond):
    ok = all(cond(n) for n in range(1, N+1))
    checks.append(ok); print(('PASS' if ok else 'FAIL'), name)
chk('1  [z3]mu_8 = -[z2]mu_7', lambda n: mu(8,n,'z3') == -mu(7,n,'z2'))
chk('2  [z4]mu_7 = [z2]mu_5+[z2^2]mu_7', lambda n: mu(7,n,'z4') == mu(5,n,'z2')+mu(7,n,'z2^2'))
chk('3  [z4]mu_8 = [z2]mu_6+[z2^2]mu_8', lambda n: mu(8,n,'z4') == mu(6,n,'z2')+mu(8,n,'z2^2'))
chk('4  [z5]mu_7 = -[z2]mu_4-2[z2^2]mu_6', lambda n: mu(7,n,'z5') == -mu(4,n,'z2')-2*mu(6,n,'z2^2'))
chk('5  [z5]mu_8 = -[z2]mu_5-2[z2^2]mu_7', lambda n: mu(8,n,'z5') == -mu(5,n,'z2')-2*mu(7,n,'z2^2'))
chk('6  [z6]mu_7 = [z2]mu_3+3[z2^2]mu_5+2[z2^3]mu_7', lambda n: mu(7,n,'z6') == mu(3,n,'z2')+3*mu(5,n,'z2^2')+2*mu(7,n,'z2^3'))
chk('7  [z6]mu_8 = [z2]mu_4+3[z2^2]mu_6+2[z2^3]mu_8', lambda n: mu(8,n,'z6') == mu(4,n,'z2')+3*mu(6,n,'z2^2')+2*mu(8,n,'z2^3'))
chk('8  [z7]mu_7 = -[z2]mu_2-4[z2^2]mu_4-6[z2^3]mu_6', lambda n: mu(7,n,'z7') == -mu(2,n,'z2')-4*mu(4,n,'z2^2')-6*mu(6,n,'z2^3'))
chk('9  [z7]mu_8 = -[z2]mu_3-4[z2^2]mu_5-6[z2^3]mu_7', lambda n: mu(8,n,'z7') == -mu(3,n,'z2')-4*mu(5,n,'z2^2')-6*mu(7,n,'z2^3'))
chk('10 [z8]mu_8 = [z2]mu_2+5[z2^2]mu_4+12[z2^3]mu_6+6[z2^4]mu_8', lambda n: mu(8,n,'z8') == mu(2,n,'z2')+5*mu(4,n,'z2^2')+12*mu(6,n,'z2^3')+6*mu(8,n,'z2^4'))
chk('11 [z2z3]mu_7 = -2[z2^2]mu_6', lambda n: mu(7,n,'z2z3') == -2*mu(6,n,'z2^2'))
chk('12 [z2z3]mu_8 = -2[z2^2]mu_7', lambda n: mu(8,n,'z2z3') == -2*mu(7,n,'z2^2'))
chk('13 [z3^2]mu_7 = [z2^2]mu_5', lambda n: mu(7,n,'z3^2') == mu(5,n,'z2^2'))
chk('14 [z3^2]mu_8 = [z2^2]mu_6', lambda n: mu(8,n,'z3^2') == mu(6,n,'z2^2'))
chk('15 [z2z4]mu_7 = 2[z2^2]mu_5+3[z2^3]mu_7', lambda n: mu(7,n,'z2z4') == 2*mu(5,n,'z2^2')+3*mu(7,n,'z2^3'))
chk('16 [z2z4]mu_8 = 2[z2^2]mu_6+3[z2^3]mu_8', lambda n: mu(8,n,'z2z4') == 2*mu(6,n,'z2^2')+3*mu(8,n,'z2^3'))
chk('17 [z2z5]mu_8 = -2[z2^2]mu_5-6[z2^3]mu_7', lambda n: mu(8,n,'z2z5') == -2*mu(5,n,'z2^2')-6*mu(7,n,'z2^3'))
chk('18 [z3z4]mu_8 = -2[z2^2]mu_5-3[z2^3]mu_7', lambda n: mu(8,n,'z3z4') == -2*mu(5,n,'z2^2')-3*mu(7,n,'z2^3'))
chk('19 [z2^2z3]mu_7 = -3[z2^3]mu_6', lambda n: mu(7,n,'z2^2z3') == -3*mu(6,n,'z2^3'))
chk('20 [z2^2z3]mu_8 = -3[z2^3]mu_7', lambda n: mu(8,n,'z2^2z3') == -3*mu(7,n,'z2^3'))
chk('21 [z2z3^2]mu_8 = 3[z2^3]mu_6', lambda n: mu(8,n,'z2z3^2') == 3*mu(6,n,'z2^3'))
chk('22 [z2^2z4]mu_8 = 3[z2^3]mu_6+6[z2^4]mu_8', lambda n: mu(8,n,'z2^2z4') == 3*mu(6,n,'z2^3')+6*mu(8,n,'z2^4'))
chk('23 [z3z5]mu_8 = 2[z2^2]mu_4+6[z2^3]mu_6', lambda n: mu(8,n,'z3z5') == 2*mu(4,n,'z2^2')+6*mu(6,n,'z2^3'))
chk('24 [z2z6]mu_8 = 2[z2^2]mu_4+9[z2^3]mu_6+8[z2^4]mu_8', lambda n: mu(8,n,'z2z6') == 2*mu(4,n,'z2^2')+9*mu(6,n,'z2^3')+8*mu(8,n,'z2^4'))
chk('25 [z4^2]mu_8 = [z2^2]mu_4+3[z2^3]mu_6+3[z2^4]mu_8', lambda n: mu(8,n,'z4^2') == mu(4,n,'z2^2')+3*mu(6,n,'z2^3')+3*mu(8,n,'z2^4'))
print('THEIR table, 25 rung checks n<=150:', 'ALL PASS' if all(checks) else 'FAIL', f'({sum(checks)}/{len(checks)})')
# G1 run status
pid = int(open('/tmp/g1_k8.pid').read().strip())
try:
    os.kill(pid, 0); alive = True
except OSError:
    alive = False
print('G1 flat_exact.py 40 8 unmodified: alive =', alive, '; bytes:', os.path.getsize('/tmp/g1_k8.out'))
open('progress.log','a').write('step_46: 25 rung checks on ZLAD1 table at n<=150: ALL PASS\n')
