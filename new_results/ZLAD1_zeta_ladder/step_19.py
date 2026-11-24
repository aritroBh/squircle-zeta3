import json
from fractions import Fraction as Fr
t = json.load(open('mu_exact_N8_K6.json'))
# canonical checks
ok2 = all(Fr(t['2'][str(n)].get('z2','0'))*(4*n*n-1)==1 for n in range(1,9))
ok3 = all(Fr(t['3'][str(n)].get('z3','0'))*(4*n*n-1)==-1 for n in range(1,9))
print('z2@k2:', ok2, ' z3@k3:', ok3)
# compare with G1 flat_exact N100 K3
ref = json.load(open('<project>/research/10_frontier/G1_zeta3_any_angle/flat_exact_N100_K3.json'))
keymap = {'(0, 0, 0)':'1','(1, 0, 0)':'z2','(0, 1, 0)':'z3','(0, 0, 1)':'z4','(2, 0, 0)':'z2^2'}
mismatch = 0
for k in ['0','1','2','3']:
    for n in range(1,9):
        for mk, mv in ref[k][str(n)].items():
            mine = t[k][str(n)].get(keymap[mk],'0')
            if Fr(mine) != Fr(mv):
                mismatch += 1
                print('MISMATCH', k, n, mk, mv, mine)
print('k<=3 vs G1 flat_exact: mismatches =', mismatch)
# G2 K4 json
ref4 = json.load(open('<project>/research/08_theory/G2_M4alpha/flat_exact_N100_K4.json'))
mm = 0
for n in range(1,9):
    for mk, mv in ref4['4'][str(n)].items():
        mine = t['4'][str(n)].get(keymap[mk],'0')
        if Fr(mine) != Fr(mv):
            mm += 1; print('MISMATCH4', n, mk, mv, mine)
print('k=4 vs G2: mismatches =', mm)
for n in [1,2,3]:
    print('mu_5,%s:'%n, t['5'][str(n)])
    print('mu_6,%s:'%n, t['6'][str(n)])
