import mpmath as mp, json
mp.mp.dps = 60
L = (mp.zeta(5, mp.mpf(1)/8) + mp.zeta(5, mp.mpf(3)/8) - mp.zeta(5, mp.mpf(5)/8) - mp.zeta(5, mp.mpf(7)/8))/8**5
print('L(chi8,5) mpmath:', mp.nstr(L, 40))
d=json.load(open('<project>/research/08_theory/M5b_w6_closed_form/Lchi8_5_1800_planner.json'))
print('planner           :', d['value'][:41])
# direct series partial check
s = mp.nsum(lambda n: (1 if n%8 in (1,3) else (-1 if n%8 in (5,7) else 0))/n**5, [1, 1000])
print('direct partial sum:', mp.nstr(s, 12))
