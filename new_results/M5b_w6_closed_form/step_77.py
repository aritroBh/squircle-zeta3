import sys, json
from flint import arb, acb, ctx
ctx.prec = 300
# chi8: primitive mod-8 character with chi(1)=chi(3)=1, chi(5)=chi(7)=-1
def Lchi(s, prec):
    ctx.prec = prec
    out = acb(0)
    for a in [1,3,5,7]:
        chi = 1 if a in (1,3) else -1
        out += chi*acb(a/8).zeta(acb(s), acb(a,8))
    return (out/acb(8)**s).real
v = Lchi(5, 300)
print(v.str(40))
d = json.load(open('<project>/research/08_theory/M5b_w6_closed_form/Lchi8_5.json'))
print(d['value'][:40])
import mpmath as mp
mp.mp.dps = 60
print('diff:', mp.nstr(abs(mp.mpf(v.mid().str(50)) - mp.mpf(d['value'][:55])), 5))
