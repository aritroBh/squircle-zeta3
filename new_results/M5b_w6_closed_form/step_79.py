from flint import arb, acb, ctx
import json, mpmath as mp
def Lchi(s, prec):
    ctx.prec = prec
    out = acb(0)
    for a in [1,3,5,7]:
        chi = 1 if a in (1,3) else -1
        out += chi*acb.zeta(acb(s), acb(a,8))
    return (out/acb(8)**s).real
v = Lchi(5, 400)
print(v.str(40))
d = json.load(open('<project>/research/08_theory/M5b_w6_closed_form/Lchi8_5.json'))
mp.mp.dps = 80
print('diff vs M4e:', mp.nstr(abs(mp.mpf(v.mid().str(70)) - mp.mpf(d['value'][:75])), 5))
