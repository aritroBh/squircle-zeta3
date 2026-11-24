from flint import acb, arb, ctx
import json, mpmath as mp
mp.mp.dps = 1900
def Lchi_flint(s, bits):
    ctx.prec = bits
    out = acb(0)
    for a, chi in [(1,1),(3,-1),(5,-1),(7,1)]:
        out += chi*acb.zeta(acb(s), acb(arb(a)/8))
    return (out/acb(8)**s).real
v = Lchi_flint(5, 6600)
print(v.str(20))
sv = v.mid().str(1960, radius=False)
d = json.load(open('<project>/research/08_theory/M5b_w6_closed_form/Lchi8_5_1800_planner.json'))
diff = abs(mp.mpf(sv) - mp.mpf(d['value']))
print('flint Hurwitz vs planner 1810 digits: diff =', mp.nstr(diff, 5))
print('flint radius:', v.rad().str(8))
assert v.imag.contains(0)
if diff < mp.mpf('1e-1805'):
    json.dump({'source':'flint Hurwitz zeta (cross-checked vs planner file to 1e-1805)','digits':1950,'value':sv},
              open('<project>/research/08_theory/M5b_w6_closed_form/Lchi8_5_2100.json','w'))
    print('saved Lchi8_5_2100.json (1950 digits)')
