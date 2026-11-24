import json, time
import mpmath as mp
from flint import acb, arb, ctx
t0=time.monotonic()
ctx.prec = 12600
out = acb(0)
for a, chi in [(1,1),(3,-1),(5,-1),(7,1)]:
    out += chi*acb.zeta(acb(5), acb(arb(a)/8))
v = (out/acb(8)**5).real
print('time %.1fs rad %s'%(time.monotonic()-t0, v.rad().str(8)))
sv = v.mid().str(3620, radius=False)
# cross-check vs 2500-digit file
old = json.load(open('Lchi8_5_2500.json'))['value']
mp.mp.dps = 2600
diff = abs(mp.mpf(sv) - mp.mpf(old))
print('Lchi3600 vs Lchi8_5_2500 absdiff:', mp.nstr(diff, 8))
json.dump({'source':'flint Hurwitz zeta at 12600 bits (chi: +1 on n=+-1 mod 8); cross-checked vs Lchi8_5_2500','digits':3600,'value':sv}, open('Lchi8_5_3600.json','w'))
print('saved Lchi8_5_3600.json; first 60 digits:', sv[:60])
print(open('<project>/research/08_theory/M5c_w6_high_precision/w6_hp.out').read())
