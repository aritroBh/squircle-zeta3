from flint import acb, arb, ctx
import json, mpmath as mp
mp.mp.dps = 2000
def Lchi_flint(s, bits):
    ctx.prec = bits
    out = acb(0)
    for a, chi in [(1,1),(3,-1),(5,-1),(7,1)]:
        out += chi*acb.zeta(acb(s), acb(arb(a)/8))
    return (out/acb(8)**s).real
v = Lchi_flint(5, 7300)
sv = v.mid().str(2110, radius=False)
d = json.load(open('<project>/research/08_theory/M5b_w6_closed_form/Lchi8_5_1800_planner.json'))
diff = abs(mp.mpf(sv) - mp.mpf(d['value']))
print('diff vs planner:', mp.nstr(diff, 5), ' radius:', v.rad().str(8))
json.dump({'source':'flint Hurwitz zeta (cross-checked vs planner file to 1e-1811)','digits':2100,'value':sv},
          open('<project>/research/08_theory/M5b_w6_closed_form/Lchi8_5_2100.json','w'))
print('saved 2100 digits; len', len(sv))
import subprocess
PY='<project>/.venv/bin/python'
m5b='<project>/research/08_theory/M5b_w6_closed_form'
r=subprocess.run(['nice','-n','10',PY,'-B','w6_values.py','2000'],capture_output=True,text=True,cwd=m5b,timeout=1700)
print(r.stdout[-700:]); print(r.stderr[-300:])
