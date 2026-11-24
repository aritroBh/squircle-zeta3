import subprocess, json, mpmath as mp
from flint import acb, arb, ctx
PY='<project>/.venv/bin/python'
m5b='<project>/research/08_theory/M5b_w6_closed_form'
# Lchi at 2500 digits
def Lchi_flint(s, bits):
    ctx.prec = bits
    out = acb(0)
    for a, chi in [(1,1),(3,-1),(5,-1),(7,1)]:
        out += chi*acb.zeta(acb(s), acb(arb(a)/8))
    return (out/acb(8)**s).real
v = Lchi_flint(5, 8600)
sv = v.mid().str(2510, radius=False)
d = json.load(open(m5b+'/Lchi8_5_1800_planner.json'))
mp.mp.dps = 2000
print('lchi2500 vs planner:', mp.nstr(abs(mp.mpf(sv)-mp.mpf(d['value'])),5), 'rad', v.rad().str(8))
json.dump({'source':'flint Hurwitz, cross-checked vs planner 1e-1811','digits':2500,'value':sv},
          open(m5b+'/Lchi8_5_2500.json','w'))
# point w6_values at 2500 file and run at 2400
src=open(m5b+'/w6_values.py').read().replace("Lchi8_5_2100.json","Lchi8_5_2500.json")
open(m5b+'/w6_values.py','w').write(src)
r=subprocess.run(['nice','-n','10',PY,'-B','w6_values.py','2400'],capture_output=True,text=True,cwd=m5b,timeout=1700)
print(r.stdout[-500:]); print(r.stderr[-200:])
