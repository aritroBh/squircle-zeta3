
import os, subprocess
WORK = 'research/08_theory/D13s_fix_recertify/D13s_work'
FIX = os.path.join(WORK, 'fixed')
code = r'''
# H3 lower-box: certify H'(t) = t^3 Hp3(t) >= 0 on [0, 0.13] via arb box eval of the
# manifestly regular Hp3 (no cancellation); then H(t) >= H(0) = 0 => truth >= h3-box lower.
from flint import arb, ctx
ctx.dps = 60
rt2 = arb(2).sqrt()
def iv(lo, hi):
    lo = arb(lo).lower(); hi = arb(hi).upper()
    return arb((lo+hi)/2, (hi-lo)/2)
def Hp3(Tb):
    den = 1+rt2*Tb+Tb**2
    q = ((1+rt2*Tb)/den).sqrt()
    w = 1/(den*(1+q))
    Gp = (rt2/3)*(rt2+2*Tb)/den**2
    hp = rt2*w**2*Tb + rt2*w*Tb*(2+rt2*Tb)/(q*den**2)
    return Gp - hp
# single box over [0, 0.13] plus finer subdivisions to display margin
for (a,b) in [(0.0, 0.13), (0.0, 0.01), (0.01, 0.05), (0.05, 0.13)]:
    v = Hp3(iv(a,b))
    print(f"Hp3([{a},{b}]) = {v}  lower = {arb(v.lower())}")
v = Hp3(iv(0.0, 0.13))
assert arb(v.lower()) > 0, "Hp3 not certified positive"
print("CERTIFIED: Hp3(t) > 0 on [0, 0.13]  =>  H increasing  =>  H >= 0  (h3 lower box contains truth for t <= 0.13;")
print("certified domain N>=100 has t <= y(1+y^2) at y=sqrt(pi/200):", end=" ")
y = (arb.pi()/200).sqrt(); print(arb((y*(1+y*y)).upper()), "< 0.128 < 0.13 OK)")
'''
open(os.path.join(FIX,'certify_H3.py'),'w').write(code)
PY='<home>)/.venv/bin/python'
r = subprocess.run([PY, 'certify_H3.py'], cwd=FIX, capture_output=True, text=True, timeout=600)
print(r.stdout); print(r.stderr[-1000:] if r.stderr else "")
