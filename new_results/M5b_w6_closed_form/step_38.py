import os
m5b='<project>/research/08_theory/M5b_w6_closed_form'
src=open(os.path.join(m5b,'m5b_precision.py')).read()
old="""    tail_b, tail_a = arb(0), arb(0)
    for i, j, ts in T:
        pb = sum((mag(c/(vb-poles[z])) for z, c in ts), arb(0))
        pa = sum((mag(c/(va-poles[z])) for z, c in ts), arb(0))
        tail_b += pb*moment(Cb, i)*Mb**j
        tail_a += pa*La**i*moment(Ca, j)"""
new="""    def pow_int(x, k):  # flint quirk: 0-centered ball ** int = nan; use repeated multiply
        out = arb(1)
        for _ in range(k):
            out *= x
        return out
    tail_b, tail_a = arb(0), arb(0)
    for i, j, ts in T:
        pb = sum((mag(c/(vb-poles[z])) for z, c in ts), arb(0))
        pa = sum((mag(c/(va-poles[z])) for z, c in ts), arb(0))
        tail_b += pb*moment(Cb, i)*pow_int(Mb, j)
        tail_a += pa*pow_int(La, i)*moment(Ca, j)"""
assert old in src
src=src.replace(old,new)
open(os.path.join(m5b,'m5b_precision.py'),'w').write(src)
print('fixed')
# rerun quick test (checkpoint exists for 600)
import subprocess
PY='<project>/.venv/bin/python'
r=subprocess.run(['nice','-n','10',PY,'-B','m5b_precision.py','w6','600'],capture_output=True,text=True,cwd=m5b,timeout=900)
print(r.stdout[-1200:]); print(r.stderr[-800:])
