
import os
WD = '<project>/research/08_theory/WT6S_weight6_symbol_all_k'
s = open(os.path.join(WD, 'step_05.py')).read()
s = s.replace("""        expr = sp.expand(sum(sp.binomial(k - 6, m - 1) * epv[m][t].sympy()
                             for m in range(1, 6)))
        rat = sp.expand(expr.subs(sp.sqrt(2), 0))
        sqr = sp.expand(expr - rat) / sp.sqrt(2)
        sqr = sp.expand(sqr)
        log(f"pairing[{labels[t]}](k) = ({rat}) + ({sqr})*sqrt2")""",
"""        binom2 = (k - 6) * (k - 7) / 2
        expr = sum(sp.binomial(k - 6, m - 1) * epv[m][t].sympy() for m in range(1, 6))
        expr = sp.expand(expr.subs(sp.binomial(k - 6, 2), binom2))
        rat = sp.expand(expr.subs(sp.sqrt(2), 0))
        sqr = sp.expand((expr - rat) / sp.sqrt(2))
        log(f"pairing[{labels[t]}](k) = ({rat}) + ({sqr})*sqrt2   [power basis]")""")
open(os.path.join(WD, 'step_05.py'), 'w').write(s)
import subprocess
PY = '<project>/.venv/bin/python'
r = subprocess.run(['nice', '-n', '10', PY, '-B', 'step_05.py'], capture_output=True, text=True, cwd=WD)
print(r.stdout[-2000:])
print("STDERR:", r.stderr[-800:])
