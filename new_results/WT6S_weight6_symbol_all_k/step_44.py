
import os
WD = '<project>/research/08_theory/WT6S_weight6_symbol_all_k'
s = open(os.path.join(WD, 'step_05.py')).read()
s = s.replace("""        poly = [sum(sp.binomial(kv - 6, m - 1) * epv[m][t] for m in range(1, 6))
                for t in range(4)]""",
"""        poly = [sum((epv[m][t] * int(sp.binomial(kv - 6, m - 1)) for m in range(1, 6)), ZEROQ)
                for t in range(4)]""")
s = s.replace("""        expr = sum(sp.binomial(k - 6, m - 1) * (epv[m][t].p + sp.sqrt(2) * epv[m][t].q)
                   for m in range(1, 6))
        expr = sp.expand(expr)
        rat = sp.expand(expr).subs(sp.sqrt(2), 0)
        sqr = sp.Rational(str(sp.expand(expr).coeff(sp.sqrt(2), 1)))
        log(f"pairing[{labels[t]}](k) = ({rat}) + ({sqr})*sqrt2")""",
"""        expr = sp.expand(sum(sp.binomial(k - 6, m - 1) * epv[m][t].sympy()
                             for m in range(1, 6)))
        rat = sp.expand(expr.subs(sp.sqrt(2), 0))
        sqr = sp.expand(expr - rat) / sp.sqrt(2)
        sqr = sp.expand(sqr)
        log(f"pairing[{labels[t]}](k) = ({rat}) + ({sqr})*sqrt2")""")
open(os.path.join(WD, 'step_05.py'), 'w').write(s)
import subprocess
PY = '<project>/.venv/bin/python'
r = subprocess.run(['nice', '-n', '10', PY, '-B', 'step_05.py'], capture_output=True, text=True, cwd=WD)
print(r.stdout)
print("STDERR:", r.stderr[-2000:])
