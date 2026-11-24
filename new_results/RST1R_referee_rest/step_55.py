code = open('step_11.py').read()
code = code.replace("""s = sp.symbols('s'); k = sp.symbols('k', integer=True, positive=True)
expr1 = -s**k/(s+2)**(k+1) - s**(k-1)/(s+2)**k + 2*s**(k-1)*(s+1)/(s+2)**(k+1)
d1 = sp.powsimp(sp.expand(sp.together(expr1)))
print("(a1) A_k^[2] identity difference:", d1)
assert d1 == 0
# also spot-check at several integer k,s
import random
random.seed(1)
for kk in [1,2,7,50]:
    for sv in [sp.Rational(3,7), sp.Integer(5), sp.I]:
        assert sp.simplify(expr1.subs({k:kk, s:sv})) == 0
print("(a1) also spot-checked at integer k and rational/imaginary s: 0")""",
"""s = sp.symbols('s'); k = sp.symbols('k', integer=True, positive=True)
# numerator after clearing (s+2)^{k+1}: -s^k - s^{k-1}(s+2) + 2 s^{k-1}(s+1)
num = sp.expand(-s**k - s**(k-1)*(s+2) + 2*s**(k-1)*(s+1))
print("(a1) cleared numerator:", num)
assert num == 0   # -s^k - s^k - 2 s^{k-1} + 2 s^k + 2 s^{k-1} = 0
# numeric spot-check at many integer k and complex s (mpmath, 80 dps)
import mpmath as _mp
_mp.mp.dps = 80
import random
random.seed(1)
worst = _mp.mpf(0)
for _ in range(50):
    kk = random.randint(1, 500)
    sv = _mp.mpc(random.uniform(-5,5), random.uniform(-5,5))
    lhs = -sv**kk/(sv+2)**(kk+1) - sv**(kk-1)/(sv+2)**kk
    rhs = -2*sv**(kk-1)*(sv+1)/(sv+2)**(kk+1)
    worst = max(worst, abs(lhs-rhs)/max(1,abs(rhs)))
print(f"(a1) mpmath spot-check 50 random (k,s): max rel diff = {_mp.nstr(worst,3)}")
assert worst < _mp.mpf('1e-60')""")
open('step_11.py','w').write(code)
import subprocess
PY='<project>/.venv/bin/python'
r=subprocess.run(['nice','-n','10',PY,'step_11.py'],capture_output=True,text=True,timeout=3000)
print("rc:",r.returncode); print(r.stdout); print(r.stderr[-1500:])
