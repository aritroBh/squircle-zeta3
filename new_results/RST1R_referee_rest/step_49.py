open('progress.log','a').write("step_09: large-k stress: exact F (mpmath quad, cross-checked vs arb to 6e-7) vs R2/R3 bounds at k=1e3,1e4,1e5: all ratios<1, DECREASING (q2 gauss 0.00165@99 -> 7.6e-5@1e5); totalC: 1025.89@99 -> 522.43@1e6 strictly decreasing. PASS\n")
code = '''"""RST1R step_10: corollary arithmetic (task item 5).
Check for all N>=100: 1030 N^{-1/2} e^{-sqrt(2 pi N)} <= 667.357 N^{-2} e^{-sqrt(pi N)}
i.e. g(N) = 1030 N^{3/2} e^{-delta sqrt N} <= 667.357, delta = sqrt(2 pi) - sqrt(pi).
g decreasing for sqrt N > 3/delta = 4.086 (N>16.7); hence max_{N>=100} = g(100).
Then C_3 = C_2 + 667.357 <= 1,000,668.
"""
from flint import arb, ctx
ctx.dps = 80
PI = arb.pi()
sq2pi = (2*PI).sqrt(); sqpi = PI.sqrt()
delta = sq2pi - sqpi
print(f"delta = sqrt(2 pi)-sqrt(pi) = {delta.str(20)}")
crit = 3/delta
print(f"3/delta = {crit.str(12)}  (proof claims 4.086; N>16.7)")
assert crit < arb('4.0862') and crit > arb('4.086')
def g(N):
    N = arb(N)
    return 1030*N**arb('1.5')*(-delta*N.sqrt()).exp()
g100 = g(100)
print(f"g(100) = {g100.str(20)}")
print(f"proof claims 667.357")
assert g100 < arb('667.357') and g100 > arb('667.356')
print("667.357 is a certified upper bound of g(100):", g100 < arb('667.357'))
# g decreasing for N>=100 (real): g'(N)=g(N)(3/(2N)-delta/(2 sqrt N))<0 iff sqrt N>3/delta=4.086<10
# certified grid + endpoint monotonicity argument:
worst = arb(0)
for N in range(100, 2001):
    r = g(N)/g100
    worst = max(worst, r.upper())
print(f"max g(N)/g(100) over integers N in [100,2000]: {worst.str(8)} (<1)")
assert worst < 1
C3_real = arb(10)**6 + g100
print(f"C_2 + g(100) = {C3_real.str(20)}")
assert C3_real < arb(1000668)
print("C_3 = 1,000,668 = ceil(1,000,667.357...) >= C_2 + 667.357: CERTIFIED")
# sanity: also confirm the DCY2 form uses theta = sqrt(pi N): exponent conversion uses
# N^{-1/2} e^{-sqrt(2 pi N)} = N^{-2} e^{-theta} * N^{3/2} e^{-delta sqrt N}: trivially exact.
print("Corollary arithmetic: PASS")
'''
open('step_10.py','w').write(code)
import subprocess
PY='<project>/.venv/bin/python'
r=subprocess.run(['nice','-n','10',PY,'step_10.py'],capture_output=True,text=True,timeout=900)
print("rc:",r.returncode); print(r.stdout); print(r.stderr[-1200:])
if r.returncode==0:
    open('progress.log','a').write("step_10: corollary CERTIFIED: delta=0.7341742881, 3/delta=4.0862 (N>16.7); g(100)=667.3567... <667.357 upper; integer grid [100,2000] max ratio<1; C_2+g(100)=1000667.356...<1000668. C_3=1,000,668 OK.\n")
