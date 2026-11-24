"""RST1R step_10: corollary arithmetic (task item 5).
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
assert crit < arb('4.0863') and crit > arb('4.086')
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
for N in range(101, 2001):
    r = g(N)/g100
    worst = max(worst, r.upper())
print(f"max g(N)/g(100) over integers N in (100,2000]: {worst.str(8)} (<1)")
assert worst < 1
C3_real = arb(10)**6 + g100
print(f"C_2 + g(100) = {C3_real.str(20)}")
assert C3_real < arb(1000668)
print("C_3 = 1,000,668 = ceil(1,000,667.357...) >= C_2 + 667.357: CERTIFIED")
# sanity: also confirm the DCY2 form uses theta = sqrt(pi N): exponent conversion uses
# N^{-1/2} e^{-sqrt(2 pi N)} = N^{-2} e^{-theta} * N^{3/2} e^{-delta sqrt N}: trivially exact.
print("Corollary arithmetic: PASS")
