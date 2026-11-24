import subprocess
r = subprocess.run(['bash','-c','''cd <project>/research/08_theory/RST1_rest_exponent && <project>/.venv/bin/python -c "
from flint import arb, ctx
ctx.dps=60
PI=arb.pi(); sq2pi=(2*PI).sqrt(); sqpi=PI.sqrt()
# max over N>=100 of h(N) = 1030 N^{3/2} e^{-(sq2pi-sqpi) sqrt N}
# h log-derivative: 3/(2N) - (sq2pi-sqpi)/(2 sqrt N) < 0 for sqrt N > 3/(sq2pi-sqpi)=4.086, i.e. N>16.7; so max at N=100
d = sq2pi-sqpi
print('d = sqrt(2pi)-sqrt(pi) =', d.str(18))
N=arb(100)
h100 = 1030*N**arb('1.5')*(-d*N.sqrt()).exp()
print('h(100) =', h100.str(18))
# check monotonic decrease on grid then value
import math
mx=max(1030*n**1.5*math.exp(-float(d)*math.sqrt(n)) for n in range(100,2000))
print('grid max (100..2000):', mx)
C2 = arb(10)**6
print('C3 = C2 + 1030*max =', (C2+h100).str(18))
"'''], capture_output=True, text=True, timeout=600)
print(r.stdout, r.stderr[-300:])
