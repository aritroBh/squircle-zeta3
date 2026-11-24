import os
os.chdir('<project>/research/08_theory/RST1_rest_exponent')
import subprocess
# finish convexity check from step_18b (was cut off)
r = subprocess.run(['bash','-c','''cd <project>/research/08_theory/RST1_rest_exponent && nice <project>/.venv/bin/python -c "
from flint import arb, ctx
ctx.dps=50
PI=arb.pi(); SQ2=arb(2).sqrt(); a=PI/(2*SQ2)
def Q(i,k):
    i,k=arb(i),arb(k)
    return (a*(i+1)/(2*SQ2)).exp()*(-(PI*(k-i)*(i+1)).sqrt()).exp()
bad=0
for k in [99,150,300,400]:
    Q2=Q(2,k)
    for i in range(2,k//2+1):
        if Q(i,k)>Q2: bad+=1
print('violations Q_i<=Q_2 (i in [2,k/2]):', bad)
bad2=0
for k in range(26,401):
    if Q(k//2,k) > Q(13,k): bad2+=1
print('violations Q_13>=Q_(k/2) for k in [26,400]:', bad2)
def h(k):
    k=arb(k)
    return (PI*(k/2)*(k/2+1)).sqrt()-(14*PI*(k-13)).sqrt()-a*(k/2-12)/(2*SQ2)
print('h(99) =', h(99).str(10), '; h increasing on [99,300):', all(h(k+1)>h(k) for k in range(99,300)))
"'''], capture_output=True, text=True, timeout=1200)
print(r.stdout, r.stderr[-500:])
