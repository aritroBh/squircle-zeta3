import subprocess
r = subprocess.run(['bash','-c','''cd <project>/research/08_theory/RST1_rest_exponent && <project>/.venv/bin/python -c "
from flint import arb, ctx
ctx.dps=50
PI=arb.pi(); SQ2=arb(2).sqrt(); a=PI/(2*SQ2)
def Q(i,k):
    i,k=arb(i),arb(k)
    return (a*(i+1)/(2*SQ2)).exp()*(-(PI*(k-i)*(i+1)).sqrt()).exp()
viol=[k for k in range(26,401) if Q(k//2,k)>Q(13,k)]
print('violation k-range:', viol[0], '..', viol[-1], 'count', len(viol))
print('any violation for k>=99?', [k for k in viol if k>=99])
"'''], capture_output=True, text=True, timeout=300)
print(r.stdout, r.stderr[-300:])
