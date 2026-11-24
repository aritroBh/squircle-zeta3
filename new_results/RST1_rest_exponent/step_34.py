import os
os.chdir('<project>/research/08_theory/RST1_rest_exponent')
step = r'''"""RST1 step_18b: certified lemma validation, saddle-aware quadrature."""
from flint import arb, acb, ctx
ctx.dps = 110
PI = arb.pi(); SQ2 = arb(2).sqrt(); sqpi = PI.sqrt()
a = PI/(2*SQ2)
W = acb(1,1)/SQ2

def FH_exact(j, q, lam):
    jj = arb(j); qq = arb(q)
    rL=(arb(j)*arb(q)/(SQ2*lam)).sqrt()
    def f(r, analytic):
        s = arb(1)/2 + r*W
        H = sqpi/4*(s/2).gamma()*(s/2+arb(3)/2).rgamma()
        return (-lam*r).exp()*abs(H)*abs(s)**jj/abs(s+qq)**(jj+1)
    value = acb(0)
    R = 20*rL+100
    pts=[0]
    for fr in ['0.25','0.5','0.75','1','1.5','2.5','5','10']:
        p = rL*arb(fr)
        if p < R: pts.append(p)
    pts.append(R)
    for lo,hi in zip(pts,pts[1:]):
        piece = acb.integral(f, lo, hi, abs_tol=arb("1e-130"), rel_tol=arb("1e-100"), eval_limit=200000)
        assert piece.is_finite(), (j,q,lo,hi,piece)
        value += piece
    xR = arb(1)/2 + R/SQ2
    tail = sqpi/SQ2 * xR**arb('-1.5')*SQ2/qq*(-lam*R).exp()/lam
    return value.real + arb(0, tail.upper())

def FH_gauss(j, q, lam):
    j,q = arb(j),arb(q)
    rL=(j*q/(SQ2*lam)).sqrt(); c0=1+q/SQ2; Sp=(2*SQ2*lam*j*q).sqrt()
    assert rL>c0
    CL=2**arb('1.5')/(3*lam); c2=8*SQ2*j*q/(4*rL+c0)**3
    mid=((rL-c0)/(2*SQ2))**arb('-1.5')*(2*PI/c2).sqrt()
    return sqpi/q*(CL*(-Sp/4).exp()+mid+(SQ2*rL)**arb('-1.5')/lam)*(lam*c0/2).exp()*(-Sp).exp()

def FH_box(j,q,lam):
    j,q=arb(j),arb(q)
    rL=(j*q/(SQ2*lam)).sqrt(); c0=1+q/SQ2; Sp=(2*SQ2*lam*j*q).sqrt()
    return 2*SQ2/q*(2*rL*(lam*c0/2).exp()+1/lam)*(-Sp).exp()

print("== Lemma R2/R3 validation (exact/bound) ==")
cases = [(99,2,a),(98,2,a),(400,2,a),(87,13,a),(97,3,a),(91,9,a),(99,1,3*a),(400,1,3*a),(50,5,a),(10,2,a),(200,6,a),(1,50,a)]
worst=arb(0)
for (j,q,lam) in cases:
    Fe=FH_exact(j,q,lam); Fg=FH_gauss(j,q,lam); Fb=FH_box(j,q,lam)
    rg=(Fe/Fg); rb=(Fe/Fb)
    worst=max(worst, rg.upper(), rb.upper())
    print(f" j={j} q={q} lam/a={(lam/a).str(4)}: exact={Fe.str(8)} gauss ratio={rg.str(6)} box ratio={rb.str(6)}", flush=True)
print(f"worst exact/bound ratio: {worst.str(8)} (<1 required)", flush=True)

print("== convexity lemma ==")
def Q(i,k):
    i,k=arb(i),arb(k)
    return (a*(i+1)/(2*SQ2)).exp()*(-(PI*(k-i)*(i+1)).sqrt()).exp()
bad=0
for k in [99,150,300,400]:
    Q2=Q(2,k)
    for i in range(2,k//2+1):
        if Q(i,k)>Q2: bad+=1
print(f"violations Q_i<=Q_2 (i in [2,k/2]): {bad}", flush=True)
bad2=0
for k in range(26,401):
    if Q(k//2,k) > Q(13,k): bad2+=1
print(f"violations Q_13>=Q_(k/2) for k in [26,400]: {bad2}", flush=True)
# h(k) = sqrt(pi(k/2)(k/2+1)) - sqrt(14 pi(k-13)) - a(k/2-12)/(2 sqrt2)  >=0 for k>=99 (Q_13 vs Q_k/2)
def h(k):
    k=arb(k)
    return (PI*(k/2)*(k/2+1)).sqrt()-(14*PI*(k-13)).sqrt()-a*(k/2-12)/(2*SQ2)
print("h(99) =", h(99).str(10), "; h increasing on grid:", all(h(k+1)>h(k) for k in range(99,300)), flush=True)
'''
open('step_18b.py','w').write(step)
import subprocess
r = subprocess.run(['bash','-c','cd <project>/research/08_theory/RST1_rest_exponent && nice <project>/.venv/bin/python step_18b.py'], capture_output=True, text=True, timeout=3000)
print(r.stdout)
print(r.stderr[-1500:] if r.stderr else '')
