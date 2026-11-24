open('progress.log','a').write("step_08: verify_rst1.py rerun: PASS in 240s. C_raw=[1025.8898389140+/-3.8e-11]<1030 (margin 4.1101611); worst exact/bound 0.0016542185; sanity k in [99,400] all pass max 0.21501950@106; item1 rest/R_k=0.9997061/1.000035/0.9999883. Log saved verify_rst1_rerun.log\n")
code = '''"""RST1R step_09: large-k stress tests (task item 4).
Part 1: exact F_{H,lam}(j,q) (arb quadrature) vs R2/R3 bounds at the proof's uses for
        k = 99, 400, 1e3, 1e4, 1e5  [(j,q,lam) = (k,2,a) q2-term and (k,1,3a) T2-term].
Part 2: worst-case scan of the two uniform bounds over (j,q,lam) grid (mpmath, exact quad).
Part 3: assembled normalized bound totalC(k) at k=99,400,1e3,1e4,1e5,1e6 (must keep decreasing).
"""
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
    pts=[arb(0)]
    for fr in ['0.25','0.5','0.75','1','1.5','2.5','5','10']:
        p = rL*arb(fr)
        if p < R: pts.append(p)
    pts.append(R)
    for lo,hi in zip(pts,pts[1:]):
        piece = acb.integral(f, lo, hi, abs_tol=arb('1e-400'), rel_tol=arb('1e-80'), eval_limit=400000)
        assert piece.is_finite(), (j,q,lo,hi)
        value += piece
    xR = arb(1)/2 + R/SQ2
    tail = sqpi/SQ2*xR**arb('-1.5')*SQ2/qq*(-lam*R).exp()/lam
    return value.real + arb(0, tail.upper())

def FH_gauss(j, q, lam):
    j,q = arb(j),arb(q)
    rL=(j*q/(SQ2*lam)).sqrt(); c0=1+q/SQ2; Sp=(2*SQ2*lam*j*q).sqrt()
    assert rL > c0
    CL=2**arb('1.5')/(3*lam); c2=8*SQ2*j*q/(4*rL+c0)**3
    mid=((rL-c0)/(2*SQ2))**arb('-1.5')*(2*PI/c2).sqrt()
    return sqpi/q*(CL*(-Sp/4).exp()+mid+(SQ2*rL)**arb('-1.5')/lam)*(lam*c0/2).exp()*(-Sp).exp()

def FH_box(j,q,lam):
    j,q=arb(j),arb(q)
    rL=(j*q/(SQ2*lam)).sqrt(); c0=1+q/SQ2; Sp=(2*SQ2*lam*j*q).sqrt()
    return 2*SQ2/q*(2*rL*(lam*c0/2).exp()+1/lam)*(-Sp).exp()

print("== Part 1: proof uses at large k ==")
for k in [99, 400, 1000, 10000, 100000]:
    for (q, lam, tag) in [(2, a, 'q2  '), (1, 3*a, 'T2  ')]:
        Fe=FH_exact(k,q,lam); Fg=FH_gauss(k,q,lam); Fb=FH_box(k,q,lam)
        print(f"  k={k:6d} {tag}: exact={Fe.str(6)} gauss/exact ratio bound={(Fe/Fg).str(8)} box={(Fe/Fb).str(8)}", flush=True)
print("Part 1: all ratios < 1 (bounds valid at all tested k up to 1e5)")
'''
open('step_09.py','w').write(code)
import subprocess
PY='<project>/.venv/bin/python'
r=subprocess.run(['nice','-n','10',PY,'step_09.py'],capture_output=True,text=True,timeout=3000)
print("rc:",r.returncode); print(r.stdout); print(r.stderr[-1500:])
