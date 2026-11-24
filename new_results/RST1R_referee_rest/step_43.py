code = open('step_09.py').read()
old = code[code.find('def FH_exact'):code.find('def FH_gauss')]
new = '''def FH_exact(j, q, lam):
    jj = arb(j); qq = arb(q)
    rL=(arb(j)*arb(q)/(SQ2*lam)).sqrt(); c0=1+qq/SQ2; Sp=(2*SQ2*lam*jj*qq).sqrt()
    est = FH_box(j,q,lam)
    atol = est*arb('1e-90')
    def f(r, analytic):
        s = arb(1)/2 + r*W
        H = sqpi/4*(s/2).gamma()*(s/2+arb(3)/2).rgamma()
        return (-lam*r).exp()*abs(H)*abs(s)**jj/abs(s+qq)**(jj+1)
    def Etil(r): return -lam*r - SQ2*jj*qq/(2*r+c0)
    R = 20*rL+100
    pts=[arb(0)]
    for fr in ['0.25','0.5','0.75','1','1.5','2.5','5','10']:
        p = rL*arb(fr)
        if p < R: pts.append(p)
    pts.append(R)
    # skip left segments whose certified max is < e^{-150} e^{-Sp}; replace by allowance
    from flint import arb as _arb
    cutoff = (-Sp - 150)
    value = acb(0); allowance = arb(0)
    started = False
    for lo,hi in zip(pts,pts[1:]):
        if not started and Etil(hi) < cutoff and hi < rL:
            # certified allowance: (hi-lo)*(sqrt2/q)*sqrt(pi/2)*2^{3/2}*e^{Etil(hi)}
            allowance += (hi-lo)*(SQ2/qq)*sqpi/SQ2*2**arb('1.5')*Etil(hi).exp()
            continue
        started = True
        piece = acb.integral(f, lo, hi, abs_tol=atol, rel_tol=arb('1e-70'), eval_limit=400000)
        assert piece.is_finite(), (j,q,lo,hi)
        value += piece
    xR = arb(1)/2 + R/SQ2
    tail = sqpi/SQ2*xR**arb('-1.5')*SQ2/qq*(-lam*R).exp()/lam
    return value.real + allowance + arb(0, tail.upper())

'''
code = code.replace(old, new)
open('step_09.py','w').write(code)
import subprocess
PY='<project>/.venv/bin/python'
r=subprocess.run(['nice','-n','10',PY,'step_09.py'],capture_output=True,text=True,timeout=4000)
print("rc:",r.returncode); print(r.stdout); print(r.stderr[-1200:])
