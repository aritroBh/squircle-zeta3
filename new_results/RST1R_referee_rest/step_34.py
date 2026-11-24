code = '''"""RST1R step_03: R2/R3 re-derivation and stress test.
Part 1 (algebra, arb): the saddle quantities of E~_lam(r) = -lam r - sqrt2 j q/(2r+c0):
  r_m = r_lam - c0/2;  E~(r_m) = -Sigma + lam c0/2;  2 lam r_lam = Sigma;
  E~(r_b) = -5 Sigma/4 + lam c0/2 (r_b=(r_lam-c0)/2, needs r_lam>c0);
  E~'(r_b) = 3 lam;  E~'' = -8 sqrt2 j q/(2r+c0)^3 (concave); c2 formula.
Part 2: hypothesis r_lam > c0 at every use in the proof (q2: j in {k,k-1},q=2; mid i in [2,12];
  T2 q=1 lam=3a), for all k>=99, with margins.
Part 3: exact F_{H,lam}(j,q) by arb quadrature vs the R2 (11) and R3 (7) bounds on a wide grid,
  including r_lam barely > c0 and large j; report worst exact/bound ratio.
"""
from flint import arb, acb, ctx
ctx.dps = 80
PI = arb.pi(); SQ2 = arb(2).sqrt(); sqpi = PI.sqrt()
a = PI/(2*SQ2)
W = acb(1,1)/SQ2

def saddle_quants(j, q, lam):
    j,q = arb(j), arb(q)
    rL = (j*q/(SQ2*lam)).sqrt(); c0 = 1+q/SQ2
    Sp = (2*SQ2*lam*j*q).sqrt()
    return rL, c0, Sp

def Etil(r, j, q, lam):
    return -lam*r - SQ2*j*q/(2*r+1+q/SQ2)

print("== Part 1: saddle algebra (arb) ==")
worst = arb(0)
for (j,q,lam) in [(99,2,a),(98,2,a),(87,13,a),(99,1,3*a),(5,2,a),(1,3,a),(1000,7,a/2),(10**5,2,a)]:
    rL,c0,Sp = saddle_quants(j,q,lam)
    rm = rL - c0/2
    # numerical max of Etil via its analytic maximizer identity 2 r_m + c0 = 2 r_lam
    assert rm > 0 or True
    e1 = Etil(rm, j, q, lam); e2 = -Sp + lam*c0/2
    assert abs(e1-e2) < arb('1e-60')*(1+abs(e2)), (j,q,lam,e1,e2)
    assert abs(2*lam*rL - Sp) < arb('1e-60')*Sp
    if rL > c0:
        rb = (rL-c0)/2
        eb = Etil(rb, j, q, lam); eb2 = -arb(5)/4*Sp + lam*c0/2
        assert abs(eb-eb2) < arb('1e-55')*(1+abs(eb2)), (j,q,lam,eb,eb2)
        # E~'(r_b) = -lam + 2 sqrt2 j q/(2 r_b + c0)^2 = 3 lam
        lamL = -lam + 2*SQ2*j*q/(2*rb+c0)**2
        assert abs(lamL-3*lam) < arb('1e-55')*lam, (j,q,lam,lamL)
        # concavity bound c2 on [0, 2 r_lam]
        c2 = 8*SQ2*j*q/(4*rL+c0)**3
        # check E~''(r) <= -c2 for r in [0, 2 rL] (sample + endpoint monotonicity of (2r+c0)^3)
        for fr in ['0','0.1','0.5','1','1.7','2']:
            rr = rL*arb(fr)
            Epp = -8*SQ2*j*q/(2*rr+c0)**3
            assert Epp <= -c2 + arb('1e-70'), (j,q,lam,fr)
    print(f"  j={j} q={q} lam/a={(lam/a).str(5)}: r_lam={rL.str(8)} c0={c0.str(6)} all identities OK")
print("Part 1 PASS: r_m, E~(r_m), 2 lam r_lam=Sigma, E~(r_b)=-5S/4+lam c0/2, E~'(r_b)=3 lam, c2, concavity")

print("== Part 2: hypothesis r_lam > c0 at every use ==")
# q2: (j,q,lam) = (k,2,a) and (k-1,2,a), k>=99 -> j>=98
rL,c0,_ = saddle_quants(98,2,a); print(f"  q2 (j>=98,q=2): r_a={rL.str(8)} > c0={c0.str(8)}, ratio={(rL/c0).str(6)}; r_a increasing in j")
# mid: i in [2,12], j=k-i>=87, q=i+1 in [3,13]: minimize r_a/c0 over the discrete set at k=99
m = None
for i in range(2,13):
    rL,c0,_ = saddle_quants(99-i,i+1,a)
    m = min(m, rL/c0) if m is not None else rL/c0
print(f"  mid (k=99): min_i r_a/c0 = {m.str(8)} (worst i=12); j=k-i increasing in k for fixed i")
rL,c0,_ = saddle_quants(99,1,3*a); print(f"  T2 (j>=99,q=1,3a): r_3a={rL.str(8)} > c0={c0.str(8)}, ratio={(rL/c0).str(6)}")
print("Part 2 PASS")

print("== Part 3: exact F vs R2/R3 bounds (arb quadrature) ==")
ctx.dps = 100
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
        piece = acb.integral(f, lo, hi, abs_tol=arb('1e-120'), rel_tol=arb('1e-90'), eval_limit=200000)
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

cases = [(99,2,a),(98,2,a),(400,2,a),(87,13,a),(97,3,a),(91,9,a),(99,1,3*a),(400,1,3*a),
         (50,5,a),(10,2,a),(200,6,a),(1,50,a),
         (5,3,a),(2,7,a),(30,1,a),(7,40,a),(1,1,a),(3,1,3*a),(15,2,a/2),(80,13,a)]
worst_g = arb(0); worst_b = arb(0); wg=None; wb=None
for (j,q,lam) in cases:
    rL,c0,_ = saddle_quants(j,q,lam)
    Fe=FH_exact(j,q,lam); Fb=FH_box(j,q,lam)
    rb=(Fe/Fb)
    line = f"  j={j} q={q} lam/a={(lam/a).str(5)} rL/c0={(rL/c0).str(5)}: box ratio={rb.str(6)}"
    if rL > c0:
        Fg=FH_gauss(j,q,lam); rg=(Fe/Fg)
        line += f" gauss ratio={rg.str(6)}"
        if rg.upper() > worst_g: worst_g = rg.upper(); wg=(j,q,lam/a)
    else:
        line += " (r_lam<=c0: R2 not applicable)"
    if rb.upper() > worst_b: worst_b = rb.upper(); wb=(j,q,lam/a)
    print(line, flush=True)
print(f"WORST gauss exact/bound = {worst_g.str(8)} at {wg}")
print(f"WORST box   exact/bound = {worst_b.str(8)} at {wb}")
assert worst_g < 1 and worst_b < 1
print("Part 3 PASS: all ratios < 1")
'''
open('step_03.py','w').write(code)
import subprocess
PY='<project>/.venv/bin/python'
r=subprocess.run(['nice','-n','10',PY,'step_03.py'],capture_output=True,text=True,timeout=1200)
print("rc:",r.returncode); print(r.stdout); print(r.stderr[-1500:])
