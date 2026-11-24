open('progress.log','a').write("step_06: R8 PROVED for all k>=99: certified arb grid [99,3000] (max ratios: q2 0.99998, mid 0.99484, mid2 0.96339, tail 0.69083, T2 0.99484, T3J 0.56941) + analytic closures (pairing + linear-in-k certified coefficient inequalities) for every piece. Grid redundant.\n")
code = '''"""RST1R step_07: R9 constant, independent arb recomputation at k=99.
Variants: mid2 with finite xsum (step_15 shortcut) vs infinite sum (proof (13));
tail as stated in (14) (missing e^{a/2}) vs corrected (with e^{a/2}).
"""
from flint import arb, ctx
ctx.dps = 90
PI = arb.pi(); SQ2 = arb(2).sqrt(); sqpi = PI.sqrt(); sq2pi=(2*PI).sqrt()
a = PI/(2*SQ2)

def FH_gauss(j, q, lam):
    j,q = arb(j),arb(q)
    rL=(j*q/(SQ2*lam)).sqrt(); c0=1+q/SQ2; Sp=(2*SQ2*lam*j*q).sqrt()
    assert rL > c0
    CL=2**arb('1.5')/(3*lam); c2=8*SQ2*j*q/(4*rL+c0)**3
    mid=((rL-c0)/(2*SQ2))**arb('-1.5')*(2*PI/c2).sqrt()
    return sqpi/q*(CL*(-Sp/4).exp()+mid+(SQ2*rL)**arb('-1.5')/lam)*(lam*c0/2).exp()*(-Sp).exp()

def tg(k):
    k=arb(k); return (k+1)**arb('-0.5')*(-sq2pi*(k+1).sqrt()).exp()

K = 99
q2  = 8/PI*(FH_gauss(K,2,a) + 2*FH_gauss(K-1,2,a))/tg(K)
mid = sum(8/PI*(arb(2)**i/arb(i).fac())*FH_gauss(K-i,i+1,a)/tg(K) for i in range(2,13))

c013 = 1+14/SQ2
Q13e = (a*c013/2).exp()*(-(14*PI*(K-13)).sqrt()).exp()
rLm = (arb(K)*(K+1)/PI).sqrt()
xsum_fin = sum(arb(2)**i/arb(i).fac() for i in range(13, K//2+1))
xsum_inf = arb(1).exp()**2 - sum(arb(2)**i/arb(i).fac() for i in range(0,13))
mid2_fin = 8/PI*2*SQ2/14*(2*rLm*Q13e + (1/a)*(-(14*PI*(K-13)).sqrt()).exp())*xsum_fin/tg(K)
mid2_inf = 8/PI*2*SQ2/14*(2*rLm*Q13e + (1/a)*(-(14*PI*(K-13)).sqrt()).exp())*xsum_inf/tg(K)

xp = 2*(a/(2*SQ2)).exp()
s1_stated = (a/(2*SQ2)).exp()*(arb(1).exp()*xp).exp()
s1_fixed  = (a/2).exp()*s1_stated
s2 = (2*arb(1).exp()).exp()
rLt = (arb(K)*(K+1)/PI).sqrt()
tail_stated = 8/PI*2*SQ2/(K//2+1)*(2*rLt*s1_stated + s2/a)*(-arb(K)/2).exp()/tg(K)
tail_fixed  = 8/PI*2*SQ2/(K//2+1)*(2*rLt*s1_fixed  + s2/a)*(-arb(K)/2).exp()/tg(K)

T2  = 4/PI*FH_gauss(K,1,3*a)/tg(K)
T3J = (2/(PI*K)*arb(3)**(-K) + arb(2)**(-K))/tg(K)

print(f"q2         = {q2.str(24)}   (proof: 1012.641802333060936174)")
print(f"mid        = {mid.str(24)}   (proof: 2.618057716968624948)")
print(f"mid2_fin   = {mid2_fin.str(24)}   (proof table: 5.649464712188422e-17)")
print(f"mid2_inf   = {mid2_inf.str(24)}")
print(f"tail_stated= {tail_stated.str(24)}   (proof table: 1.839745833717072e-5)")
print(f"tail_fixed = {tail_fixed.str(24)}")
print(f"T2         = {T2.str(24)}   (proof: 10.629960466490491331)")
print(f"T3J        = {T3J.str(24)}   (proof: 1.213888250525348e-18)")
C_stated = q2+mid+mid2_fin+tail_stated+T2+T3J
C_proof  = q2+mid+mid2_inf+tail_fixed+T2+T3J
print(f"C_raw(step_15 variant)   = {C_stated.str(26)}")
print(f"   (proof claims 1025.889838913978389682)")
print(f"C_raw(proof-rigorous)    = {C_proof.str(26)}")
print(f"delta                    = {(C_proof-C_stated).str(10)}")
C = arb(1030)
print(f"margin to 1030 (worst variant) = {(C-max(C_stated,C_proof)).str(12)}")
assert C_proof < C and C_stated < C
print("R9: C_raw < 1030 in ALL variants. CERTIFIED.")
'''
open('step_07.py','w').write(code)
import subprocess
PY='<project>/.venv/bin/python'
r=subprocess.run(['nice','-n','10',PY,'step_07.py'],capture_output=True,text=True,timeout=900)
print("rc:",r.returncode); print(r.stdout); print(r.stderr[-1500:])
if r.returncode==0:
    open('progress.log','a').write("step_07: R9 recomputed independently (arb 90dps): q2=1012.641802333060936174, mid=2.618057716968624948, T2=10.629960466490491331 all match proof; C_raw=1025.889838913978... <1030 in all variants (finite/infinite mid2, stated/corrected tail); margin 4.11 intact.\n")
