import os, subprocess
os.chdir("<project>/research/08_theory/DCY2_decay_first_correction")
partb = '''"""step_07_partb.py — part (b): |Model - T1_complex| <= C_b N^{-2} e^{-theta}.

Model = (2/pi) e^{phi(s0)} H(s0) omega sum_p p_p mu_p^inf,
T1    = -i (sqrt pi/N) e^{-theta+i theta} (1+varpi1/sqrt N).
Write  Model = -i y e^{-theta+i theta} * ModelRatio,  T1 = -i y e^{-theta+i theta} * T1ratio:
  ModelRatio = (2/pi) e^{S(y)} H(s0) omega sum p_p mu_p,
  T1ratio    = (2y/sqrt pi)(1+varpi1/sqrt N).
Then |Model - T1| * N^2 e^theta = y N^2 |ModelRatio - T1ratio| = (pi^2/(4 y^3)) |ModelRatio - T1ratio|.
No e^{-theta} factors anywhere: interval-N friendly.
"""
from flint import arb, acb, ctx

PI = arb.pi(); SQ2 = arb(2).sqrt()
OMEGA = acb(1,1)/SQ2; WB = OMEGA.conjugate()
bC = PI*WB/2
varpi1 = acb( (24+12*PI+PI**2)/(48*PI.sqrt()), (24-12*PI+PI**2)/(48*PI.sqrt()) )

def y_of(N): return (PI/(2*arb(N))).sqrt()

def S_of(N):
    y = y_of(N)
    return sum((-1)**m * PI*y**(m-2)*WB**m/(2*m) for m in range(3, 12))

def H_of(N):
    s0 = OMEGA/y_of(N)
    return PI.sqrt()/2 * (s0/2+1).gamma() * (s0/2+arb(3)/2).rgamma()

def ell_of(N, m):
    y = y_of(N); N = arb(N)
    if m == 1: return (N-1)*y - N*y/(1+y*WB) + acb(0,1)*PI*OMEGA/2
    return (-1)**(m+1) * y**m/m * ((N-1) - N/(1+y*WB)**m)

def mu_p_of(N, p):
    y = y_of(N)
    l1, l2 = ell_of(N,1), ell_of(N,2)
    t0 = -l1/(2*l2)
    base = (-l2*t0**2).exp()*(-PI/l2).sqrt()
    tot = acb(0)
    for m in range(0, p//2+1):
        df = arb(1)
        for k in range(1, 2*m+1, 2):
            df *= k
        tot += arb(arb(p).fac()/(arb(2*m).fac()*arb(p-2*m).fac())) * t0**(p-2*m) * df / (-2*l2)**m
    return base*tot

def momentsum_of(N):
    y = y_of(N)
    l3,l4,l5 = ell_of(N,3), ell_of(N,4), ell_of(N,5)
    s0 = OMEGA/y
    H0 = H_of(N)
    # h_j via polygamma: H^{(m)}/(m! H) * omega^m, H'/H = (1/2)(psi(s/2)-psi(3/2+s/2))
    # H^{(m)} = sum binom(m-1,k) L_{k+1} H^{(m-1-k)}, L_j = 2^{-(j+1)}(psi^{(j-1)}(s/2)-psi^{(j-1)}(3/2+s/2))
    Ls = []
    for j in range(1, 6):
        Ls.append( arb(2)**(-(j)) * ((s0/2).digamma() if j==1 else (s0/2).polygamma(j-1))
                   - arb(2)**(-(j)) * ((s0/2+arb(3)/2).digamma() if j==1 else (s0/2+arb(3)/2).polygamma(j-1)) )
    Hd = [H0]
    for m in range(1, 4):
        Hd.append(sum(arb(arb(m-1).fac()/(arb(k).fac()*arb(m-1-k).fac()))*Ls[k]*Hd[m-1-k] for k in range(m)))
    h1 = OMEGA*Hd[1]/H0
    h2 = OMEGA**2*Hd[2]/(2*H0)
    h3 = OMEGA**3*Hd[3]/(6*H0)
    P = {}
    def add(p,c): P[p] = P.get(p, acb(0))+c
    add(0,1); add(1,h1); add(2,h2); add(3,h3)
    add(3,l3); add(4,l4); add(5,l5)
    add(6,l3**2/2); add(7,l3*l4); add(9,l3**3/6)
    add(4,h1*l3); add(5,h1*l4); add(5,h2*l3); add(7,h1*l3**2/2)
    return sum(cp*mu_p_of(N,p) for p,cp in P.items())

def ModelRatio_of(N):
    return (2/PI)*S_of(N).exp()*H_of(N)*OMEGA*momentsum_of(N)

def T1ratio_of(N):
    return (2*y_of(N)/PI.sqrt())*(1+varpi1/arb(N).sqrt())

def D_of(N):
    y = y_of(N)
    return (PI**2/(4*y**3))*abs(ModelRatio_of(N) - T1ratio_of(N))

if __name__ == "__main__":
    ctx.dps = 40
    for N in [400, 1000, 2000, 3000, 10000, 100000]:
        print(f"N={N:7d}  D(N) = {D_of(N).str(10)}")
    # interval cells
    for lo, hi in [(400,408),(1000,1020),(3000,3060),(10000,10200),(100000,102000)]:
        N = arb((lo+hi)/2, (hi-lo)/2)
        d = D_of(N)
        print(f"[{lo},{hi}]  D upper = {d.upper().str(10)}")
'''
open("step_07_partb.py","w").write(partb)
import subprocess
r = subprocess.run(["<project>/.venv/bin/python","step_07_partb.py"],
                   capture_output=True, text=True, cwd=".", timeout=600)
print(r.stdout, r.stderr[-2500:] if r.returncode else "")
