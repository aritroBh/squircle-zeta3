"""RST1 step_03: numerical anatomy of the exact rest M_k - M_main.
- exact M_k from D13f_sweep/mk_deg9.tsv (arb balls, midpoints at dps=90)
- M_main = (2/pi) Im int_ray H E Q_k ds  (mpmath quadrature, CX6/DCY2 formula)
- rest = M_k - M_main: magnitude vs e^{-sqrt(2 pi N)}, phase, and the q=2 model
  R_k = -(4/pi) Im int_Gamma H E * 2 s^{k-1}(s+1)/(s+2)^{k+1} ds
"""
import mpmath as mp, re
mp.mp.dps = 90
pi=mp.pi; j=1j; omega=mp.expjpi(mp.mpf(1)/4)

def load():
    out={}
    for line in open("<project>/research/08_theory/D13f_sweep/mk_deg9.tsv"):
        p=line.rstrip("\n").split("\t")
        if not p[0].isdigit(): continue
        mid=re.match(r"\[([-0-9.eE+]+)",p[1]).group(1)
        out[int(p[0])]=mp.mpf(mid)
    return out
rows=load()

def H(s): return mp.beta(s/2, mp.mpf(3)/2)/2

def Mmain(k):
    N=k+1
    rs=mp.sqrt(2*N/pi)
    f=lambda r: mp.expj(pi*(r*omega)/2)*(r*omega/(1+r*omega))**N*H(r*omega)/(r*omega)*omega
    pts=[0,rs/8,rs/2,rs,2*rs,4*rs,10*rs,mp.inf]
    return 2/pi*mp.quad(f,pts).imag

def Rmodel(k):
    """q=2 model of the rest: -(4/pi) Im int_ray H E s^{k-1}(s+1)/(s+2)^{k+1} ds"""
    rs=2*mp.sqrt(k/pi)
    def f(r):
        s=r*omega
        return H(s)*mp.expj(pi*s/2)*s**(k-1)*(s+1)/(s+2)**(k+1)*omega
    pts=[0,rs/8,rs/2,rs,2*rs,4*rs,10*rs,mp.inf]
    return -4/pi*mp.quad(f,pts).imag

sq2pi=mp.sqrt(2*pi)
print(" k     N   |rest|          log|rest|+sq(2piN)   |rest|*N e^{sq(2piN)}   rest/R_model   phase(rest)-sq(2piN) (mod 2pi)")
for k in [100,120,140,160,180,200,220,240,260,280,300,320,340,360,380,400]:
    N=k+1
    mm=Mmain(k); rest=rows[k]-mm; Rm=Rmodel(k)
    lg=mp.log(abs(rest))+sq2pi*mp.sqrt(N)
    ph=mp.fmod(mp.arg(rest)-sq2pi*mp.sqrt(N),2*pi)
    print(f"{k:4d} {N:4d}  {mp.nstr(abs(rest),6):12s}  {mp.nstr(lg,8):12s}   {mp.nstr(abs(rest)*N*mp.exp(sq2pi*mp.sqrt(N)),6):10s}   {mp.nstr(rest/Rm,5):10s}   {mp.nstr(ph,6)}")
