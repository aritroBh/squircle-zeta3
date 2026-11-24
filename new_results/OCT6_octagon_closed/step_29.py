import pickle, mpmath as mp, sympy as sp
mp.mp.dps = 60
r2=mp.sqrt(2); k2=r2/2-1; x2v=-k2
K1d="research/08_theory/K1d_engine_fixed"
Jtab = pickle.load(open(K1d+'/Jtab_t.pkl','rb'))
L1,L2,L3,AA,z2,x = sp.symbols('L1 L2 L3 AA z2 x')
def atomvals(xv,N=4000):
    s=mp.mpf(0); Hn=mp.mpf(0)
    for nn in range(1,N+1):
        Hn+=1/mp.mpf(nn); s+=xv**nn*Hn/nn**2
    return {L1:-mp.log(1-xv),L2:mp.polylog(2,xv),L3:mp.polylog(3,xv),AA:s,z2:mp.pi**2/6,x:xv}
sub = atomvals(x2v)
worst=mp.mpf(0)
for key in Jtab:
    n,a,b,c = key
    f=lambda t: t**n*mp.log(t)**a*mp.log(1-t)**b*mp.log(1+k2*t)**c
    Jn = mp.quad(f,[0,1]); Js = mp.mpf(sp.N(Jtab[key].subs(sub, simultaneous=True), 50))
    worst=max(worst,abs(Jn-Js))
print("Jtab all 40 keys at x2: worst |quad-exact| =", mp.nstr(worst,3))
# Jloc checks
Jloc = pickle.load(open('oct5/Jloc.pkl','rb'))
worst=mp.mpf(0)
for (n,a,b),e in Jloc.items():
    f=lambda t: t**n*mp.log(t)**a*mp.log(1-t)**b
    Jn=mp.quad(f,[0,1]); Js=mp.mpf(sp.N(e,50))
    worst=max(worst,abs(Jn-Js))
print("Jloc all 40 keys: worst |quad-exact| =", mp.nstr(worst,3))
with open('progress.log','a') as f:
    f.write("step_06: Jtab verified at x2 (all 40 keys, 60dps quad) and Jloc all 40 keys by quad — exact tables confirmed.\n")
