import sympy as sp, pickle, os, mpmath as mp, time
wd = "<project>/research/08_theory/DEC2_decagon_exact_assembly"
os.chdir(wd)
mp.mp.dps = 100
s5 = mp.sqrt(5); phi = (1+s5)/2
x1 = -1/phi; x4 = 1/phi**2
Jtab = pickle.load(open('Jtab_dec2.pkl','rb'))
L1,L2,L3,AA,z2,xs = sp.symbols('L1 L2 L3 AA z2 x')
def atom_vals(xv):
    return {L1: -mp.log(1-xv), L2: mp.polylog(2,xv), L3: mp.polylog(3,xv),
            AA: mp.nsum(lambda m: xv**m*mp.harmonic(m)/m**2, [1, mp.inf]),
            z2: mp.pi**2/6, xs: xv}
def quad(tp, kap):
    n,a,b,c = tp
    return mp.quad(lambda tt: tt**n*mp.log(tt)**a*mp.log(1-tt)**b*mp.log(1+kap*tt)**c, [0,1])
t0=time.time(); worst1=worst4=mp.mpf(0)
for side,(xv,kap) in [(1,(x1,1/phi)), (4,(x4,-1/phi**2))]:
    av = atom_vals(xv)
    for tp,e in Jtab.items():
        val = mp.mpf(sp.N(e.subs(av, simultaneous=True), 100))
        q = quad(tp, kap)
        d = abs(val-q)
        if side==1: worst1=max(worst1,d)
        else: worst4=max(worst4,d)
print("Jtab at x1: worst |table-quad| =", mp.nstr(worst1,3))
print("Jtab at x4: worst |table-quad| =", mp.nstr(worst4,3), " (%.0fs, 80 values)"%(time.time()-t0))
with open('progress.log','a') as f:
    f.write(f"step_08: Jtab validated at x1,x4 vs quadrature (100 dps): worst {mp.nstr(worst1,2)}, {mp.nstr(worst4,2)}\n")
