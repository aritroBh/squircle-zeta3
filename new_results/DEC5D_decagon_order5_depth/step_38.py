
import pickle, mpmath as mp
mp.mp.dps = 200
r120 = pickle.load(open("p5_run120.pkl","rb"))['res']
r170 = pickle.load(open("p5_run170.pkl","rb"))['res']
v120 = [mp.mpf(x) for x in r120]; v170 = [mp.mpf(x) for x in r170]
for k in range(6):
    d = abs(v120[k]-v170[k])
    print("P^(%d): |dps120 - dps170| = %s  (agreement ~%d digits)" % (k, mp.nstr(d,3), int(-mp.log10(d)) if d>0 else 999))

sqrt5 = mp.sqrt(5); phi = (1+sqrt5)/2
L = mp.sqrt((25-11*sqrt5)/10)
print("\nP(1) - 10L =", mp.nstr(abs(v170[0]-10*L), 3))
# P''' vs 60 L S at high precision
ft = pickle.load(open("<project>/research/08_theory/DEC1_decagon/final_table_pq.pkl","rb"))['table_pq']
l2 = mp.log(2); lf = mp.log(phi); z2 = mp.pi**2/6; z3 = mp.zeta(3)
Li3p = mp.polylog(3, 1/phi); Li3m = mp.polylog(3, -1/phi)
vals = {'l2^3': l2**3, 'l2^2*lf': l2**2*lf, 'l2^2': l2**2, 'l2*lf^2': l2*lf**2,
 'l2*lf': l2*lf, 'l2*z2': l2*z2, 'l2': l2, 'lf^3': lf**3, 'lf^2': lf**2,
 'lf*z2': lf*z2, 'lf': lf, 'z2': z2, 'z3': z3, 'L3(-1/phi)': Li3m, 'L3(1/phi)': Li3p}
S = sum((mp.mpf(str(p)) + mp.mpf(str(q))*sqrt5)*vals[k] for k,(p,q) in ft.items())
print("P''' - 60LS =", mp.nstr(abs(v170[3]-60*L*S), 3))
# also check C_5 zeta(3) coefficient claim (thm:Cm): C_m = 48 m sin^4(pi/2m) cos(pi/2m) (3 sin^2(pi/2m)-1)
C5 = 48*5*mp.sin(mp.pi/10)**4*mp.cos(mp.pi/10)*(3*mp.sin(mp.pi/10)**2-1)
print("\nC_5 (thm:Cm) =", mp.nstr(C5, 40))
C5b = mp.sqrt(10-2*sqrt5)*(105*sqrt5-255)/32
print("C_5 (DEC1 form) =", mp.nstr(C5b, 40), " diff:", mp.nstr(abs(C5-C5b),3))
# save P5/L for later
P5overL = v170[5]/L
print("\nP5/L =", mp.nstr(P5overL, 120))
pickle.dump({'P5overL': mp.nstr(P5overL, 180), 'P5': mp.nstr(v170[5],180), 'Pvals': [mp.nstr(v,180) for v in v170]}, open("p5_gold.pkl","wb"))
