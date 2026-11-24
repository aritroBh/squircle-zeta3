import pickle, mpmath as mp, time
mp.mp.dps=150
out=pickle.load(open('gamma_pt1_w12.pkl','rb'))
c12=out[12]
CLASS=c12[9]+mp.mpf(14)/9*c12[10]
print("W12 irreducible coords at Gamma point (a1=0.01):")
print(f"  z2^2 z(3,5)   = {mp.nstr(c12[7],15)}")
print(f"  z2 z(3,7)     = {mp.nstr(c12[8],15)}")
print(f"  CLASS z39+14/9 z57 = {mp.nstr(CLASS,15)}")
print(f"  z(2,1,1,8)    = {mp.nstr(c12[11],15)}")
m=max(abs(c12[7]),abs(c12[8]),abs(CLASS),abs(c12[11]))
print(f"  max|irreducible coord| = {mp.nstr(m,10)}  (scale max|coord|={mp.nstr(max(abs(x) for x in c12),4)})")
def log_progress(line):
    with open('progress.log','a') as f: f.write(time.strftime('%H:%M:%S ')+line+'\n')
log_progress("Stage4 W12: traced Gamma at a1=0.01 (Q11 crosses 0 at a2~2.7413). Gamma point (a1,a2,a4,a5)=(0.01,2.74131,2.74174,1.04607). W12 mpf engine (150d): controls W8 z(3,5),W10 z(3,7),W11 z3z35 & f11 all vanish ~1e-72 (confirms on Gamma). W12 irreducibles ALL NONZERO: z2^2z35=-46.72, z2z37=-250.29, CLASS(z39+14/9 z57)=+919.27, z(2,1,1,8)=+389.34. => this Gamma point NOT irreducible-free. Answer trending NO.")
print("recorded")
