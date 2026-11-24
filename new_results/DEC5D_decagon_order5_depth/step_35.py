
import mpmath as mp
mp.mp.dps = 60
sqrt5 = mp.sqrt(5)
phi = (1+sqrt5)/2
L = mp.sqrt((25-11*sqrt5)/10)
print("10L        =", mp.nstr(10*L, 55))
print("P(1)       = 2.008114158862272798697976726327421423306")
P1 = mp.mpf("2.008114158862272798697976726327421423306")
print("diff P(1)-10L:", mp.nstr(abs(P1-10*L), 4))

# 15-term S from DEC1 final_table_pq.pkl
import pickle
base = "<project>/research/08_theory/"
ft = pickle.load(open(base+"DEC1_decagon/final_table_pq.pkl","rb"))['table_pq']
l2 = mp.log(2); lf = mp.log(phi); z2 = mp.pi**2/6; z3 = mp.zeta(3)
Li3p = mp.polylog(3, 1/phi); Li3m = mp.polylog(3, -1/phi)
vals = {
 'l2^3': l2**3, 'l2^2*lf': l2**2*lf, 'l2^2': l2**2, 'l2*lf^2': l2*lf**2,
 'l2*lf': l2*lf, 'l2*z2': l2*z2, 'l2': l2, 'lf^3': lf**3, 'lf^2': lf**2,
 'lf*z2': lf*z2, 'lf': lf, 'z2': z2, 'z3': z3, 'L3(-1/phi)': Li3m, 'L3(1/phi)': Li3p}
S = mp.mpf(0)
for k,(p,q) in ft.items():
    S += (mp.mpf(str(p)) + mp.mpf(str(q))*sqrt5)*vals[k]
P3_cf = 60*L*S
print("60*L*S     =", mp.nstr(P3_cf, 55))
print("P'''(1)    = -4.776200180882238416632022453444885836964")
P3 = mp.mpf("-4.776200180882238416632022453444885836964")
print("diff:", mp.nstr(abs(P3-P3_cf), 4))
