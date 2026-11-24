import pickle
import mpmath as mp
mp.mp.dps = 800
W8 = pickle.load(open('wparts_dec_d800.pkl','rb'))
r5 = mp.sqrt(5); l2 = mp.log(2); lf = mp.log((1+r5)/2); z2 = mp.zeta(2)
Wt1_exact = 1200*l2 + (3000-600*r5)*lf
Wt2_exact = (mp.mpf(-1425)/16 - mp.mpf(5415)/16*r5)*z2 - 2400*l2**2 + (-12000+2400*r5)*l2*lf + (mp.mpf(-86925)/4 + mp.mpf(17325)/2*r5)*lf**2
print("Wt0 numeric:", mp.nstr(mp.mpf(W8['Wt0']), 3), "(exact 0)")
print("Wt1: |numeric - exact| =", mp.nstr(abs(mp.mpf(W8['Wt1']) - Wt1_exact), 3))
print("Wt2: |numeric - exact| =", mp.nstr(abs(mp.mpf(W8['Wt2']) - Wt2_exact), 3))
open('progress.log','a').write("step_34: EXACT Wt0=0, Wt1=1200*log2+(3000-600*sqrt5)*logphi, Wt2=(-1425/16-5415/16*sqrt5)*z2-2400*log2^2+(-12000+2400*sqrt5)*log2*logphi+(-86925/4+17325/2*sqrt5)*logphi^2; vs 800-digit: Wt1 %s, Wt2 %s\n" % (
    mp.nstr(abs(mp.mpf(W8['Wt1']) - Wt1_exact), 2), mp.nstr(abs(mp.mpf(W8['Wt2']) - Wt2_exact), 2)))
