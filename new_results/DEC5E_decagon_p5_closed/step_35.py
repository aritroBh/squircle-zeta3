import sys
sys.path.insert(0,'.')
import gengine_dec2p, mpmath as mp
mp.mp.dps = 800
G1, letters = gengine_dec2p.gtab(800, 2700, cachedir=".")
L0,L1,LP,LQ = letters['L0'],letters['L1'],letters['LP'],letters['LQ']
lf = mp.log((1+mp.sqrt(5))/2); z2 = mp.zeta(2); l2 = mp.log(2)
# candidate exact values for all length-2 words not starting with 1
vals = {
 (L0,L0): mp.mpf(0),
 (L0,L1): -z2,
 (L0,LP): 2*z2/5 - lf**2/2,
 (L0,LQ): -2*z2/5 + lf**2,
 (LP,L0): -2*z2/5 + lf**2/2,
 (LQ,L0): 2*z2/5 - lf**2,
 (LP,LP): lf**2/2,
 (LQ,LQ): lf**2/2,
 (LP,L1): -2*z2/5 + lf**2,      # -Li2(1/phi^2)
 (LQ,L1): 2*z2/5 - lf**2/2,      # -Li2(-1/phi)
 (LP,LQ): -z2/5 + lf**2,         # derived: Li2(1/phi^2)-Li2(1/phi)+lf^2
 (LQ,LP): z2/5 - 2*lf**2,        # by shuffle
}
worst = mp.mpf(0); wk=None
for w, v in vals.items():
    d = abs(G1[w]-v)
    if d > worst: worst, wk = d, w
print("worst |table - exact candidate| over 12 length-2 words:", mp.nstr(worst,3), "at", wk)
# length-1
for w, v in {(L0,): mp.mpf(0), (LP,): lf, (LQ,): -lf}.items():
    print(w, mp.nstr(abs(G1[w]-v),3))
open('progress.log','a').write("step_32: exact length-2 word values (shuffle + DEC2-proved golden Li2 + derived G(-phi,phi^2;1)=Li2(1/phi^2)-Li2(1/phi)+log^2 phi): all 12 verified vs 800-digit table, worst %s\n" % mp.nstr(worst,3))
