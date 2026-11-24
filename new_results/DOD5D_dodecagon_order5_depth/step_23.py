import pickle, time
from mpmath import mp, mpf, nstr, sqrt
# Sanity vs known values
mp.dps = 80
P1_engine = mpf("1.66440849890435430350045161841510269829689839051162011648807")
L6 = sqrt(26-15*sqrt(3))
print("P(1) vs 12*L6 residual:", nstr(abs(P1_engine - 12*L6), 4))
P3_engine = mpf("-5.20324042951766621029582970427243461666280533434661632110814")
P3_630 = pickle.load(open("<project>/research/08_theory/DOD1_dodecagon/P3_630.pkl","rb"))
print(type(P3_630), str(P3_630)[:80])
