import pickle
from mpmath import mp, mpf, nstr, sqrt
mp.dps = 190
res = pickle.load(open("p5_170.pkl","rb"))["res"]
P = [mpf(x) for x in res]
L6 = sqrt(26-15*sqrt(3))
print("P(1) vs 12*L6 resid:", nstr(abs(P[0]-12*L6), 4))
P3_630 = pickle.load(open("<project>/research/08_theory/DOD1_dodecagon/P3_630.pkl","rb"))
print("P'''(1) vs P3_630 resid:", nstr(abs(P[3]-P3_630['P3']), 4))
print("P5 120 vs 170 agree:", nstr(abs(P[5]-mpf('-62.5372477549863953530912928220580357268476927741037326268814581807469573142287452821422072155921812848246065294308796281')),4))
# u normalization: u = sqrt(2-sqrt3) = (sqrt6-sqrt2)/2 ; L6 = u*(2-sqrt3)
u = sqrt(2-sqrt(3))
print("L6/u =", nstr(L6/u, 30), " 2-sqrt3 =", nstr(2-sqrt(3),30))
print("u - (sqrt6-sqrt2)/2 =", nstr(abs(u-(sqrt(6)-sqrt(2))/2), 4))
print("P5/u =", nstr(P[5]/u, 120))
