import os
work = "<project>/research/08_theory/AV5b_layer_theorem"
os.chdir(work)
import mpmath as mp, json
mp.mp.dps = 60
from scipy.interpolate import CubicSpline

# m=4 pieces from extracted profiles
prof = [(mp.mpf(s), mp.mpf(a), mp.mpf(b), mp.mpf(c)) for s,a,b,c in json.load(open("ckpt_prof_m4.json"))]
m = 4
H = mp.sin(mp.pi/(2*m)); alpha = mp.tan(mp.pi/(2*m))
a1 = prof[-1][1]; a2 = prof[-1][2]; a3 = prof[-1][3]
sf = [float(p[0]) for p in prof]
g1v = [float(p[1]-a1) for p in prof]; g2v = [float(p[2]-a2) for p in prof]
sp1 = CubicSpline(sf, g1v, bc_type='natural'); sp2 = CubicSpline(sf, g2v, bc_type='natural')
g1 = lambda s: mp.mpf(str(sp1(float(s)))); g2 = lambda s: mp.mpf(str(sp2(float(s))))
Ig1 = 2*mp.quad(g1, [0, 20])
Ig1sq = 2*mp.quad(lambda s: g1(s)**2/2, [0, 20])
Isg1 = 2*mp.quad(lambda s: alpha*s*g1(s), [0, 20])
Ig2 = 2*mp.quad(g2, [0, 20])
P1 = 8*(mp.sqrt(2)-1)*mp.sqrt(2-mp.sqrt(2))
A3 = mp.mpf("-4.39449630527605215997911011380609991382237745")/6
c1 = a1*P1
c2 = (a2 + a1**2/2)*P1 + 2*m*H*Ig1
c3 = (a3 + a1*a2 + a1**3/6)*P1 + 2*m*H*(a1*Ig1 + Ig2 + Ig1sq + Isg1)
A3_lemma = c1 - 2*c2 + c3
print("m=4: A3 exact =", mp.nstr(A3, 30))
print("m=4: A3 lemma =", mp.nstr(A3_lemma, 30))
print("residual:", mp.nstr(A3 - A3_lemma, 6))
# also back out int g2 as an identity check
Ig2_back = (A3 - P1*(a1 - 2*(a2+a1**2/2) + a3 + a1*a2 + a1**3/6))/(2*m*H) - (a1*Ig1 + Ig1sq + Isg1 - 2*Ig1)
print("int g2 (extracted) =", mp.nstr(Ig2, 20))
print("int g2 (back-out)  =", mp.nstr(Ig2_back, 20))
