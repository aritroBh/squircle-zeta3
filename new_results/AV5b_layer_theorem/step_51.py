import os
work = "<project>/research/08_theory/AV5b_layer_theorem"
os.chdir(work)
step = r'''
import mpmath as mp
import json
mp.mp.dps = 70
prof = [(mp.mpf(s), mp.mpf(a), mp.mpf(b), mp.mpf(c)) for s,a,b,c in json.load(open("ckpt_prof_m2.json"))]
ss = [p[0] for p in prof]; ds = ss[1]-ss[0]
def trapz(vals):
    return mp.fsum((vals[i]+vals[i+1])/2*ds for i in range(len(vals)-1))
g1 = [p[1] for p in prof]; g2 = [p[2] for p in prof]; g3 = [p[3] for p in prof]
A = trapz([g3[i] + g1[i]**3/6 for i in range(len(prof))])
B = trapz([abs(prof[i][0])*g1[i]**2/2 for i in range(len(prof))])
D = trapz([-(prof[i][0]**2/2)*g1[i] for i in range(len(prof))])
print("A = int(g3 + g1^3/6) =", mp.nstr(A, 25))
print("B = int |s| g1^2/2 =", mp.nstr(B, 25))
print("D = -int s^2 g1/2 =", mp.nstr(D, 25))
print("H*(A+B+D) =", mp.nstr((A+B+D)/mp.sqrt(2), 25))
# analytic pieces
print("int_R g3 (analytic -7z4/8 *?):")
g3v = 2*mp.quad(lambda s: (mp.mpf(2)/3)*s**3*(mp.tanh(s)-1), [0, mp.inf])
print("  int_0inf g3 =", mp.nstr(g3v/2, 25), "  -7 zeta4/16 =", mp.nstr(-7*mp.zeta(4)/16, 25))
print("  int_R g3 =", mp.nstr(g3v, 25))
print("int_0inf g1^3 =", mp.nstr(mp.quad(lambda s: mp.log(1+mp.exp(-2*s))**3, [0, mp.inf]), 25))
print("int_0inf s g1^2 =", mp.nstr(mp.quad(lambda s: s*mp.log(1+mp.exp(-2*s))**2, [0, mp.inf]), 25))
print("int_0inf s^2 g1 =", mp.nstr(mp.quad(lambda s: s**2*mp.log(1+mp.exp(-2*s)), [0, mp.inf]), 25), " 7z4/32 =", mp.nstr(7*mp.zeta(4)/32, 25))
'''
with open("av5b_debug.py","w") as f: f.write(step)
exec(step)
