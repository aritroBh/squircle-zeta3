import os
work = "<project>/research/08_theory/AV5b_layer_theorem"
os.chdir(work)
step = r'''
import mpmath as mp
import json
mp.mp.dps = 70

def load_prof(mm):
    return [(mp.mpf(s), mp.mpf(a), mp.mpf(b), mp.mpf(c)) for s,a,b,c in json.load(open(f"ckpt_prof_m{mm}.json"))]
def trapz(ss, ds, vals):
    return mp.fsum((vals[i]+vals[i+1])/2*ds for i in range(len(ss)-1))

for mm in [3, 4]:
    prof = load_prof(mm)
    ss = [p[0] for p in prof]; ds = ss[1]-ss[0]
    H = mp.sin(mp.pi/(2*mm)); alpha = mp.tan(mp.pi/(2*mm))
    a1 = prof[-1][1]  # D1 at s=20
    a2 = prof[-1][2]
    g1 = [p[1]-a1 for p in prof]; g2 = [p[2]-a2 for p in prof]
    I_g2   = trapz(ss, ds, g2)
    I_g1sq = trapz(ss, ds, [x**2/2 for x in g1])
    I_sg1  = trapz(ss, ds, [(alpha*abs(ss[i])+a1)*g1[i] for i in range(len(ss))])
    c3 = H*(I_g2 + I_g1sq + I_sg1)
    z3 = mp.zeta(3)
    print(f"m={mm}: a1={mp.nstr(a1,18)} a2={mp.nstr(a2,14)}")
    print(f"   c3_lay = {mp.nstr(c3, 25)}   c3/zeta3 = {mp.nstr(c3/z3, 20)}")
    C_m = 48*mm*mp.sin(mp.pi/(2*mm))**4*mp.cos(mp.pi/(2*mm))*(3*mp.sin(mp.pi/(2*mm))**2-1)
    print(f"   C_m/12m target = {mp.nstr(C_m/12/mm, 20)}  c3 - C_m z3/12m = {mp.nstr(c3 - C_m*z3/(12*mm), 8)}")
'''
with open("av5b_c3m4.py","w") as f: f.write(step)
exec(step)
