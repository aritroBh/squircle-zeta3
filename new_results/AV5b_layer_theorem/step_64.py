import os
work = "<project>/research/08_theory/AV5b_layer_theorem"
os.chdir(work)
exec(open("av5b_locmodel.py").read().split('def extract')[0])
mp.mp.dps = 50
Zo = Loc(3)
tau = mp.mpf('0.006'); p = 1/(1-tau)
for s in [mp.mpf('2'), mp.mpf('6'), mp.mpf('12')]:
    u = tau*s
    # fine grid max of val(t)
    vals = []
    Nt = 400
    for k in range(1, Nt):
        t = mp.mpf(k)/Nt
        z = Zo.zs(t, p)
        xt0 = Zo.v0[0]+t*Zo.dv[0]; xt1 = Zo.v0[1]+t*Zo.dv[1]
        mx = (Zo.v0[0]+Zo.v1[0])/2; my = (Zo.v0[1]+Zo.v1[1])/2
        Td = (Zo.dv[0]/Zo.L, Zo.dv[1]/Zo.L)
        ell = Td[0]*(xt0-mx)+Td[1]*(xt1-my)
        vals.append(((Zo.H+z)*mp.cos(u) + ell*mp.sin(u), t))
    best = max(vals)
    h = best[0]
    m_ = Zo.tent(u)
    print(f"s={mp.nstr(s,4)}: u={mp.nstr(u,7)} h_loc={mp.nstr(h,14)} tent={mp.nstr(m_,14)} D={mp.nstr(mp.log(h/m_)/tau,14)} t*={mp.nstr(best[1],6)}")
    # also the h(u) from class
    print("   class h:", mp.nstr(Zo.h(u, p), 14))
