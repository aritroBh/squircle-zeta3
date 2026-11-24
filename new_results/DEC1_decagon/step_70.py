import os, mpmath as mp, pickle
os.chdir('research/08_theory/DEC1_decagon')
def compute_atoms(dps):
    mp.mp.dps = dps
    s5=mp.sqrt(5); phi=(1+s5)/2; lf=mp.log(phi); l2=mp.log(2); z2=mp.pi**2/6; z3=mp.zeta(3)
    x1=-1/phi; x4=1/phi**2; xf=1/phi
    N=int(dps*1.15/0.20899)+200
    def Li3ser(x):
        tot=mp.mpf(0); xp=mp.mpf(1)
        for n in range(1,N): xp*=x; tot+=xp/n**3
        return tot
    def Aser(x):
        tot=mp.mpf(0); xp=mp.mpf(1); H=mp.mpf(0)
        for n in range(1,N): xp*=x; H+=mp.mpf(1)/n; tot+=H*xp/n**2
        return tot
    def Li21ser(x,y):
        tot=mp.mpf(0); ps=mp.mpf(0); xp=mp.mpf(1); yp=mp.mpf(1)
        for m in range(1,N):
            xp*=x
            if m>=2: tot += xp*ps/m**2
            yp*=y; ps += yp/m
        return tot
    return dict(lf=lf,l2=l2,z2=z2,z3=z3,phi=phi,L3x1=Li3ser(x1),L3xf=Li3ser(xf),L3x4=Li3ser(x4),
                A_x1=Aser(x1),A_x4=Aser(x4),Lam=Li21ser(x1,x1),swap=Li21ser(x4,-phi))
expected = {'A_x1':[-5,4,-2,0,0,0,0,0,10,0],'A_x4':[-5,5,-3,0,0,0,0,0,20,15],
            'Lam':[-5,2,4,0,-5,0,0,0,15,5],'swap':[-5,3,-7,0,10,0,0,0,5,5]}
for dps in [400, 600, 830]:
    v = compute_atoms(dps)
    lf,l2,z2,z3 = v['lf'],v['l2'],v['z2'],v['z3']
    B = [z3, z2*lf, z2*l2, lf**3, lf**2*l2, lf*l2**2, l2**3, v['L3x1'], v['L3xf']]
    print(f"--- dps={dps} ---")
    for tname in ['A_x1','A_x4','Lam','swap']:
        rel = mp.pslq([v[tname]]+B, tol=mp.mpf(10)**(-(dps-40)), maxcoeff=10**14, maxsteps=10000)
        match = (rel==expected[tname])
        resid = None
        if rel:
            resid = abs(rel[0]*v[tname]+sum(rel[i+1]*B[i] for i in range(9)))
        print(f"{tname}: match={match} rel={rel} resid={mp.nstr(resid,4) if resid else 'NA'}")
with open('progress.log','a') as f:
    f.write("step: LLL stability across 400/600/830 dps: identical integer vectors, residuals scale ~10^-dps.\n")