import os, mpmath as mp, pickle, time
os.chdir('research/08_theory/DEC1_decagon')
mp.mp.dps = 830
s5=mp.sqrt(5); phi=(1+s5)/2; lf=mp.log(phi); l2=mp.log(2); z2=mp.pi**2/6; z3=mp.zeta(3)
x1=-1/phi; x4=1/phi**2; xf=1/phi
N=4200
def Li3ser(x):
    tot=mp.mpf(0); xp=mp.mpf(1)
    for n in range(1,N):
        xp*=x; tot+=xp/n**3
    return tot
def Aser(x):
    tot=mp.mpf(0); xp=mp.mpf(1); H=mp.mpf(0)
    for n in range(1,N):
        xp*=x; H+=mp.mpf(1)/n; tot+=H*xp/n**2
    return tot
def Li21ser(x,y):
    tot=mp.mpf(0); ps=mp.mpf(0); xp=mp.mpf(1); yp=mp.mpf(1)
    for m in range(1,N):
        xp*=x
        if m>=2: tot += xp*ps/m**2
        yp*=y; ps += yp/m
    return tot
t0=time.time()
vals={}
vals['L3_x4']=Li3ser(x4); print("L3(1/phi^2) done", time.time()-t0)
vals['L3_x1']=Li3ser(x1)
vals['L3_xf']=Li3ser(xf)
vals['A_x1']=Aser(x1); print("A done", time.time()-t0)
vals['A_x4']=Aser(x4)
vals['Lam']=Li21ser(x1,x1); print("Lam done", time.time()-t0)
vals['swap']=Li21ser(x4,-phi); print("swap done", time.time()-t0)
# tail estimates: terms ~ phi^-m poly: next term magnitude
print("tail check Lam term:", mp.nstr(abs((x1)**N*( -mp.log(1-x1))/N**2),6))
for k,v in vals.items(): print(k, mp.nstr(v,40))
# cross-checks vs 100-dps quad values
print("Lam check:", mp.nstr(abs(vals['Lam']-mp.mpf('-0.0512616305864763273812632603657')),8))
print("swap check:", mp.nstr(abs(vals['swap']-mp.mpf('-0.0634178216821646121508597')),8))
vals.update({'lf':lf,'l2':l2,'z2':z2,'z3':z3,'phi':phi})
pickle.dump(vals, open('atoms830.pkl','wb'))
with open('progress.log','a') as f:
    f.write(f"step: 830-dps atoms: L3 at {{1/phi^2,-1/phi,1/phi}}, A(x1),A(x4), Lam=Li21(x1,x1), swap=Li21(x4,-phi). saved atoms830.pkl\n")