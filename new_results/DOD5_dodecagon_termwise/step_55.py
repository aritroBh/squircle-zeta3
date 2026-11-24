import sys, os, pickle, time
base = "<project>/research/08_theory/"
sys.path.insert(0, base+"DOD5_dodecagon_termwise/")
from dod4_lib import *
t0=time.time()
mp.mp.dps = 620
atoms = pickle.load(open(D1+"atoms620.pkl","rb"))
Xs = {1: P(3,-2,3), 2: P(1,-1), 5: P(3,-1,3), 6: P(2,-1,2)}
IDX=[1,2,5,6]
Li2R = lambda z: mp.re(mp.polylog(2, mp.mpc(z)))
Li3R = lambda z: mp.re(mp.polylog(3, mp.mpc(z)))

# 1) A-atom Nielsen: A(x) = Li3(x) + 1/2[ln|x|ln^2(1-x) + 2 ln(1-x) Li2(1-x) - 2 Li3(1-x) + 2 z3]
z3 = mp.zeta(3)
worst = mp.mpf(0)
for i in IDX:
    x = atoms['X'][i]
    cf = Li3R(x) + (mp.log(abs(x))*mp.log(abs(1-x))**2 + 2*mp.log(abs(1-x))*Li2R(1-x) - 2*Li3R(1-x) + 2*z3)/2
    d = abs(cf - atoms['Ax'][i]); worst = max(worst,d)
print("A-atom Nielsen vs atoms['Ax'], worst:", mp.nstr(worst,3))

# 2) T3p formula
def Fnie(z):
    return mp.log(abs(z))*mp.log(abs(1-z))**2 + 2*mp.log(abs(1-z))*Li2R(1-z) - 2*Li3R(1-z)
def T3p_formula(p, q):
    gam = q/(q-p); dlta = gam*(1-p)
    A = mp.log(abs((p-q)/p))
    return ( A*A*(mp.log(1)-mp.log(1-p)) - 2*A*(Li2R(gam)-Li2R(dlta)) + (Fnie(gam)-Fnie(dlta)) )/p
worst = mp.mpf(0)
for i in IDX:
    for j in IDX:
        if i==j: continue
        p, q = atoms['X'][i], atoms['X'][j]
        d = abs(T3p_formula(p,q) - atoms['T3p'][(i,j)]); worst = max(worst,d)
print("T3p_formula vs atoms['T3p'] (12 pairs), worst:", mp.nstr(worst,3))

# 3) exact Li3 args of T3p: (2/p)[Li3(1-dl) - Li3(1-gam)]
gam={}; dl={}
for i in IDX:
    for j in IDX:
        if i==j: continue
        gam[(i,j)] = pdiv(Xs[j], psub(Xs[j], Xs[i]))
        dl[(i,j)]  = pmul(gam[(i,j)], psub(ONE, Xs[i]))
pts_T3p = set()
for i in IDX:
    for j in IDX:
        if i==j: continue
        pts_T3p.add(psub(ONE, gam[(i,j)])); pts_T3p.add(psub(ONE, dl[(i,j)]))
print("distinct T3p Li3 args (1-gam, 1-dl):", len(pts_T3p))
for z in sorted(pts_T3p, key=pfloat):
    f = factor_sunit(z)
    print("   %-14s ~%+.8f  sunit=%s"%(str(z), pfloat(z), f!='NOT_SUNIT' and f is not None and f[0]!='NOT_SUNIT'))
# also x_i and 1-x_i
pts_AX = set(Xs.values()) | {psub(ONE,Xs[i]) for i in IDX}
print("x_i, 1-x_i all S-units:", all(factor_sunit(z)[0]!='NOT_SUNIT' for z in pts_AX))
pickle.dump({'pts_T3p':sorted(pts_T3p,key=pfloat), 'pts_AX':sorted(pts_AX,key=pfloat)}, open('dod5_sd1_li3pts.pkl','wb'))
print("elapsed %.0fs"%(time.time()-t0))
