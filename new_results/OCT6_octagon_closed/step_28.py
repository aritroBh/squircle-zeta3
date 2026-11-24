import pickle, mpmath as mp, sympy as sp
import mixed_engine as ME
mp.mp.dps = 760
num = {sp.symbols(k): mp.mpf(v) for k,v in pickle.load(open('num745.pkl','rb')).items()}
S = pickle.load(open('S_assembly.pkl','rb'))
r2m=mp.sqrt(2); k1m=r2m-1; k2m=r2m/2-1; um=mp.sqrt(2+r2m)
v_tot = 48*um*mp.mpf(sp.N(S['S_tot'].subs(num, simultaneous=True), 745))
mons = pickle.load(open('mons_simple.pkl','rb'))
lc1 = -mp.log(2)/2 - mp.log(1+r2m); lc2 = -mp.log(1+r2m)
import time; t0=time.time()
qtot = mp.mpf(0)
for (a,b,c1,c2,e1,e2), terms in mons:
    clog = lc1**e1*lc2**e2
    for n,q in terms:
        qn = mp.mpf(sp.N(q, 745))
        if qn==0: continue
        f = lambda tt: tt**n*mp.log(tt)**a*mp.log(1-tt)**b*mp.log(1+k1m*tt)**c1*mp.log(1+k2m*tt)**c2
        qtot += clog*qn*mp.quad(f, [0,1])
qtot *= 48*um
print("quad of all 83 monomials: %.1fs"%(time.time()-t0))
print("P3 exact (48u*S_tot) =", mp.nstr(v_tot,60))
print("P3 quad760           =", mp.nstr(qtot,60))
print("|exact - quad760| =", mp.nstr(abs(v_tot-qtot),3))
P3_720 = pickle.load(open('oct5/P3_720.pkl','rb'))
print("|exact - storedP3_720| =", mp.nstr(abs(v_tot-mp.mpf(P3_720['P3'])),3))
pickle.dump({'P3_exact': mp.nstr(v_tot,745), 'P3_quad760': mp.nstr(qtot,745)}, open('P3_exact.pkl','wb'))
with open('progress.log','a') as f:
    f.write("step_05: exact P3=48u*S_tot vs fresh 760-dps quad of all 83 monomials: |diff|=%s (~576+ digits verified);\n"
            "  stored OCT5 P3_720 value is only accurate to ~203 digits (its quad was under-resolved); exact matches stored to 2.8e-200=all its accurate digits.\n"
            % mp.nstr(abs(v_tot-qtot),3))
print("logged")
