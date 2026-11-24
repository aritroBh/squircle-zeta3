import pickle, mpmath as mp, sympy as sp
mp.mp.dps = 745
r2=mp.sqrt(2); k1=r2-1; k2=r2/2-1; x1=-k1; x2=-k2
l2=mp.log(2); lu=mp.log(1+r2); z2=mp.pi**2/6; z3=mp.zeta(3)
reps = {'m1':-1, 'mr2h':-r2/2, 'x1':x1, 'x2':x2, 'k1':k1, 'h':mp.mpf('0.5'), 'd':2-r2, 'r2h':r2/2}
Li2v = {nm: mp.polylog(2,a) for nm,a in reps.items()}
Li3v = {nm: mp.polylog(3,a) for nm,a in reps.items()}
# closed-form checks
print("Li2(-1)+z2/2:", mp.nstr(abs(Li2v['m1']+z2/2),3))
print("Li2(1/2)-(z2/2-l2^2/2):", mp.nstr(abs(Li2v['h']-(z2/2-l2**2/2)),3))
print("Li3(-1)+3z3/4:", mp.nstr(abs(Li3v['m1']+3*z3/4),3))
print("Li3(1/2)-(7z3/8-z2*l2/2+l2^3/6):", mp.nstr(abs(Li3v['h']-(7*z3/8-z2*l2/2+l2**3/6)),3))
vals = {'1':mp.mpf(1),'l2':l2,'lu':lu,'z2':z2,'z3':z3}
for nm in reps: vals['Li2_'+nm]=Li2v[nm]; vals['Li3_'+nm]=Li3v[nm]
# weight-3 products
vals['z2l2']=z2*l2; vals['z2lu']=z2*lu
vals['l2^3']=l2**3; vals['l2^2lu']=l2**2*lu; vals['l2lu^2']=l2*lu**2; vals['lu^3']=lu**3
# weight-2 products (for Li2 relations)
vals['l2^2']=l2**2; vals['l2*lu']=l2*lu; vals['lu^2']=lu**2
# Li2(rep)*log products (weight 3)
for nm in reps:
    if nm in ('m1','h'): continue
    vals['Li2_'+nm+'*l2']=Li2v[nm]*l2; vals['Li2_'+nm+'*lu']=Li2v[nm]*lu
# A atoms + Lam + swap
num = {sp.symbols(k): mp.mpf(v) for k,v in pickle.load(open('num745.pkl','rb')).items()}
import mixed_engine as ME
vals['AA_1']=num[ME.AA_1]; vals['AA_2']=num[ME.AA_2]
vals['Lam']=num[ME.Lam]; vals['swap']=num[ME.swap]
pickle.dump(vals, open('basis745.pkl','wb'))
print("basis745 saved;", len(vals), "named values at 745 digits")
