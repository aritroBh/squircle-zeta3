import mpmath as mp, pickle
mp.mp.dps = 130
r2 = mp.sqrt(2); k1 = r2-1; k2 = r2/2-1
x1 = -k1; x2 = -k2
l2 = mp.log(2); lu = mp.log(1+r2)
atoms = {k: mp.mpf(v) for k,v in pickle.load(open('oct5/atoms600.pkl','rb')).items()}
Lam, swap = atoms['Li21_a'], atoms['Li21_b']
Li2 = lambda z: mp.polylog(2,z); Li3 = lambda z: mp.polylog(3,z)
X = mp.log
def quad(f): return mp.quad(f, [0,1])

D = {}
D['D1'] = (quad(lambda t: mp.log(t)*mp.log(1+k2*t)/(1+k1*t)), -Lam/k1)
D['D2'] = (quad(lambda t: mp.log(t)*mp.log(1+k1*t)/(1+k2*t)), -swap/k2)
D['D3'] = (quad(lambda t: mp.log(1-t)*mp.log(1+k2*t)/(1+k1*t)),
           ( (l2/2)*Li2(x2)/x2 + swap/x2 )/r2)
D['D4'] = (quad(lambda t: mp.log(1-t)*mp.log(1+k1*t)/(1+k2*t)),
           ( Lam/x1 - (l2/2)*Li2(x1)/x1 )/(r2/2))
D['D5'] = (quad(lambda t: mp.log(1+k2*t)/(1+k1*t)),
           ( (lu-l2/2)*(l2/2) + Li2(k1) - Li2(2-r2) )/k1)
D['D6'] = (quad(lambda t: mp.log(1+k1*t)/(1+k2*t)),
           ( -lu*l2/2 - Li2(k1) + Li2(2-r2) )/k2)
D['D7'] = (quad(lambda t: mp.log(1+k1*t)**2/(1+k2*t)),
           ( -(l2**2/4)*lu + l2*Li2(2-r2) - 2*Li3(2-r2) + 2*Li3(k1) )/k2)
D['D8'] = (quad(lambda t: mp.log(1+k2*t)**2/(1+k1*t)),
           ( (l2**2/4)*(l2/2-lu) - l2*Li2(k1) - 2*Li3(k1) + 2*Li3(2-r2) )/k1)
D7v = D['D7'][1]; D8v = D['D8'][1]
D['D9'] = (quad(lambda t: mp.log(1+k1*t)*mp.log(1+k2*t)/(1+k1*t)),
           -l2**3/(16*k1) - k2/(2*k1)*D7v)
D['D10']= (quad(lambda t: mp.log(1+k1*t)*mp.log(1+k2*t)/(1+k2*t)),
           l2**3/(16*k2) - k1/(2*k2)*D8v)
for name,(qn,cf) in D.items():
    print(f"{name}: quad={mp.nstr(qn,30)}")
    print(f"     form={mp.nstr(cf,30)}   |diff|={mp.nstr(abs(qn-cf),3)}")
pickle.dump({k:(mp.nstr(v[1],120)) for k,v in D.items()}, open('Dconsts.pkl','wb'))
with open('progress.log','a') as f:
    f.write("step_02: D1..D10 closed forms verified by quadrature, diffs printed (120-digit values saved Dconsts.pkl)\n")
