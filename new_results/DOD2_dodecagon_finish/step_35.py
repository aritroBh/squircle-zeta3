import mpmath as mp, pickle, time
mp.mp.dps = 620
t0=time.time()
D1 = '../DOD1_dodecagon/'
atoms = pickle.load(open(D1+'atoms620.pkl','rb'))
X=atoms['X']; T3a=atoms['T3a']
g = pickle.load(open('gdswaps620.pkl','rb'))
SwG, SwD = g['SwG'], g['SwD']
B = pickle.load(open('basis102_620.pkl','rb'))
Q = pickle.load(open('basis_quotient.pkl','rb'))
keep = Q['keep']
Bn = [B['names'][i] for i in keep]; Bv = [B['vals'][i] for i in keep]
R = T3a[(1,2,5)] - (SwG[(1,2,5)]-SwD[(1,2,5)])/X[1]
print("pslq R1_25 vs B_ind (80 elts)...")
rel = mp.pslq([R]+Bv, tol=mp.mpf('1e-580'), maxcoeff=10**9, maxsteps=40000)
print("pslq done t=%ds"%(time.time()-t0))
if rel is None:
    print("PSLQ: NO relation maxcoeff 1e9")
else:
    resid = rel[0]*R + sum(rel[i+1]*Bv[i] for i in range(len(Bv)))
    print("resid:", mp.nstr(abs(resid),4), "target coef:", rel[0], "maxc:", max(abs(x) for x in rel))
    print([(Bn[i],rel[i+1]) for i in range(len(Bv)) if rel[i+1]!=0])
