import mpmath as mp, pickle, time
mp.mp.dps = 620
t0=time.time()
D1 = '../DOD1_dodecagon/'
atoms = pickle.load(open(D1+'atoms620.pkl','rb'))
X = atoms['X']; T3a = atoms['T3a']
g = pickle.load(open('gdswaps620.pkl','rb'))
SwG, SwD = g['SwG'], g['SwD']
B = pickle.load(open('basis102_620.pkl','rb'))
Q = pickle.load(open('basis_quotient.pkl','rb'))
keep54 = [i for i in Q['keep'] if i < 54]
Bn = [B['names'][i] for i in keep54]; Bv = [B['vals'][i] for i in keep54]
R = T3a[(2,5,6)] - (SwG[(2,5,6)]-SwD[(2,5,6)])/X[2]
print("pslq R(2,5,6) vs B54_ind (32 elts), t=%ds"%(time.time()-t0))
rel = mp.pslq([R]+Bv, tol=mp.mpf('1e-570'), maxcoeff=10**10, maxsteps=30000)
print("pslq done t=%ds"%(time.time()-t0))
if rel is None:
    print("PSLQ: NO relation maxcoeff 1e10")
else:
    resid = rel[0]*R + sum(rel[i+1]*Bv[i] for i in range(len(Bv)))
    print("resid:", mp.nstr(abs(resid),4), "tcoef:", rel[0], "maxc:", max(abs(x) for x in rel))
    print([(Bn[i],rel[i+1]) for i in range(len(Bv)) if rel[i+1]!=0])
