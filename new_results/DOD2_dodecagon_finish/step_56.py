import pickle
Q = pickle.load(open('basis_quotient.pkl','rb'))
B = pickle.load(open('basis102_620.pkl','rb'))
# do elim expressions reference e-point columns (index >= 54)?
for pv, expr in Q['elim'].items():
    hi = [j for j in expr if j >= 54]
    if hi:
        print("pivot", pv, B['names'][pv], "references", [(j, B['names'][j]) for j in hi])
print("---")
keep54 = [i for i in Q['keep'] if i < 54]
print("B54 elements kept:", len(keep54), [B['names'][i] for i in keep54])
