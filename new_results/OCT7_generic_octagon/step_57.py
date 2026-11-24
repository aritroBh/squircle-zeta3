import pickle
d = pickle.load(open('lam_lll_450.pkl','rb'))
for k, v in d['results'].items():
    print(k, v[0])
    if v[0]=='REDUCES':
        import ast
        rel = ast.literal_eval(v[1])
        nz = [(i,c) for i,c in enumerate(rel) if c != 0]
        print("   nonzero entries:", nz)
        names = pickle.load(open('basis_names.pkl','rb'))['basis_names']
        for i,c in nz:
            print("   ", c, "*", names[i] if i < len(names) else "TARGET")
