import pickle
for f in ["basis_indep.pkl","basis_names_exact.pkl","basis_names.pkl"]:
    try:
        d = pickle.load(open(f,"rb"))
        print("="*20, f, type(d), len(d) if hasattr(d,'__len__') else '')
        print(repr(d)[:2500])
    except Exception as e:
        print(f,"ERR",e)
