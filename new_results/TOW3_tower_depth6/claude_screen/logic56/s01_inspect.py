import pickle, sympy as sp, sys
sys.dont_write_bytecode=True
for f in ['e_fits.pkl','e61_closed.pkl','arho56_simplified.pkl','Sderiv_fits.pkl']:
    o = pickle.load(open(f,'rb'))
    print(f, type(o), (list(o.keys()) if isinstance(o,dict) else ''))
    if isinstance(o,dict):
        for k,v in o.items():
            print('  ',k, str(v)[:300])
