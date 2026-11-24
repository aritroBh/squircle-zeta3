import pickle
base = '<project>/research/08_theory'
work = f'{base}/TOW3_tower_depth6'
for f in ['L_simplified.pkl','rho_closed_forms.pkl','e_fits.pkl','Sderiv_fits.pkl']:
    obj = pickle.load(open(f'{work}/{f}','rb'))
    print('='*70); print(f, type(obj))
    if isinstance(obj, dict):
        for k in sorted(obj.keys(), key=str):
            v = str(obj[k])
            print(' ', k, ':', v[:200])
