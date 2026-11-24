import pickle, json
for f in ['fit_wt3.pkl','fit_wt4.pkl','fit_R_lchi.pkl']:
    d = pickle.load(open(f,'rb'))
    print('='*20, f, 'height', d['height'])
    rel = d.get('rel')
    print(json.dumps({k: (str(v[0]), str(v[1])) for k, v in sorted(rel.items())}, indent=0))
