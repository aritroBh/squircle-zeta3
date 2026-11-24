import pickle
for f in ['../OCT6_octagon_closed/final_tables.pkl','../OCT6_octagon_closed/S_final.pkl']:
    try:
        d = pickle.load(open(f,'rb'))
        print('='*20, f, type(d))
        if isinstance(d, dict):
            for k in list(d)[:10]: print('  key:', k, '->', str(d[k])[:200])
    except Exception as e:
        print(f, 'ERR', e)
