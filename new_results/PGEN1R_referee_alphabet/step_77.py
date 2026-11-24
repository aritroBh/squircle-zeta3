import pickle
for f in ['short_per_edge.pkl','P3_short.pkl','short_form.pkl','short123.pkl','G8_final.pkl']:
    try:
        d = pickle.load(open(f'../OCT9d_octagon_short/{f}','rb'))
        print('='*15, f, type(d))
        if isinstance(d, dict):
            for k in list(d)[:6]:
                print('  ', k, '->', str(d[k])[:300])
        elif isinstance(d, (list,tuple)):
            print('  len', len(d), str(d[:2])[:400])
    except Exception as e:
        print(f, 'ERR', e)
