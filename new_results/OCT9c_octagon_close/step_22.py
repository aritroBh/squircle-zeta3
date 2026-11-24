import pickle
d='../OCT9b_octagon_debug/'
lam = pickle.load(open(d+'lam_cache.pkl','rb'))
print(type(lam), list(lam.keys())[:20] if isinstance(lam,dict) else len(lam))
for k,v in (lam.items() if isinstance(lam,dict) else []):
    print(k, '->', str(v)[:80])
