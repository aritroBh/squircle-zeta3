import pickle
d = pickle.load(open('cands29_reduced.pickle','rb'))
fw = d['freewords']; cands = d['cands']
print(fw[:5])
c0 = cands[0]
print(type(c0))
if isinstance(c0, dict):
    print(list(c0.keys()))
    for k,v in c0.items():
        print(k, type(v), str(v)[:200])
