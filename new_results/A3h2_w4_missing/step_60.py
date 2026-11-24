import pickle
d = pickle.load(open('cands29_reduced.pickle','rb'))
print(list(d.keys()))
for k,v in d.items():
    print(k, type(v), (len(v) if hasattr(v,'__len__') else v))
