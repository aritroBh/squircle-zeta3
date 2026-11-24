import pickle
d = pickle.load(open('cands29_reduced.pickle','rb'))
print(type(d))
if isinstance(d, list):
    print(len(d))
    print(type(d[0]))
    print(d[0] if not isinstance(d[0], dict) else list(d[0].keys()))
