import pickle
d = pickle.load(open("test_indep_450.pkl","rb"))
print(type(d), list(d.keys()) if isinstance(d,dict) else len(d))
s = repr(d)
print(s[:4000])
