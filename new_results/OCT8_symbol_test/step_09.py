import os, pickle
os.chdir("<project>/research/08_theory/OCT7_generic_octagon")
d = pickle.load(open("test_indep_450.pkl","rb"))
print(type(d), list(d.keys()) if isinstance(d,dict) else len(d))
print(repr(d)[:4000])
