import pickle
base='<project>'
d=pickle.load(open(base+'/research/08_theory/A3f4_glanois_basis/five_relations.pickle','rb'))
print(type(d))
if isinstance(d,dict): 
    for k,v in d.items(): print(k, type(v), str(v)[:300])
elif isinstance(d,(list,tuple)):
    print(len(d)); print(str(d[0])[:500])
