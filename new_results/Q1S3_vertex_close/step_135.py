import pickle
names=pickle.load(open('q3_names.pkl','rb'))['names']
bad=[(i,type(v)) for i,(k,v) in enumerate(names) if not isinstance(v,tuple)]
print("non-tuple values:", bad[:5], " total:", len(bad))
print("entry 272:", names[272], " entry 300:", names[300][:2])
