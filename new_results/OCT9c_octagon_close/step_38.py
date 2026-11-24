import pickle
p3 = pickle.load(open('../OCT9b_octagon_debug/P3_assembled.pkl','rb'))
print(type(p3))
for k,v in (p3.items() if isinstance(p3,dict) else enumerate(p3)):
    print(k, str(v)[:110])
