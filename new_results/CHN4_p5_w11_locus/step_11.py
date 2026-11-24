import pickle
d = pickle.load(open('test_pt_w1112.pkl','rb'))
print(type(d))
if isinstance(d, dict):
    for k,v in d.items():
        print("key:", k, "type:", type(v))
        print("   ", v)
else:
    print(d)
