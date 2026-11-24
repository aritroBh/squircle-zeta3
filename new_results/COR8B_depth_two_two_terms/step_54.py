
import pickle
try:
    fs = pickle.load(open('../COR8_depth_two_any_angle/full_sols_v3.pkl','rb'))
    print(type(fs))
    if isinstance(fs, dict):
        for k,v in list(fs.items())[:10]:
            print(k, type(v))
    elif isinstance(fs, (list,tuple)):
        print("len", len(fs))
        for x in fs[:5]: print(type(x), str(x)[:200])
except Exception as e:
    print("ERR", e)
