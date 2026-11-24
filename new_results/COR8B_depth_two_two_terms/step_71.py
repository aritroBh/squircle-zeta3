
import pickle
s2 = pickle.load(open('../COR8_depth_two_any_angle/symbols_2kept.pkl','rb'))
for k in s2:
    v = s2[k]
    print("KEY", k, "tensors:", len(v))
    for kk in list(v.keys())[:6]:
        print("   ", kk, v[kk])
