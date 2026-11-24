
import pickle
fs = pickle.load(open('../COR8_depth_two_any_angle/full_sols_v3.pkl','rb'))
sols = fs['sols']; alln = fs['alln']
print("sols keys:", list(sols.keys())[:10])
print("alln len:", len(alln), "first items:", alln[:5])
# inspect structure of one solution
for k in list(sols.keys())[:3]:
    v = sols[k]
    print("key", k, "type", type(v), "len", len(v) if hasattr(v,'__len__') else '?')
    print(str(v)[:500])
