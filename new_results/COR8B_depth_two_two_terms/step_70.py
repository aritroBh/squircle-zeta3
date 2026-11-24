
import pickle
for f in ['dt_syms.pkl','dt_syms_2kept.pkl','symbols_2kept.pkl']:
    d = pickle.load(open('../COR8_depth_two_any_angle/'+f,'rb'))
    print("="*60); print(f, type(d))
    if isinstance(d, dict):
        for k in list(d.keys())[:8]:
            print("  KEY:", k, type(d[k]))
    elif isinstance(d, list):
        print("len", len(d), str(d[:3])[:400])
