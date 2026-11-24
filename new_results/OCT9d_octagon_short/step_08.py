import pickle, sys
base = "<project>/research/08_theory/"
O9B = base+"OCT9b_octagon_debug"
sys.path.insert(0, O9B)
s123 = pickle.load(open(O9B+'/sol_B123.pkl','rb'))
print(type(s123), list(s123.keys()) if isinstance(s123, dict) else len(s123))
for k in s123:
    print("edge", k, "keys:", list(s123[k].keys()))
    for kk, vv in s123[k].items():
        if kk == 'sol':
            print("  sol: len", len(vv), "first few:", vv[:5])
        elif kk == 'Wlo':
            print("  Wlo: type", type(vv), "str head:", str(vv)[:300])
        else:
            print(" ", kk, "=", str(vv)[:300])
