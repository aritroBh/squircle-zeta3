import pickle, os
wts = '<project>/research/08_theory/WTS_weight_k_structure/theory'
for k in [6,7,8]:
    K = pickle.load(open(f'{wts}/kern_{k}.pkl','rb'))
    print(f"--- kern_{k}.pkl: keys = {sorted(K.keys())}")
    sample_key = sorted(K.keys())[-1]
    print("   sample key", sample_key, "->", {kk: (type(v).__name__) for kk,v in K[sample_key].items()})
    # print structure of one entry
    e = K[sample_key]
    print("   'real' poles:", list(e['real'].keys())[:8], "..." if len(e['real'])>8 else "")
    print("   'ipos':", e['ipos'], " 'ineg':", e['ineg'])
