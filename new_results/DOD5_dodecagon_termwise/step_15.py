import pickle
base = "<project>/research/08_theory/"
D1, D2, D3 = base+"DOD1_dodecagon/", base+"DOD2_dodecagon_finish/", base+"DOD3_dodecagon_symbol/"
rl = pickle.load(open(D3+"R_li3_content.pkl","rb"))
print("R_li3_content type:", type(rl), rl.keys() if hasattr(rl,'keys') else len(rl))
import json
def show(o, depth=0, maxd=2):
    pre = "  "*depth
    if isinstance(o, dict):
        for k in list(o)[:4]:
            print(pre, repr(k)[:60], "->", type(o[k]))
            if depth<maxd: show(o[k], depth+1, maxd)
    elif isinstance(o, (list,tuple)):
        print(pre, type(o), "len", len(o), "first:", repr(o[0])[:120] if o else "")
    else:
        print(pre, repr(o)[:200])
show(rl)
rv = pickle.load(open(D3+"R_all12_vals620.pkl","rb"))
print("\nR_all12_vals620:", type(rv))
show(rv, maxd=1)
rs = pickle.load(open(D3+"R_all12_sym_str.pkl","rb"))
print("\nR_all12_sym_str:", type(rs))
show(rs, maxd=1)
