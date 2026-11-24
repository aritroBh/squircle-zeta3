
import pickle
sols2 = pickle.load(open('classical_sols_tracked.pkl','rb'))
seen = pickle.load(open('li_arg_pairs.pkl','rb'))
# inverse maps: (u,x) -> arg string for Li's; (u,) -> arg for Li1
gx = {}
g1 = {}
for (u,x), arg in seen.items():
    gx.setdefault((u,x), arg)
for (u,x), arg in seen.items():
    g1.setdefault((u,), arg)
def fmt_name(name):
    kind = name[0]
    if kind=='Li4': return f"Li4({gx.get(name[1], name[1])})"
    if kind=='L3L1': return f"Li3({gx.get(name[1], name[1])})*Li1({g1.get(name[2], name[2])})"
    if kind=='L2L2': return f"Li2({gx.get(name[1], name[1])})*Li2({gx.get(name[2], name[2])})"
    if kind=='L2L1L1': return f"Li2({gx.get(name[1], name[1])})*Li1({g1.get(name[2], name[2])})*Li1({g1.get(name[3], name[3])})"
    if kind=='L1^4': return "*".join(f"Li1({g1.get(k,k)})" for k in name[1])
lines = []
for e in (-1,0,1):
    lines.append(f"\n=== coefficient of c^({e}) ===")
    for name, cf in sorted(sols2[e].items(), key=lambda kv: str(kv[0])):
        lines.append(f"  {cf}  {fmt_name(name)}")
txt = "\n".join(lines)
open('classical_identity.txt','w').write(txt)
print(txt[:4500])
