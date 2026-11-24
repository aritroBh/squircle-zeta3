import os, pickle, json
from fractions import Fraction as Fr
os.chdir('<project>/research/08_theory/LEAD1_leading_genfun')
base = '<project>'
theory = os.path.join(base, 'research', '08_theory', 'WTS_weight_k_structure', 'theory')
kern = {}
for K in range(1, 19):
    fn = os.path.join(theory, f'kern_{K}.pkl')
    with open(fn,'rb') as f: kern[K] = pickle.load(f)
POLE = {(0,0):"0",(1,1):"a",(-1,-1):"-a",(-1,1):"b",(1,-1):"-b"}
def c(k,key,r):
    P = kern.get(k,{}).get(key)
    if P is None or r not in P['real']: return (Fr(0),Fr(0))
    z = P['real'][r]; return (z[0],z[1])
out = {}
for d in range(1,9):
    keys = sorted({key for K in kern for key in kern[K] if sum(key)==d}, key=lambda kk:-kk[0])
    for key in keys:
        poles = sorted({r for K in kern for r in kern[K].get(key,{'real':{}})['real']})
        rec = {}
        for r in poles:
            seq = [c(K,key,r) for K in range(d+1,19)]
            D=[]; s=list(seq)
            for m in range(1,d+1):
                D.append(s[0]); s=[(a2[0]-a1[0],a2[1]-a1[1]) for a1,a2 in zip(s,s[1:])]
            ok = all(t==(Fr(0),Fr(0)) for t in s)
            rec[POLE[r]] = {"Dtop":[str(D[d-1][0]),str(D[d-1][1])], "poly_ok": ok, "npts": len(seq)}
        out[f"{d}|{key[0]},{key[1]}"] = rec
json.dump(out, open('Dtop.json','w'), indent=1)
# cross-check d<=7 and 8|(8,0) against genfun.json
gj = json.load(open(os.path.join(theory,'genfun.json')))
bad = []
for name, rec in gj.items():
    d = int(name.split('|')[0])
    for pole, pr in rec.items():
        mine = out[name][pole]['Dtop']
        theirs = pr['D'][d-1]
        if mine != theirs: bad.append((name,pole,mine,theirs))
print("mismatches vs genfun.json:", bad)
print("all poly_ok:", all(pr['poly_ok'] for rec in out.values() for pr in rec.values()))
print("keys:", list(out.keys()))
with open('progress.log','a') as f:
    f.write('step_08: extracted top Newton coeffs D_d (d=1..8, all keys, k=d+1..18) -> Dtop.json; matches genfun.json for d<=7 + 8|(8,0)\n')
