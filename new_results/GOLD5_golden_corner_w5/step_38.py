import json
d = json.load(open('../G2_M4alpha/flat_exact_N640_K4.json'))
for k in d:
    sub = d[k]
    ks = list(sub.keys())
    print(k, "n-range:", ks[:3], "...", ks[-3:], "count:", len(ks))
    # sample
    n0 = ks[0]
    print("   sample n=", n0, ":", str(sub[n0])[:200])
f = open('../G2_M4alpha/flat_exact_N640_K4.out').read()
print(f[:3000])
