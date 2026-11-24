import json, re
SRC = 'research/08_theory/K4c_3d_edge_law/'
tr = open(SRC+'transcript.json').read()
print("transcript length:", len(tr))
for pat in ['3.2688325650667', '3.5164954371452986', 'q_B', '1.0757649304908']:
    idxs = [m.start() for m in re.finditer(re.escape(pat), tr)]
    print(f"\nPattern {pat!r}: {len(idxs)} occ at {idxs[:3]}")
    if idxs:
        i = idxs[0]
        print("  ctx:", tr[max(0,i-150):i+120].replace('\\n','\n'))
