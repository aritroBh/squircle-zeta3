import json, re
tr = open(SRC+'transcript.json').read()
print("transcript length:", len(tr))
# Find occurrences of key strings
for pat in ['3.2688325650667', '3.5164954371', 'q_B', '1.0757649304908']:
    idxs = [m.start() for m in re.finditer(re.escape(pat), tr)]
    print(f"\nPattern {pat!r}: {len(idxs)} occurrences at {idxs[:5]}")
    if idxs:
        i = idxs[0]
        print("  context:", tr[max(0,i-200):i+200].replace('\\n','\n'))
