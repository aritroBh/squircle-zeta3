import re
paper = '<project>/paper/ieee/main.tex'
txt = open(paper).read()
print("LEN:", len(txt))
# Find "Diamond End" section
i = txt.find('Diamond End')
print(txt[i-200:i+200] if i>=0 else "not found 'Diamond End'")
# search for delta(N_3), (2,1), invariant
for pat in [r'delta\(N', r'\\delta\(N', r'\(2,1\)', r'invariant']:
    idxs=[m.start() for m in re.finditer(pat, txt)]
    print(f"pattern {pat!r}: {len(idxs)} hits at {idxs[:8]}")
