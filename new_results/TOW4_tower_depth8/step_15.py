import os, re
base = '<project>/research/08_theory'
with open(base + '/CX8_tow2_review/report.md') as f:
    txt = f.read()
# find sections about master equation and x-indexing
for kw in ['master equation', 'x-index', 'x index', '7.', 'item 7', 'write-up corrections']:
    idxs = [m.start() for m in re.finditer(re.escape(kw), txt)]
    print(kw, idxs[:5])
print(txt[7000:13000])
