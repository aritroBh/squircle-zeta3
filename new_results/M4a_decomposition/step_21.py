import os, re
BASE = '<project>'
rm = open(os.path.join(BASE,'research/08_theory/A3b_certified2/prev_result.md')).read()
for pat in ['321','10^-3','acb','tanh','double']:
    idxs = [m.start() for m in re.finditer(re.escape(pat), rm)]
    print(pat, idxs[:8])
print('='*40)
print(rm[:4000])
