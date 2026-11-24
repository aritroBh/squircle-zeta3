import os, sys, pickle
Z3 = 'Zeta' + '(3)'
ROOT = os.path.join('<home>/Downloads/Research', Z3)
WORK = os.path.join(ROOT, 'research/08_theory/M4b_w4pp_exact')
os.chdir(WORK)
sys.path.insert(0, WORK)
import m4b_common as mc
rows = mc.load_T0()
from collections import Counter
fam = Counter(l[1][0] for l, _ in rows)
print('T0 families:', dict(fam))
# sample labels of each family
seen = {}
for l, _ in rows:
    if l[1][0] not in seen:
        seen[l[1][0]] = l
for k, v in seen.items():
    print(k, v)
# check y nonzero counts vs S
Y = pickle.load(open('stage2_y.pkl','rb'))
S = pickle.load(open(os.path.join(ROOT,'research/08_theory/A3h5_exact_certificate/stage2_rows_S.pkl'),'rb'))['S']
used = [j for j in range(5042) if Y['rat'][j] or Y['sq'][j]]
print('used rows:', len(used))
