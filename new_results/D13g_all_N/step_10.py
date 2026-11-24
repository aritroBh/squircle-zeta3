
import os
BASE='research/08_theory'
lines = open(os.path.join(BASE,'D13f_sweep','mk_deg9.tsv')).read().strip().split('\n')
print("rows:", len(lines), "first idx:", lines[0][:20], "last idx:", lines[-1][:30])
# parse: k, Mkmain arb string, third column
import re
rows = []
for ln in lines:
    parts = ln.split('\t')
    rows.append((int(parts[0]), parts[1], parts[2] if len(parts)>2 else None))
print("cols in row 5:", len(lines[5].split('\t')))
for k in [99, 100, 150, 250, 399]:
    r = rows[k] if rows[k][0]==k else [r for r in rows if r[0]==k][0]
    print(r[0], r[1][:60], '|', (r[2] or '')[:40])
