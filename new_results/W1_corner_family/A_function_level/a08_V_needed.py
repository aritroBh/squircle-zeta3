"""a08_V_needed.py -- which non-nice Li4 arguments are forced in the classical part of V - Theta_V
(symbol mod products).  Reuses a07's setup."""
import sys
exec(open('a07_V_feasibility.py').read().split("pickle.dump({'adm'")[0])
import sympy as sp
non = [i for i in range(len(cols)) if i not in nice_cols]
print('non-nice columns:', len(non))
r0 = Mn.rank()
need = []
for i in non:
    Mi = M[:, nice_cols + [i]]
    if Mi.row_join(bvec).rank() == Mi.rank():
        need.append(i)
print('single non-nice columns that complete the nice set:')
for i in need:
    print('   ', cols[i][0], ' -> value at 1/sqrt2:', cols[i][1], ' (numeric %.6f)' % float(cols[i][1]))
# which non-nice point VALUES occur at all among admissible args with value in (-1,1)
vals = {}
for i in non:
    v = cols[i][1]
    vals.setdefault(v, []).append(str(cols[i][0]))
print('non-nice point values present:')
for v, fs in sorted(vals.items(), key=lambda t: float(t[0])):
    print('   %-28s %10.6f  %d args' % (v, float(v), len(fs)))
