# Planner: print the exact zero-relations (no T) found with the extra atoms, and which extra atoms they use.
import sys, runpy
sys.argv = ['x', '700', '700']
g = runpy.run_path('planner_d21_extra_atoms.py', run_name='planner_d21')
withT, good = g['r1']
names = g['base_names']
for a, res in good:
    used = [(names[k], a[k]) for k in range(len(a)) if a[k]]
    ex = [u for u in used if u[0].startswith(('z3', 'T3'))]
    print('maxcoef', max(map(abs, a)), 'nterms', len(used), 'residual', g['mp'].nstr(res, 3))
    print('   extra atoms used:', ex)
    print('   kernel combos used:', [u for u in used if u[0].startswith('K')])
