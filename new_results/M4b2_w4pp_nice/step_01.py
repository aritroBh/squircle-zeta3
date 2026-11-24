import os
base = '<project>'
# list relevant dirs
for d in ['research/08_theory/M4b_w4pp_exact', 'research/08_theory/M4a_decomposition']:
    p = os.path.join(base, d)
    print('===', d, os.path.isdir(p))
    if os.path.isdir(p):
        for f in sorted(os.listdir(p)):
            print('   ', f, os.path.getsize(os.path.join(p,f)))
