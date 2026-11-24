import os
for d in ['K1d_engine_fixed','K1_nonlocal_law','K1b_nonlocal_finish']:
    p = 'research/08_theory/%s/result.md'%d
    if os.path.exists(p):
        c = open(p).read()
        print("="*75)
        print("###", d, "len", len(c))
        print(c[:5500])
