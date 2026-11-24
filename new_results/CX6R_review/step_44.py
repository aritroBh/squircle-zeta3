import subprocess, os
for d in ('D13s_fix_recertify','D13u_close_h3'):
    p=f'<project>/research/08_theory/{d}'
    print("==== ", d, sorted(os.listdir(p))[:20])
    for cand in ('result.md','RESULT.md','summary.md'):
        fp = os.path.join(p, cand)
        if os.path.exists(fp):
            txt = open(fp).read()
            import re
            hits = [l for l in txt.splitlines() if '254' in l][:6]
            print("\n".join(hits[:6]))
            break
