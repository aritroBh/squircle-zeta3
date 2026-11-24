import os
base = '<project>'
here = base + '/research/08_theory/M4b2_w4pp_nice'
print('work dir contents:', sorted(os.listdir(here)))
print()
log = open(base+'/research/08_theory/M4b_w4pp_exact/progress.log').read()
# print first and key timing lines
for ln in log.splitlines():
    if any(k in ln for k in ['step0','solve done','EXACT','saved','rank','CRT','s=']):
        print(ln[:150])
