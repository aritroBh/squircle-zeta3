import os
dst = '<project>/research/08_theory/M4i_w5_noL'
for n in range(16, 36):
    p = f'{dst}/step_{n:02d}.py'
    if os.path.exists(p):
        head = open(p).read().strip().split('\n')[0][:95]
        # also grab a identifying line
        src = open(p).read()
        tag = ''
        for key in ['lrel_build','lrel_verify','nol_build','nol_optimize2','nol_export','verify_json','nol_cross','deps_consistency','drop_L','lrel_smoke','progress.log']:
            if key in src: tag += key+';'
        print(f'step_{n:02d}: {head}   [{tag}]')
