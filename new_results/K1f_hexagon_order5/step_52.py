base = '<project>/research/08_theory'
k1d = f'{base}/K1d_engine_fixed'
for s in [20, 26, 27]:
    src = open(f'{k1d}/step_{s:02d}.py').read()
    print('='*25, 'step_%02d'%s, '(%d chars)'%len(src), '='*25)
    print(src[:3500])
