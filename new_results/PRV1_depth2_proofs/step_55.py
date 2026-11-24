import os
for d in ['DS3_hexagon_order5','K1f_hexagon_order5','K1e_hexagon_closed']:
    p = f'research/08_theory/{d}/result.md'
    if os.path.exists(p):
        with open(p) as f: txt = f.read()
        import re
        hits = [l for l in txt.split('\n') if 'Li_21' in l or 'Li_{2,1}' in l or 'Landen' in l][:12]
        print('='*20, d)
        print('\n'.join(hits))
