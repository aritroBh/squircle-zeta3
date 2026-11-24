base = '<project>/research/08_theory'
import json, os
for f in ['K1c_engine_predict/geometry.json', 'K1e_hexagon_closed/geometry_fresh.json', 'K1d_engine_fixed/geometry.json']:
    p = f'{base}/{f}'
    if os.path.exists(p):
        print('===', f)
        print(open(p).read()[:1500])
