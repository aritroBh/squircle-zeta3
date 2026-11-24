import os
for d in ['DOD1_dodecagon','DOD2_dodecagon_finish','DOD3_dodecagon_symbol','DOD4_dodecagon_closed']:
    p = '../'+d
    print('=====', d)
    print([f for f in sorted(os.listdir(p)) if not f.startswith('step_') and not f.startswith('transcript')][:30])
