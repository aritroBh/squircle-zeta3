import os
base = os.path.dirname(os.getcwd())
for f in ['OCT8_symbol_test','DOD2_dodecagon_finish','DOD1_dodecagon']:
    p = os.path.join(base,f)
    print('='*20, f)
    for x in sorted(os.listdir(p)):
        print('   ', x)
