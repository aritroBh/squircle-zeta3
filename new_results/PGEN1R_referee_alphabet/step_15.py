import os
open('progress.log','a').write('PGEN1R start: read PGEN1 result.md, paper subsection, DEC2 result.md\n')
# explore data sources for exact closed forms
for d in ['../DEC2_decagon_exact_assembly','../OCT6_octagon_closed','../OCT5_octagon_exact',
          '../DOD5_dodecagon_termwise','../OCT9c_octagon_close','../OCT9d_octagon_short','../OCT9b_octagon_debug']:
    if os.path.exists(d):
        print(d, sorted(os.listdir(d))[:25])
    else:
        print(d, 'MISSING')
