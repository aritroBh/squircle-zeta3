with open('../PGEN1_polygon_alphabet/result.md') as f:
    txt = f.read()
print(txt[7500:])
print('='*80)
import os
print(sorted(os.listdir('../PGEN1_polygon_alphabet'))[:40])
print(sorted(os.listdir('../DEC2_decagon_exact_assembly')))
