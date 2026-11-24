
import os, json
ROOT = '<project>'
M5A = os.path.join(ROOT, 'research/08_theory/M5a_next_coefficient')
ap = json.load(open(os.path.join(M5A, 'atom_pairings.json')))
# print the ones for our 4 atoms
for kk in ['Li_33(b,1)', 'Li_33(1/sqrt2,1)', 'Li_33(b,2)', 'Li_42(b,1/sqrt2)']:
    print(kk, ap.get(kk))
print("\nAll keys with Li_33 or Li_42 at b or 1/sqrt2:")
for kk, vv in ap.items():
    if ('Li_33' in kk or 'Li_42' in kk):
        print(" ", kk, vv)
