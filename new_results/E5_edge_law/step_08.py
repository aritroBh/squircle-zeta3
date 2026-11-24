import os, json
d17='research/08_theory/D17_zeta5_polygons'
print(open(os.path.join(d17,'closed_form_J5.json')).read()[:3000])
print("="*80)
print(open(os.path.join(d17,'J5_of_h.txt')).read()[:3000])
