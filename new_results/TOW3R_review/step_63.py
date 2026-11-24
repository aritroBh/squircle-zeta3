import hashlib, os
base = '<project>/research/08_theory/'
d2 = base+'TOW2_tower_all_depths/'
d3 = base+'TOW3_tower_depth6/'
def md5(p):
    return hashlib.md5(open(p,'rb').read()).hexdigest()
for f in ['tow2lib.py','L_simplified.pkl','arho_simplified.pkl','e_fits.pkl','prove_S12.py','prove_S34.py','prove_tops.py','prove_adj34.py','prove_E3b34.py','step_48.py','step_49.py','step_52.py']:
    p2, p3 = d2+f, d3+f
    e2, e3 = os.path.exists(p2), os.path.exists(p3)
    same = md5(p2)==md5(p3) if (e2 and e3) else None
    print(f'{f:28s} TOW2:{e2} TOW3:{e3} identical:{same}')
