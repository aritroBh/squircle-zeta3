import os, json
os.chdir('<project>/research/08_theory/AV2_beukers_squircle')
from mpmath import mp, mpf, nstr
mp.dps = 120
b2 = json.load(open('b2_table.json'))
b3f = json.load(open('b3_final.json'))
b3h = json.load(open('b3_hifi.json'))
print("=== B2 table (selected, 45 digits) ===")
for p in ['1','1.5','2','3','4','6','8','12','16','24','32','48','64','96','128','256','1024','4096']:
    e = b2[p]
    print("p=%5s  %s   (err %s)" % (p, mp.nstr(mpf(e['B2']), 45), e['err']))
print()
print("=== B3 table with achieved precision ===")
for p in ['1','1.5','2','3','4','6','8','12','16','24','32','48','64','256','1024','4096','16384','262144','1048576']:
    e = b3h.get(p) or b3f.get(p)
    if e:
        d = e.get('d') or e.get('d_fine_coarse')
        dg = e.get('d_grid')
        print("p=%8s  %s  (d=%s, dgrid=%s)" % (p, mp.nstr(mpf(e['B3']), 45), d, dg))
