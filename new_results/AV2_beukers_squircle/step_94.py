import os, json
os.chdir('<project>/research/08_theory/AV2_beukers_squircle')
print(sorted(os.listdir('.')))
print()
# consolidate B3 values with error estimates
b3f = json.load(open('b3_final.json'))
b3h = json.load(open('b3_hifi.json'))
for p in ['1','1.5','2','3','4','6','8','12','16','24','32','48','64','256','1024','4096','16384','262144','1048576']:
    src = b3h.get(p) or b3f.get(p)
    if src:
        print(p, src)
