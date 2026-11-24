import shutil, os, time
work = '<project>/research/08_theory/TOW5_tower_depth10'
tow4 = '<project>/research/08_theory/TOW4_tower_depth8'
tow2 = '<project>/research/08_theory/TOW2_tower_all_depths'
town1 = '<project>/research/08_theory/TOWN1_tower_large_n'
wts = '<project>/research/08_theory/WTS_weight_k_structure/theory'

copies = [
    (tow4, 'tow2lib.py'), (tow4, 'tow3lib.py'), (tow4, 'tow4lib.py'),
    (tow4, 'L_simplified.pkl'), (tow4, 'L_parsed.pkl'), (tow4, 'check_data.py'),
    (tow4, 'arho_simplified.pkl'),
    (tow2, 'prove_S34.py'),
    (town1, 'trunc_engine.py'), (town1, 'driver.py'),
    (wts, 'fastred.py'),
]
for src, f in copies:
    shutil.copy(os.path.join(src, f), os.path.join(work, f))
    print('copied', f)

# Also copy kernel data read-only? No - read from wts directly.
# Log progress
with open(os.path.join(work,'progress.log'),'a') as f:
    f.write(f"{time.strftime('%H:%M:%S')} setup: copied tow2lib/tow3lib/tow4lib/L caches/check_data/prove_S34/trunc_engine/driver/fastred from TOW4,TOW2,TOWN1,WTS\n")
print("done")
