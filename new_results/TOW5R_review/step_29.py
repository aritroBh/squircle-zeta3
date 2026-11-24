import subprocess, os
t5 = '<project>/research/08_theory/TOW5_tower_depth10'
base = '<project>/research/08_theory'
pairs = [
 (t5+'/tow2lib.py', base+'/TOW2_tower_all_depths/tow2lib.py'),
 (t5+'/tow3lib.py', base+'/TOW3_tower_depth6/tow3lib.py'),
 (t5+'/tow4lib.py', base+'/TOW4_tower_depth8/tow4lib.py'),
 (t5+'/fastred.py', base+'/WTS_weight_k_structure/theory/fastred.py'),
 (t5+'/prove_S34.py', base+'/TOW2_tower_all_depths/prove_S34.py'),
 (t5+'/check_data.py', base+'/TOW3_tower_depth6/check_data.py'),
 (t5+'/L_simplified.pkl', base+'/TOW2_tower_all_depths/L_simplified.pkl'),
 (t5+'/L_parsed.pkl', base+'/TOW3_tower_depth6/L_parsed.pkl'),
 (t5+'/arho_simplified.pkl', base+'/TOW2_tower_all_depths/arho_simplified.pkl'),
 (t5+'/forced_F78.pkl', base+'/TOW4_tower_depth8/forced_F78.pkl'),
]
for a, b in pairs:
    if not os.path.exists(b):
        print('MISSING', b); continue
    r = subprocess.run(['cmp', '-s', a, b])
    print(('SAME ' if r.returncode == 0 else 'DIFF '), os.path.basename(a), 'vs', b.split('/')[-2])
# find TOWN1 trunc engine
for f in os.listdir(base+'/TOWN1_tower_large_n'):
    if 'trunc' in f or 'engine' in f: print('TOWN1 has', f)
r = subprocess.run(['grep','-rl','trunc_engine', base+'/TOWN1_tower_large_n/'], capture_output=True, text=True)
print(r.stdout)
