import os, shutil, json
base = '<project>/research/08_theory'
c = base + '/ZQ5c_zagier_w5_capped'
d = base + '/ZQ5d_zagier_d13_clean'
os.makedirs(d+'/out/d13', exist_ok=True)
os.makedirs(d+'/planner', exist_ok=True)
# copy libs (COPY, never edit originals)
for f in ['w5blib.py','w5pipe.py','w5verify.py']:
    shutil.copy(c+'/'+f, d+'/'+f)
shutil.copy(base+'/ZQ5_zagier_weight5/planner/planner_li5.py', d+'/planner_li5.py')
# copy data
for f in ['candidates_cap.json','kernel5_cap.json','genuine_lattice_cap.json','candidates.json',
          'vals_cap_2000.pkl','atoms5_cap_d17.pkl','lll_capd17_P1500.json','lll_capd17_P1800.json']:
    shutil.copy(c+'/out/d13/'+f, d+'/out/d13/'+f)
shutil.copy(c+'/planner/d13_row8.json', d+'/planner/d13_row8.json')
# quick sanity: full pool size, heights
J = json.load(open(d+'/out/d13/candidates.json'))
print('full d13 pool:', J.keys() if isinstance(J,dict) else len(J))
if isinstance(J, dict):
    print('meta:', {k:J[k] for k in J if k!='cands'})
    from collections import Counter
    hs = Counter(x['h'] for x in J['cands'])
    print('N cands:', len(J['cands']), 'heights:', dict(sorted(hs.items())))
with open(d+'/progress.log','a') as f:
    f.write('setup: copied w5blib/w5pipe/w5verify/planner_li5 + d13 data from ZQ5c (cap275 pool, kernel, genuine lattice, vals@2000, full 487 pool)\n')
print('files in workdir:', sorted(os.listdir(d)))
