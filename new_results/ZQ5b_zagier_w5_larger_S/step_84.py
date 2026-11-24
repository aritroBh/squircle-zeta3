import os, json, subprocess, time
from fractions import Fraction as Ff
os.chdir('<project>/research/08_theory/ZQ5b_zagier_w5_larger_S')
os.system("pkill -f step_11b_d13_lll.py"); time.sleep(1)
# check D=17 planner relation arg heights vs candidate pool
ZQ='../ZQ5_zagier_weight5'
J = json.load(open(f'{ZQ}/out/d17/candidates.json'))
R = json.load(open(f'{ZQ}/planner/d17_relation_planner.json'))
h_of = {}
for c in J['cands']:
    key = f"{c['x'][0]}+{c['x'][1]}*sqrt17"
    h_of[key] = max(abs(e) for e in c['vx'])
used_hs = sorted(h_of.get(k, '?') for k in R['li5'].keys())
import collections
print("D=17 used-arg heights:", dict(sorted(collections.Counter(used_hs).items())))
print("max used h:", max(used_hs), " (pool B was 8, h dist had up to 8)")
# also for my final form
F17 = json.load(open('d17_final.json'))
used_hs2 = sorted(h_of.get(k, '?') for k in F17['li5'].keys())
print("D=17 final-form arg heights:", dict(sorted(collections.Counter(used_hs2).items())))
