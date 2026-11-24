import json, os
import mpmath as mp
work='<project>/research/08_theory/WT4_weight4_all_k'
mp.mp.dps = 900
res = json.load(open(os.path.join(work,'hp_results.json')))
for tgt in ['D1','D2']:
    v9 = mp.mpf(res[f'{tgt}@900']['value'])
    v8 = mp.mpf(res[f'{tgt}@850']['value'])
    d = abs(v9-v8)
    print(f"{tgt}: |@900 - @850| = {mp.nstr(d,4)}  -> {mp.nstr(-mp.log10(d),1)} digits agreement (mp.mp.dps={mp.mp.dps})")
    print(f"  {tgt}(860 digits) = {res[f'{tgt}@900']['value'][:100]}...")
