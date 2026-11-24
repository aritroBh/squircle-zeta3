import json, mpmath as mp
mp.mp.dps = 1800
m5b='<project>/research/08_theory/M5b_w6_closed_form'
def load(t,b):
    return json.load(open(f'{m5b}/{t}_{b}.json'))
w656, w658 = load('w6',5600), load('w6',5800)
wp46, wp48 = load('w5pp',4600), load('w5pp',4800)
d6 = abs(mp.mpf(w656['value_full'])-mp.mpf(w658['value_full']))
d5 = abs(mp.mpf(wp46['value_full'])-mp.mpf(wp48['value_full']))
print('W6 5600 vs 5800 full-value absdiff:', mp.nstr(d6,6))
print('W5pp 4600 vs 4800 full-value absdiff:', mp.nstr(d5,6))
# canonical W6_5400.json (>=1600 digits): use 5800 run value truncated to 1650 digits
v = mp.mpf(w658['value_full'])
s1650 = mp.nstr(v, 1650)
json.dump({'source':'m5b_precision.py w6 5600+5800 (agree to 1e-1715)','digits':1650,'value':s1650,
           'radius_5800': w658['radius'], 'tail_5800': w658['tail_bound']},
          open(m5b+'/W6_5400.json','w'), indent=1)
print('saved W6_5400.json with 1650 digits')
s1350 = mp.nstr(mp.mpf(wp48['value_full']), 1350)
json.dump({'source':'m5b_precision.py w5pp 4600+4800 (agree)','digits':1350,'value':s1350},
          open(m5b+'/W5pp_1350.json','w'), indent=1)
print('saved W5pp_1350.json')
with open(m5b+'/progress.log','a') as f:
    f.write(f'2026-10-27T16:45 W6 5600/5800 agree {mp.nstr(d6,6)}; W5pp 4600/4800 agree {mp.nstr(d5,6)}; W6_5400.json + W5pp_1350.json saved\n')
