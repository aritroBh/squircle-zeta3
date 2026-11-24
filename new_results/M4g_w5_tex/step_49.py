import json
pl=json.load(open('atoms_planner.json'))
for f in ('L','lam','Li2(b)','Li5((-1+1*sqrt2))'):
    print(f, '->', pl[f][:80])
print('N5 ->', pl['N5'][:60])
print('Lchi ->', pl['Lchi8_5'][:60])
import mpmath as mp
mp.dps=470
print('direct L   :', mp.nstr(mp.log(1+mp.sqrt(2)),60))
print('direct Li2b:', mp.nstr(mp.polylog(2,mp.sqrt(2)-1),60))
