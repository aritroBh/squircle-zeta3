import re
d='research/08_theory/D13s_fix_recertify/D13s_work/fixed/'
for fn in ['step_10_J192.py','step_38.py']:
    print('='*20, fn)
    for i,l in enumerate(open(d+fn).read().split('\n'),1):
        if re.search(r'rt2\s*/|\.sqrt\(\)|\*\*\s*arb|\*\*\s*arb|2\*\*arb|\*\*"1\.5"|\*\*arb\("', l) or 'ecw' in l or 'arb(2)/pi' in l:
            print(f"{i:4d}: {l.rstrip()}")
