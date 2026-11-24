d='research/08_theory/D13s_fix_recertify/D13s_work/fixed/'
for f in ['step_18.py','planner_sweep_ratio.py']:
    print('='*30, f)
    print(open(d+f).read())
print('='*30, 'D13r result.md')
print(open('research/08_theory/D13r_lemma_review/result.md').read()[:4000])
