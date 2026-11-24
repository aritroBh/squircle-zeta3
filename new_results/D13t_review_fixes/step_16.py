d='research/08_theory/D13s_fix_recertify/D13s_work/fixed/'
for f in ['verify_TM_final.py','verify_TM_mono.py','verify_theorem_truth.py']:
    print('='*30, f)
    print(open(d+f).read())
print('='*30,'D13s progress.log')
print(open('research/08_theory/D13s_fix_recertify/progress.log').read())
