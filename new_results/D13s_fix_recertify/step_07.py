
WORK = 'research/08_theory/D13s_fix_recertify/D13s_work'
for fn in ['D13h_final_cell__step_18.py', 'D13h_final_cell__planner_sweep_ratio.py']:
    print("="*30, fn, "="*30)
    with open(WORK + '/' + fn) as f:
        print(f.read())
