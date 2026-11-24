
WORK = 'research/08_theory/D13s_fix_recertify/D13s_work'
with open(WORK + '/D13g_all_N__step_38.py') as f:
    lines = f.readlines()
for i in range(95, 156):
    print(f"{i+1:4d}: {lines[i]}", end='')
