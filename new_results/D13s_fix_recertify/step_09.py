
WORK = 'research/08_theory/D13s_fix_recertify/D13s_work'
with open(WORK + '/D13h_final_cell__step_10_J192.py') as f:
    lines = f.readlines()
print("total lines:", len(lines))
for i in range(88, 160):
    print(f"{i+1:4d}: {lines[i]}", end='')
