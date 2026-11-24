
import os
WORK = 'research/08_theory/D13s_fix_recertify/D13s_work'

def show_lines(path, a, b, label):
    print(f"===== {label}: {os.path.basename(path)} lines {a}-{b} =====")
    with open(path) as f:
        lines = f.readlines()
    for i in range(a-1, min(b, len(lines))):
        print(f"{i+1:4d}: {lines[i]}", end='')
    print()

# B1: final cell, line ~168
show_lines(os.path.join(WORK,'D13h_final_cell__step_10_J192.py'), 160, 180, "B1 candidate")
# B2: step_38.py, line ~165
show_lines(os.path.join(WORK,'D13g_all_N__step_38.py'), 155, 180, "B2 candidate")
