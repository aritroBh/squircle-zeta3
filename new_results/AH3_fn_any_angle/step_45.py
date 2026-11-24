
import os
os.chdir('<project>/research/08_theory/AH3_fn_any_angle')
print(os.listdir('.'))
print()
md = open('result.md').read()
# quick sanity: count sections and check for placeholders
for key in ['Task 1', 'Task 2', 'Task 3', 'Task 4', 'Task 5', 'Evidence table', 'cubic']:
    print(key, '->', md.count(key))
