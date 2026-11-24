import os
for f in sorted(os.listdir('.')):
    print(f, os.path.getsize(f))
print("\n--- head of step_13.py (tight-vs-coded table) ---")
print('\n'.join(open('step_13.py').read().split('\n')[:6]))
print("\n--- head of step_14.py (verbatim step_18 rerun) ---")
print('\n'.join(open('step_14.py').read().split('\n')[:6]))
