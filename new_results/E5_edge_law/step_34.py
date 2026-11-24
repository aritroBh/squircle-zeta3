import os
print(sorted(os.listdir('.')))
print(sorted(os.listdir('steps')))
if os.path.exists('step_33.py'):
    print(open('step_33.py').read()[:400])
