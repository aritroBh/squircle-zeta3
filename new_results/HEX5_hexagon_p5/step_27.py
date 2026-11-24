import os
print(sorted(f for f in os.listdir('.') if f.endswith('.py')))
print('---')
print(open('step_09.py').read()[:400])
