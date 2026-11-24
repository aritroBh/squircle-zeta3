import os
print(os.getcwd())
for f in sorted(os.listdir('.')):
    print(f)
print("---- step_01.py ----")
with open('step_01.py') as f:
    print(f.read()[:2000])
print("---- step_02.py ----")
with open('step_02.py') as f:
    print(f.read()[:2000])
