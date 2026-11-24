import os
print("cwd files:", os.listdir('.'))
if os.path.exists('step_01.py'):
    print(open('step_01.py').read())
if os.path.exists('progress.log'):
    print(open('progress.log').read())
