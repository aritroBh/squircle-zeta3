import os
print(os.getcwd())
print(os.listdir("."))
if os.path.exists("step_01.py"):
    print(open("step_01.py").read()[:2000])
