import os
cwd = os.getcwd()
# Start progress.log and result.md scaffolding
if not os.path.exists('progress.log'):
    open('progress.log','w').write("")
open('progress.log','a').write("step_01: explored COR8B/COR8R/paper; plan fixed. Writing own engine next.\n")
print(os.listdir(cwd))
