import os
print(os.listdir('.'))
print("---- progress.log ----")
with open('progress.log') as f:
    print(f.read())
