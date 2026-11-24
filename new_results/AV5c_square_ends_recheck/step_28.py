import os
print(sorted(os.listdir('.')))
if os.path.exists('progress.log'):
    print(open('progress.log').read())
