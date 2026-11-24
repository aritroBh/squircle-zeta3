import os
for f in sorted(os.listdir(".")):
    print(f, os.path.getsize(f))
