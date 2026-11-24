import os
files = sorted(os.listdir('.'), key=lambda s: (len(s), s))
print(files)
