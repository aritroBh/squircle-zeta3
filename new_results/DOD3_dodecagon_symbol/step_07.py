import os
base = os.path.dirname(os.getcwd())
with open(os.path.join(base,'DOD1_dodecagon','result.md')) as f:
    print(f.read())
