import os
base = os.path.dirname(os.getcwd())
with open(os.path.join(base,'DOD2_dodecagon_finish','result.md')) as f:
    print(f.read())
