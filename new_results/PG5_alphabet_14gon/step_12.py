import os, pickle
base = '../DOD5_dodecagon_termwise'
print(sorted(os.listdir(base)))
print(open(os.path.join(base,'result.md')).read()[:6000])
