src = '''import pickle, sys, time
sys.path.insert(0, ".")
from qs3_symbol import ZERO, ONE, g_symbol
import dod5d_grewrite as GR
rewrite = pickle.load(open("rewrite_exact.pkl","rb"))
words = set()
for r in rewrite.values(): words.update(r.keys())
words = sorted(words, key=lambda w: (len(w), w))
t0 = time.time()
tab = {}
for i, w in enumerate(words):
    tab[w] = g_symbol(ZERO, w, ONE)
    if i % 1000 == 0: print(i, "/", len(words), "%.0fs" % (time.time()-t0), flush=True)
pickle.dump(tab, open("gsym_tab.pkl","wb"))
print("saved gsym_tab.pkl:", len(tab), "words, %.0fs" % (time.time()-t0))
'''
with open("gsym_build.py","w") as f: f.write(src)
import subprocess, pickle
logf = open("gsym_build.log","w")
proc = subprocess.Popen(["<project>/.venv/bin/python", "-u", "gsym_build.py"],
                        stdout=logf, stderr=logf)
pickle.dump(proc.pid, open("gsym_pid.pkl","wb"))
print("launched gsym build, pid", proc.pid)
