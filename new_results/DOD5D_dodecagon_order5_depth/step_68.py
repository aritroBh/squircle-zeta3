import subprocess, pickle, os
src = '''import pickle, sys, time
sys.path.insert(0, ".")
import gengine6
slotwords = pickle.load(open("needed_words.pkl","rb"))
t0 = time.time()
print("building dps=130 N=320 for %d words" % len(slotwords), flush=True)
G1, LET = gengine6.gtab(slotwords, dps=130, N=320, cachedir=".")
print("built %d values in %.1fs" % (len(G1), time.time()-t0), flush=True)
# self-checks
import mpmath as mp
mp.mp.dps = 130
r3 = mp.sqrt(3)
p0 = 4+2*r3; p2 = -(1+r3)/2; p3 = (3+r3)/2
c1 = abs(G1[(0,0,0,0,2)] + mp.polylog(5, 1/p0))
print("G(0^4,p0;1) = -Li_5(1/p0): resid", mp.nstr(c1, 3), flush=True)
g_pq = G1[(2,4)]
ref = mp.quad(lambda t: mp.log(1-t/( -(1+r3)/2))/(t-(4+2*r3)), [0,1])
print("G(p0,p2;1) vs quad: resid", mp.nstr(abs(g_pq-ref),3), flush=True)
# shuffle check: G(p0)G(p2) = G(p0,p2)+G(p2,p0)
sh = G1[(2,)]*G1[(4,)] - G1[(2,4)] - G1[(4,2)]
print("shuffle: resid", mp.nstr(abs(sh),3), flush=True)
# G(0^4,3;1) = -Li_5(1/p1)
c2 = abs(G1[(0,0,0,0,3)] + mp.polylog(5, 1/(-(3+2*r3))))
print("G(0^4,p1;1) = -Li_5(1/p1): resid", mp.nstr(c2,3), flush=True)
# log: G(p3;1) = log(1-1/p3) = log(1+kappa_3)
c3 = abs(G1[(5,)] - mp.log(1+ (r3-3)/3 ))
print("G(p3;1) = log(1+k3): resid", mp.nstr(c3,3), flush=True)
c4 = abs(G1[(0,)] + 0)  # G(0;1) regularized = 0 (log 1)
print("G(0;1):", mp.nstr(G1[(0,)],3), flush=True)
'''
with open("gtab_build.py","w") as f: f.write(src)
logf = open("gtab_build.log","w")
proc = subprocess.Popen(["<project>/.venv/bin/python", "-u", "gtab_build.py"],
                        stdout=logf, stderr=logf, cwd=os.getcwd())
pickle.dump(proc.pid, open("gtab_pid.pkl","wb"))
print("launched", proc.pid)
print("--- itab progress ---")
print(open("itab_build.log").read())
