import subprocess
PY = '<project>/.venv/bin/python'
test = '''
import mpmath as mp, time
mp.mp.dps = 400
t=time.time()
Nmax = 2**14
win = []
acc = mp.mpf(0)
for n in range(1, Nmax+1):
    acc += 1/mp.mpf(n)**2
    if n > Nmax-400: win.append(acc)
print("sums done %.1fs" % (time.time()-t))
val, maxc = mp.richardson(win)
mp.mp.dps=30
print("richardson err", mp.nstr(abs(val-mp.pi**2/6),3), "maxc", mp.nstr(maxc,3))
# compare: raw last partial error
mp.mp.dps=400
print("raw partial err", mp.nstr(abs(win[-1]-mp.pi**2/6),3))
# and a smaller window at same precision
val2, maxc2 = mp.richardson(win[-40:])
mp.mp.dps=30
print("M=40 err", mp.nstr(abs(val2-mp.pi**2/6),3), "maxc", mp.nstr(maxc2,3))
'''
r = subprocess.run([PY,'-c',test], capture_output=True, text=True, timeout=600)
print(r.stdout); print(r.stderr[-500:])
