import subprocess
PY = '<project>/.venv/bin/python'
test = '''
import mpmath as mp
mp.mp.dps = 50
def neville(a, verbose=False):
    R = list(a); M = len(R); tabs=[R]
    for k in range(1, M):
        new = []
        for i in range(M-k):
            num = mp.mpf(i+1+k)**k * R[i+1] - mp.mpf(i+1)**k * R[i]
            den = mp.mpf(i+1+k)**k - mp.mpf(i+1)**k
            new.append(num/den)
        R = new; tabs.append(R)
    return tabs
# synthetic: a_n = 1 + 1/n, exact limit 1
a = [1+1/mp.mpf(n) for n in range(1, 9)]
tabs = neville(a)
for k, t in enumerate(tabs): print(k, [mp.nstr(x, 8) for x in t[:4]])
'''
r = subprocess.run([PY,'-c',test], capture_output=True, text=True, timeout=300)
print(r.stdout)
