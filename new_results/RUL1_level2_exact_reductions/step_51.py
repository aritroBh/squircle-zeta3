import subprocess
PY = '<project>/.venv/bin/python'
test = '''
import mpmath as mp
def neville(a):
    R = list(a); M = len(R)
    for k in range(1, M):
        new = []
        for i in range(M-k):
            num = mp.mpf(i+1+k)**k * R[i+1] - mp.mpf(i+1)**k * R[i]
            den = mp.mpf(i+1+k)**k - mp.mpf(i+1)**k
            new.append(num/den)
        R = new
    return R[0]
for dps in [60, 100, 150, 200]:
    mp.mp.dps = dps
    M = 150
    partials = []
    acc = mp.mpf(0)
    for n in range(1, M+1):
        acc += 1/mp.mpf(n)**2
        partials.append(acc)
    v = neville(partials)
    print("dps=%3d M=%3d zeta2 err %s" % (dps, M, mp.nstr(abs(v-mp.pi**2/6),3)))
# try smaller M
for M in [40, 60, 80, 100]:
    mp.mp.dps = 120
    partials = []
    acc = mp.mpf(0)
    for n in range(1, M+1):
        acc += 1/mp.mpf(n)**2
        partials.append(acc)
    v = neville(partials)
    print("dps=120 M=%3d zeta2 err %s" % (M, mp.nstr(abs(v-mp.pi**2/6),3)))
'''
r = subprocess.run([PY,'-c',test], capture_output=True, text=True, timeout=600)
print(r.stdout); print(r.stderr[-300:])
