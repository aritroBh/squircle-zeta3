import subprocess
PY = '<project>/.venv/bin/python'
diag = r'''
import mpmath as mp
mp.mp.dps = 60
# plain zeta(2) partial sums at 2^j
for K in [10, 13]:
    partials = []
    acc = mp.mpf(0)
    targets = [2**j for j in range(3, K+1)]
    ti = 0
    for n in range(1, 2**K+1):
        acc += 1/mp.mpf(n)**2
        if ti < len(targets) and n == targets[ti]:
            partials.append(acc); ti += 1
    val, err = mp.richardson(partials)
    print("K=%d  zeta2 err %s  est %s" % (K, mp.nstr(abs(val-mp.pi**2/6),3), mp.nstr(err,3)))
# what about alternating zeta(1,-1) with Euler average
partials = []
acc = mp.mpf(0)
K=13
targets = [2**j for j in range(3, K+1)]
ti = 0
for n in range(1, 2**K+1):
    acc += (-1)**n/mp.mpf(n)
    if ti < len(targets) and n == targets[ti]:
        acc2 = acc + (-1)**(n+1)/(2*mp.mpf(n+1))
        partials.append(acc2); ti += 1
val, err = mp.richardson(partials)
print("log2 err %s est %s" % (mp.nstr(abs(val+mp.log(2)),3), mp.nstr(err,3)))
# shanks instead
partials = []
acc = mp.mpf(0)
ti = 0
for n in range(1, 2**K+1):
    acc += 1/mp.mpf(n)**2
    if ti < len(targets) and n == targets[ti]:
        partials.append(acc); ti += 1
val, err = mp.shanks(partials)
print("shanks zeta2 err %s" % mp.nstr(abs(val-mp.pi**2/6),3))
'''
r = subprocess.run([PY,'-c',diag], capture_output=True, text=True)
print(r.stdout); print(r.stderr[-500:])
