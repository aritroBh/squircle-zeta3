import os, time
pid = int(open('/tmp/g1_k8.pid').read().strip())
try:
    os.kill(pid, 0); alive = True
except OSError:
    alive = False
out = open('/tmp/g1_k8.out').read()
print('alive:', alive, '| output bytes:', len(out))
print(out[-2500:])
