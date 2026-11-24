import os, time, subprocess
pid = int(open('/tmp/g1_k8.pid').read().strip())
# check if process alive via kill(pid, 0)
try:
    os.kill(pid, 0)
    alive = True
except OSError:
    alive = False
print('alive:', alive)
out = open('/tmp/g1_k8.out').read()
print('--- output so far ---')
print(out[-2000:])
