import time, os, json
def st(pidf, outf):
    pid = int(open(pidf).read().strip())
    try:
        os.kill(pid, 0); alive = True
    except OSError:
        alive = False
    return pid, alive, os.path.getsize(outf)
# wait up to 6 more minutes for K=5 (its series step should take ~5-10 min at K=5)
for i in range(12):
    time.sleep(30)
    _, a5, b5 = st('/tmp/g1_k5.pid', '/tmp/g1_k5.out')
    _, a8, b8 = st('/tmp/g1_k8.pid', '/tmp/g1_k8.out')
    print(f'[+{(i+1)*30}s] K=5 alive={a5} bytes={b5} | K=8 alive={a8} bytes={b8}', flush=True)
    if b5 > 0 and b8 > 0:
        break
    if b5 > 0 and not a5:
        break
