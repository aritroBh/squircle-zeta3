import time, os
# give the background runs some time; check every 30s up to ~4 min
for i in range(8):
    time.sleep(30)
    def st(pidf, outf):
        pid = int(open(pidf).read().strip())
        try:
            os.kill(pid, 0); alive = True
        except OSError:
            alive = False
        return alive, os.path.getsize(outf)
    a5, b5 = st('/tmp/g1_k5.pid', '/tmp/g1_k5.out')
    a8, b8 = st('/tmp/g1_k8.pid', '/tmp/g1_k8.out')
    print(f'[{(i+1)*30}s] K=5: alive={a5} bytes={b5} | K=8: alive={a8} bytes={b8}', flush=True)
    if not a5 and b5 > 0 and not a8 and b8 > 0:
        break
