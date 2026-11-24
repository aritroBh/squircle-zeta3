import os, signal
for pidf in ['/tmp/g1_k5.pid', '/tmp/g1_k8.pid']:
    try:
        pid = int(open(pidf).read().strip())
        os.kill(pid, signal.SIGKILL)
        print('killed', pid)
    except Exception as e:
        print(pidf, e)
open('progress.log','a').write('step_53: unmodified G1 at K=8 died (OOM, 0 output) and K=5 exceeded ~25 min; killed. Bonus check stands at K=4 (exact match). Required k=7,8 recomputation done with my own code.\n')
print('done')
