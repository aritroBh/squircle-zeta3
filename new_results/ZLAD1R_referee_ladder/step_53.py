import os, time, json, subprocess
g1dir = '<project>/research/10_frontier/G1_zeta3_any_angle'
# check for any file I may have caused in the under-review dir
suspects = [f for f in os.listdir(g1dir) if 'K8' in f or 'K5' in f or f.startswith('flat_exact_N40_K8')]
print('suspicious files in G1 dir:', suspects)
for f in suspects:
    print(' ', f, os.path.getsize(os.path.join(g1dir, f)), 'mtime', time.ctime(os.path.getmtime(os.path.join(g1dir, f))))
# K=8 run died: check exit by inspecting if a flat_exact_N40_K8.json exists there
# K=5: wait up to 5 more minutes
for i in range(10):
    time.sleep(30)
    pid = int(open('/tmp/g1_k5.pid').read().strip())
    try:
        os.kill(pid, 0); alive = True
    except OSError:
        alive = False
    b = os.path.getsize('/tmp/g1_k5.out')
    print(f'[+{(i+1)*30}s] K=5 alive={alive} bytes={b}', flush=True)
    if not alive:
        break
