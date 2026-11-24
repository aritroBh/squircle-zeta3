import subprocess, os, time
r = subprocess.run(['ps', '-p', open('/tmp/g1_k8.pid').read().strip(), '-o', 'pid,etime,command'], capture_output=True, text=True)
print(r.stdout, r.stderr)
print(open('/tmp/g1_k8.out').read()[-1500:] if os.path.exists('/tmp/g1_k8.out') else 'no out')
