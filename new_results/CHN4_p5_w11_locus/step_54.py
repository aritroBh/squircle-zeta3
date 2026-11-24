import os, subprocess
print("out:", open('run_proof_div.out').read() if os.path.exists('run_proof_div.out') else "none")
pid=open('run_proof_div.pid').read().strip()
r=subprocess.run(['ps','-p',pid,'-o','pid,etime,%cpu,command'],capture_output=True,text=True)
print(r.stdout, r.stderr)
