import os, time
print(open('verify_gcd.log').read())
pid=int(open('verify_gcd_pid.txt').read().strip())
print('pid alive:', os.path.exists(f'/proc/{pid}') if os.path.exists('/proc') else 'unknown (macOS)')
import subprocess
r=subprocess.run(['ps','-p',str(pid),'-o','pid,etime,command'],capture_output=True,text=True)
print(r.stdout, r.stderr)
