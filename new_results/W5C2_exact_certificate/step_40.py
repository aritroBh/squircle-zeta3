import subprocess, os, time
os.chdir('<project>/research/08_theory/W5C2_exact_certificate')
st1 = os.stat('step_01.out')
print('mtime', time.ctime(st1.st_mtime), 'size', st1.st_size, 'now', time.ctime())
time.sleep(45)
st2 = os.stat('step_01.out')
print('after 45s: mtime', time.ctime(st2.st_mtime), 'size', st2.st_size)
print(subprocess.run(['bash','-c','ps -A -o pid,command 2>/dev/null | head'], capture_output=True, text=True).stdout)