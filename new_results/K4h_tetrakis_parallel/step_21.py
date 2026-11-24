import subprocess, os, time
CWD = os.getcwd()
probe = '''
import os, time
open("detach_probe.out","w").write("pid=%d ppid=%d sid=%d nice=%d\\n" % (os.getpid(), os.getppid(), os.getsid(0), os.nice(0)))
time.sleep(2)
'''
open("detach_probe.py","w").write(probe)
VENV = "<project>/.venv/bin/python3"
DETACH = "<project>/tools/detach.py"
out = open("detach_probe_stdout.txt","w")
cmd = ["nohup","nice","-n","10",VENV,DETACH,VENV,os.path.join(CWD,"detach_probe.py")]
p = subprocess.Popen(cmd, stdout=out, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL, cwd=CWD)
print("probe launcher pid:", p.pid)
time.sleep(4)
print("probe out:", open("detach_probe.out").read())
print("launcher alive?", subprocess.run(["ps","-p",str(p.pid)],capture_output=True,text=True).stdout)
print("my sid:", os.getsid(0), "my nice:", os.nice(0))
