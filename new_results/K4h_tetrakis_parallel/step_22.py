import os, time
print("tool nice:", os.nice(0), "sid:", os.getsid(0))
os.remove("detach_probe.py"); os.remove("detach_probe.out"); os.remove("detach_probe_stdout.txt")
time.sleep(600)
print(open("val96_out.txt").read()[-1500:])
