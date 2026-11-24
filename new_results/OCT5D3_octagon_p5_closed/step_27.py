import os, pickle, sys, multiprocessing
print("cpus:", multiprocessing.cpu_count())
try:
    import psutil
    print("RAM GB:", psutil.virtual_memory().total/1e9, "available:", psutil.virtual_memory().available/1e9)
except ImportError:
    print(os.popen("sysctl -n hw.memsize").read())
    print(os.popen("vm_stat | head -5").read())
