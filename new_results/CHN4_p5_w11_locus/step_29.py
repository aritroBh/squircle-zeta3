import time, os
time.sleep(45)
if os.path.exists('f11_evals.log'):
    lines=open('f11_evals.log').read().strip().split('\n')
    print(f"{len(lines)} points done. Last few:")
    for l in lines[-3:]: print(" ", l)
else:
    print("log not yet created")
print("out file:", open('run_f11_evals.out').read() if os.path.exists('run_f11_evals.out') else 'none')
