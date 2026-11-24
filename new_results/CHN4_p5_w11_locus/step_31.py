import time, os
time.sleep(120)
lines=open('f11_evals.log').read().strip().split('\n')
print(f"{len(lines)} points done.")
for l in lines[-2:]: print(" ", l)
