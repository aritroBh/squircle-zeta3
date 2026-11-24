# step_37: symbol-level proof s=3 (n=9) and s=4 (n=14) -- exact flint ranks
exec(open("step_36.py").read().split("run_proof(1, verbose=True)")[0])
import time
t0 = time.time()
run_proof(3)
print("s=3 time:", round(time.time()-t0,1), "s")
