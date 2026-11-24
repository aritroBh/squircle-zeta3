# step_40: symbol-level proof s=4 (n=14 letters)
exec(open("symbolproof_engine.py").read())
import time
t0 = time.time()
run_proof(4)
print("s=4 time:", round(time.time()-t0,1), "s")
