import multiprocessing as mpz, time, os
def f(x):
    return x*x
t0 = time.time()
try:
    with mpz.Pool(4) as pool:
        out = pool.map(f, range(10))
    print("pool OK in %.1fs:" % (time.time()-t0), out)
except Exception as e:
    print("FAILED:", e)
print("start method:", mpz.get_start_method())
print("cpu count:", os.cpu_count())
