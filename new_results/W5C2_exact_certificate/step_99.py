import os, signal, subprocess, time
os.chdir('<project>/research/08_theory/W5C2_exact_certificate')
try:
    os.kill(int(open('step_03_test.pid').read()), signal.SIGKILL); print('killed test')
except Exception as e: print(e)
# patch dixon to free Amod promptly
s = open('w5c2_dixon.py').read()
old = """    t0 = time.monotonic()
    Amod = flint.nmod_mat((A % p).astype(object).tolist(), p)
    log(f'dixon: nmod A built {time.monotonic()-t0:.1f}s; inverting')
    t0 = time.monotonic()
    Ainv = np.array(Amod.inv().tolist(), np.int64)
    log(f'dixon: nmod inv done {time.monotonic()-t0:.1f}s')
    A_nmod = flint.nmod_mat((A % p).astype(object).tolist(), p)  # reuse? not needed
    del A_nmod"""
new = """    t0 = time.monotonic()
    Amod = flint.nmod_mat((A % p).astype(object).tolist(), p)
    log(f'dixon: nmod A built {time.monotonic()-t0:.1f}s; inverting')
    t0 = time.monotonic()
    Ainv_m = Amod.inv()
    del Amod
    Ainv = np.array(Ainv_m.tolist(), np.int64)
    del Ainv_m
    import gc; gc.collect()
    log(f'dixon: nmod inv done {time.monotonic()-t0:.1f}s')"""
assert old in s
open('w5c2_dixon.py','w').write(s.replace(old, new))
print('dixon patched')
# check step_03.py content once more
print(open('step_03.py').read()[:600])