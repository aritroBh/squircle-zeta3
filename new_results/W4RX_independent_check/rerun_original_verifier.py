"""W4RX: RERUN of the original standalone verifier A3h5_exact_certificate/verify_D.py, with its
log redirected to stdout (its own log file is not written).  This is a rerun, not an
independent check."""
import sys, time
sys.dont_write_bytecode = True
sys.path.insert(0, '<project>/research/08_theory/A3h5_exact_certificate')
import verify_D as V
V.log = lambda s: print(time.strftime('%Y-%m-%d %H:%M:%S ') + '[original verify_D] ' + s, flush=True)
t = time.time(); V.main(); print(f'elapsed {time.time()-t:.1f}s')
