import sys, time, os
os.chdir('<project>/research/08_theory/DEC5E_decagon_p5_closed')
sys.path.insert(0,'.')
import gengine_dec2p, mpmath as mp
# time one shard at N=380, dps=120, in-process
t0=time.time()
LETTERS, INVP, FRAT, POWH, POWL = gengine_dec2p._setup(120, 380)
print("setup %.2fs" % (time.time()-t0))
t0=time.time()
res = gengine_dec2p._shard(((LETTERS[0], LETTERS[2], LETTERS[3]), 120, 380))
print("one shard %.2fs, %d entries" % (time.time()-t0, len(res)))
