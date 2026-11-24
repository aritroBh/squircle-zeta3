
import sys; sys.path.insert(0,'.')
import json, time
import w5lib as W
import pipe
D=17
C = W.letters_info(D)
J, cands = pipe.load_cands(D)
t0=time.time()
K = pipe.phase_kernel(D, C, cands)
json.dump(dict(K=K), open(pipe.cand_dir(D)+"/kernel.json","w"))
print(f"kernel {time.time()-t0:.0f}s dim {len(K)}")
