
import sys; sys.path.insert(0,'.')
import json
import w5lib as W
import pipe
C = W.letters_info(12)
J, cands = pipe.load_cands(12)
K = pipe.phase_kernel(12, C, cands)
json.dump(dict(K=K), open(pipe.cand_dir(12)+"/kernel.json","w"))
print("kernel saved, dim", len(K))
