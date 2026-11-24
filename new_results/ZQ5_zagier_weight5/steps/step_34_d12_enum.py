
import sys; sys.path.insert(0,'.')
import json
import w5lib as W
import pipe
C, cands = pipe.phase_enum(12, B=8)
pipe.print_cands(12, cands, C)
K = pipe.phase_kernel(12, C, cands)
json.dump(dict(K=K), open(pipe.cand_dir(12)+"/kernel.json","w"))
