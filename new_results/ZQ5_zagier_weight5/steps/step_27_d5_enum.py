
import sys; sys.path.insert(0,'.')
import time
import w5lib as W
import pipe
C, cands = pipe.phase_enum(5, B=8)
pipe.print_cands(5, cands, C)
K = pipe.phase_kernel(5, C, cands)
import json
json.dump(dict(K=K), open(pipe.cand_dir(5)+"/kernel.json","w"))
