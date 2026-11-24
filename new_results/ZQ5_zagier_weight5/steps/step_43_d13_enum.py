
import sys; sys.path.insert(0,'.')
import json
import w5lib as W
import pipe
C, cands = pipe.phase_enum(13, B=8)
pipe.print_cands(13, cands, C)
K = pipe.phase_kernel(13, C, cands)
json.dump(dict(K=K), open(pipe.cand_dir(13)+"/kernel.json","w"))
print("kernel dim", len(K))
