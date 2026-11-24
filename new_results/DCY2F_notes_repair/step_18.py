import json
base = "<project>/research/08_theory"
with open(f"{base}/DCY2_decay_first_correction/transcript.json") as f:
    tr = f.read()
print("len:", len(tr))
# find occurrences of 1409 and res/bound
import re
for pat in ["1409", "0.00141", "res = abs", "res=abs", "abs(M - T1)", "abs(M + T1)", "abs(M-T1)"]:
    idxs = [m.start() for m in re.finditer(re.escape(pat), tr)]
    print(f"{pat!r}: {len(idxs)} occurrences")
