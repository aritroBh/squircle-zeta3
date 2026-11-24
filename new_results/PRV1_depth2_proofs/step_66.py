import os
wd = 'research/08_theory/PRV1_depth2_proofs'
steps = sorted(f for f in os.listdir(wd) if f.startswith('step_'))
print("steps:", steps[0], "...", steps[-1], " count:", len(steps))

with open(f'{wd}/result.md') as f:
    txt = f.read()
old = "`step_01.py`–`step_NN.py` (derivation trail"
new = ("`step_01.py`–`step_62.py` (derivation trail")
assert old in txt
txt = txt.replace(old, new)
with open(f'{wd}/result.md','w') as f:
    f.write(txt)
print("file list patched")
# directory listing for the record
for f in sorted(os.listdir(wd)):
    if not f.startswith('step_'):
        print("  ", f, os.path.getsize(os.path.join(wd,f)))
