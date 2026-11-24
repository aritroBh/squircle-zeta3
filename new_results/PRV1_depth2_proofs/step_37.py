import os
wd = 'research/08_theory/PRV1_depth2_proofs'
# print sizes/mtimes of steps, and contents of last few
files = sorted(f for f in os.listdir(wd) if f.startswith('step_'))
for f in files:
    st = os.stat(os.path.join(wd,f))
    print(f, st.st_size)
