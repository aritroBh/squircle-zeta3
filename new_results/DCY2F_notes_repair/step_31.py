import os, shutil, time
SRC = "<project>/research/08_theory/DCY2_decay_first_correction"
CWD = "<project>/research/08_theory/DCY2F_notes_repair"
os.makedirs(CWD, exist_ok=True)

# progress.log
with open(os.path.join(CWD, "progress.log"), "w") as f:
    f.write("DCY2F notes-repair log (apply referee R1-R6 to DCY2, re-certify)\n")

def logline(msg):
    with open(os.path.join(CWD, "progress.log"), "a") as f:
        f.write(msg.rstrip("\n") + "\n")
    print("LOG:", msg)

def read(fn):
    with open(os.path.join(SRC, fn)) as f: return f.read()
def write(fn, txt):
    with open(os.path.join(CWD, fn), "w") as f: f.write(txt)

# 1) verbatim copies (already repaired where applicable)
for fn in ["dcy2_sweep.py","dcy2_regionc.py","dcy2_cert.py","dcy2_common.py",
           "dcy2_audit.py","step_02_varpi1.py","step_03_watson.py","step_07_partb.py"]:
    write(fn, read(fn))
logline("copied verbatim: dcy2_sweep,dcy2_regionc,dcy2_cert,dcy2_common,dcy2_audit,step_02_varpi1,step_03_watson,step_07_partb")
print("done")
