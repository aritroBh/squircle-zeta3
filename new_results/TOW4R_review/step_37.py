
import os, difflib
base = "<project>"
work = base + "/research/08_theory/TOW4R_review"
e56 = open(work + "/tow3_orig/prove_E3b56.py").read()
e78 = open(work + "/tow4_copy/prove_E3b78.py").read()
a56 = open(work + "/tow3_orig/prove_adj56.py").read()
afast = open(work + "/tow4_copy/prove_adj78_fast.py").read()
print("##### diff prove_E3b56 (TOW3) -> prove_E3b78 (TOW4) #####")
for line in difflib.unified_diff(e56.splitlines(), e78.splitlines(), lineterm='', n=1):
    print(line)
