
import os, difflib
base = "<project>"
work = base + "/research/08_theory/TOW4R_review"
a56 = open(work + "/tow3_orig/prove_adj56.py").read()
afast = open(work + "/tow4_copy/prove_adj78_fast.py").read()
print("##### diff prove_adj56 (TOW3) -> prove_adj78_fast (TOW4) #####")
for line in difflib.unified_diff(a56.splitlines(), afast.splitlines(), lineterm='', n=0):
    print(line)
