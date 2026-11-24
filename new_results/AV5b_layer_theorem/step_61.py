import os
work = "<project>/research/08_theory/AV5b_layer_theorem"
os.chdir(work)
print(open("progress.log").read())
print("=== summary.md exists:", os.path.exists("summary.md"))
print(sorted(os.listdir(".")))
