import os, subprocess, textwrap

base = "<project>"
work = os.path.join(base, "research/08_theory/AV5b_layer_theorem")
os.makedirs(work, exist_ok=True)
os.chdir(work)
print("CWD:", os.getcwd())

# Look at AV5 folder structure
av5 = os.path.join(base, "research/08_theory/AV5_layer_universality")
for root, dirs, files in os.walk(av5):
    depth = root.replace(av5, "").count(os.sep)
    if depth <= 1:
        print(root.replace(av5, "AV5"), "->", sorted(files)[:20])
