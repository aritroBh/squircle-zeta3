import os
# find paper/ieee/main.tex and E5_edge_law
root = "<project>/"
for cand in ["paper/ieee/main.tex"]:
    p = os.path.join(root, cand)
    print(p, os.path.exists(p))
# find E5_edge_law
for dirpath, dirnames, filenames in os.walk(os.path.join(root, "research")):
    if "E5_edge_law" in dirpath:
        print(dirpath, len(filenames), "files")
        break
for dirpath, dirnames, filenames in os.walk(root):
    if dirpath.endswith("E5_edge_law"):
        print("FOUND:", dirpath, len(filenames))
