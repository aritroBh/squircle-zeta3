import os
base = "<project>/research/08_theory"
for d in [".", base, os.path.join(base, "Q1S_vertex_symbol"), os.path.join(base, "Q1S2_vertex_close")]:
    print("=== ", d)
    try:
        for f in sorted(os.listdir(d)):
            p = os.path.join(d, f)
            print(f, os.path.getsize(p) if os.path.isfile(p) else "<dir>")
    except Exception as e:
        print("ERR", e)
