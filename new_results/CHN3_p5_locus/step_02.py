import os
cwd = os.getcwd()
print("CWD:", cwd)
print("\n=== current dir ===")
for f in sorted(os.listdir(cwd)):
    print("  ", f)

for d in ['../CHN2_chain_planes', '../CHN1_chain_bodies']:
    print("\n=== ", d, " ===")
    try:
        for f in sorted(os.listdir(d)):
            p = os.path.join(d, f)
            sz = os.path.getsize(p) if os.path.isfile(p) else '<dir>'
            print(f"   {f}  {sz}")
    except Exception as e:
        print("ERR", e)
