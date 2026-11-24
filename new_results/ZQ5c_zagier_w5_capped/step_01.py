import subprocess, os
# First explore the directory structure
base = '<project>'
for d in ['ZQ5c_w5_D29_D13', 'ZQ5b_zagier_w5_larger_S', 'ZQ5_zagier_weight5']:
    p = os.path.join(base, d)
    if os.path.exists(p):
        print("="*70)
        print(d)
        for root, dirs, files in os.walk(p):
            # skip deep
            depth = root[len(p):].count(os.sep)
            if depth > 2: continue
            indent = '  '*depth
            print(f"{indent}{os.path.basename(root)}/")
            for f in files:
                if f.endswith(('.py','.json','.log','.md','.txt')):
                    print(f"{indent}  {f}")
    else:
        print(f"MISSING: {p}")
