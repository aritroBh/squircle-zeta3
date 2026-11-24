import os
base = '<project>/research/08_theory'
for d in ['ZQ5b_zagier_w5_larger_S', 'ZQ5_zagier_weight5']:
    p = os.path.join(base, d)
    print("="*70)
    print(d)
    for root, dirs, files in os.walk(p):
        depth = root[len(p):].count(os.sep)
        if depth > 2: continue
        indent = '  '*depth
        print(f"{indent}{os.path.basename(root)}/")
        for f in sorted(files):
            if f.endswith(('.py','.json','.log','.md','.txt')):
                sz = os.path.getsize(os.path.join(root,f))
                print(f"{indent}  {f}  ({sz} bytes)")
