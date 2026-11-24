import os
base = '<project>'
for sub in ['research/08_theory/M4a_decomposition', 'research/08_theory/M4c_w5_closed_form', 
            'research/08_theory/A3_W4_proof', 'paper/squircle_note/scripts']:
    p = os.path.join(base, sub)
    print('='*80)
    print(p)
    for root, dirs, files in os.walk(p):
        depth = root.replace(p, '').count(os.sep)
        if depth <= 2:
            indent = '  '*depth
            print(f"{indent}{os.path.basename(root)}/")
            for f in sorted(files):
                if not f.endswith('.pyc'):
                    fp = os.path.join(root, f)
                    sz = os.path.getsize(fp)
                    print(f"{indent}  {f} ({sz} bytes)")
