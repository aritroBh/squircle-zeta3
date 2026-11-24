import os
# Explore working directory and inputs
for root in ['.', '../DEC5D_decagon_order5_depth', '../DEC2_decagon_exact_assembly', '../OCT5D3_octagon_p5_closed']:
    print('='*80)
    print(root)
    try:
        for f in sorted(os.listdir(root)):
            p = os.path.join(root, f)
            if os.path.isfile(p):
                print(f'  {f}  {os.path.getsize(p)} bytes')
    except Exception as e:
        print('  ERROR:', e)
