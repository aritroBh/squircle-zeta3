
import os
for d in ['../COR4_M3alpha_classical', '../COR8_depth_two_any_angle', '../COR6_true_point_symbols', '../AH3_fn_any_angle']:
    print("="*80); print(d)
    try:
        for f in sorted(os.listdir(d)):
            print("  ", f, os.path.getsize(os.path.join(d,f)))
    except Exception as e:
        print("  ERR", e)
