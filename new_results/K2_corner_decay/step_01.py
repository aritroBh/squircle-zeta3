
import os
# Explore current directory and referenced directories
for d in ['.', '../D10c_corner_general', '../D13b_Mk_constant']:
    print("="*60)
    print("DIR:", os.path.abspath(d))
    try:
        for f in sorted(os.listdir(d)):
            print("   ", f)
    except Exception as e:
        print("   ERR", e)
