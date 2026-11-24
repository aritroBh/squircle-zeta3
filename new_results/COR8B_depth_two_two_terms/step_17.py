
import os
for p in ['../COR8_depth_two_any_angle/result.md']:
    print("="*100); print(p)
    print(open(p).read()[:8000])
