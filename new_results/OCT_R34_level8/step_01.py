import sys
import os

# Check the working directory
print("Current working directory:", os.getcwd())

# List relevant directories
base_path = 'research/08_theory'
print("\nChecking directories:")
for name in ['SQ_chi3_bridge', 'A3f4_glanois_basis', 'A3f6_level16', 'A3f3_level8_compute', 'OCT6_octagon_closed', 'CX1_octagon_depth2_proof']:
    path = os.path.join(base_path, name)
    if os.path.exists(path):
        print(f"  {name}: EXISTS")
    else:
        print(f"  {name}: NOT FOUND")