import pickle, os, json, sys
from pathlib import Path

base_path = Path('research/08_theory/')
print("Base path exists:", base_path.exists())
print("\nAvailable files:")

files_to_check = [
    'A3f3_level8_compute/weight3_with_lifts.pickle',
    'A3f3_level8_compute/k3_newHrows.pickle', 
    'A3f4_glanois_basis/five_relations.pickle'
]

for f in files_to_check:
    full_path = base_path / f
    print(f"{f}: {full_path.exists()}")
    if full_path.exists():
        try:
            with open(full_path, 'rb') as fp:
                data = pickle.load(fp)
                if isinstance(data, dict):
                    print(f"  Keys: {list(data.keys())}")
                    # Print some basic info about the data
                    for key in data:
                        if isinstance(data[key], list):
                            print(f"  {key}: list of length {len(data[key])}")
                            if len(data[key]) > 0 and isinstance(data[key][0], (list, tuple, dict)):
                                print(f"    First elem type: {type(data[key][0])}")
                        elif isinstance(data[key], (int, float)):
                            print(f"  {key}: {data[key]}")
                        else:
                            print(f"  {key}: type {type(data[key])}")
                elif isinstance(data, list):
                    print(f"  List of length {len(data)}")
                    if len(data) > 0:
                        print(f"  First elem type: {type(data[0])}")
                else:
                    print(f"  Type: {type(data)}")
        except Exception as e:
            print(f"  Error loading: {e}")

print("\nChecking for other directories:")
for item in base_path.iterdir():
    if item.is_dir():
        print(f"  {item.name}")
