import pickle
from pathlib import Path

base_path = Path('research/08_theory/')

# Check if level-16 data exists
l16_path = base_path / 'A3f6_level16/l16_w3.pickle'
print(f"Level-16 data path exists: {l16_path.exists()}")

if l16_path.exists():
    try:
        with open(l16_path, 'rb') as f:
            l16_data = pickle.load(f)
        
        print("\nLevel-16 data loaded successfully")
        print(f"Keys: {list(l16_data.keys())}")
        
        if 'cols' in l16_data:
            print(f"Number of columns (words): {len(l16_data['cols'])}")
        
        if 'rows' in l16_data:
            print(f"Number of rows: {len(l16_data['rows'])}")
            if len(l16_data['rows']) > 0:
                print(f"First row type: {type(l16_data['rows'][0])}")
        
        if 'counts' in l16_data:
            print(f"Counts: {l16_data['counts']}")
            
    except Exception as e:
        print(f"Error loading level-16 data: {e}")
        
# Also check for result.md in the same directory
result_md_path = base_path / 'A3f6_level16/result.md'
print(f"\nresult.md exists: {result_md_path.exists()}")

# Check for l16.py
l16_py_path = base_path / 'A3f6_level16/l16.py'
print(f"l16.py exists: {l16_py_path.exists()}")

if l16_py_path.exists():
    try:
        with open(l16_py_path, 'r') as f:
            content = f.read(500)  # Read first 500 chars
            print(f"\nFirst 500 chars of l16.py:\n{content}")
    except Exception as e:
        print(f"Error reading l16.py: {e}")