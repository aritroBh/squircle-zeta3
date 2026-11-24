import pickle
from pathlib import Path
from fractions import Fraction

base_path = Path('research/08_theory/')

# Load data
with open(base_path / 'A3f3_level8_compute/weight3_with_lifts.pickle', 'rb') as f:
    w3_data = pickle.load(f)

print("Examining data structure...")
print(f"Type of rows: {type(w3_data['rows'])}")
print(f"Length of rows: {len(w3_data['rows'])}")

# Examine first few rows carefully
for i in range(3):
    row = w3_data['rows'][i]
    print(f"\nRow {i}:")
    print(f"  Type: {type(row)}")
    print(f"  Length: {len(row)}")
    print(f"  Content: {row}")
    
    # Try to understand the tuple structure
    if isinstance(row, tuple):
        print(f"  Elements:")
        for j, elem in enumerate(row):
            print(f"    {j}: type={type(elem)}, content={elem}")
            if isinstance(elem, tuple) and len(elem) == 2:
                print(f"      Is tuple of length 2: {elem}")
                print(f"      Element 0 type: {type(elem[0])}, value: {elem[0]}")
                print(f"      Element 1 type: {type(elem[1])}, value: {elem[1]}")

# Let's also check what the expected rank is from the data
print(f"\nData says rank_Q: {w3_data.get('rank_Q')}")
print(f"Data says quotient: {w3_data.get('quotient')}")

# Check the actual format by trying to process a row manually
row0 = w3_data['rows'][0]
print(f"\nProcessing row0 manually:")
if isinstance(row0, tuple):
    for elem in row0:
        if isinstance(elem, tuple) and len(elem) == 2:
            word, coeff = elem
            print(f"  Word: {word}, Coeff: {coeff}, Coeff type: {type(coeff)}")
