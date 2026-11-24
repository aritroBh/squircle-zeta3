import pickle
from pathlib import Path
from fractions import Fraction
import numpy as np

base_path = Path('research/08_theory/')

# Load weight3_with_lifts data
print("Loading weight3_with_lifts.pickle...")
with open(base_path / 'A3f3_level8_compute/weight3_with_lifts.pickle', 'rb') as f:
    w3_data = pickle.load(f)

print(f"Columns (words): {len(w3_data['cols'])}")
print(f"Rows: {len(w3_data['rows'])}")
print(f"Rank over Q (from data): {w3_data.get('rank_Q', 'Not found')}")
print(f"Quotient: {w3_data.get('quotient', 'Not found')}")

# Examine a few rows
print("\nFirst few rows:")
for i in range(min(3, len(w3_data['rows']))):
    row = w3_data['rows'][i]
    if isinstance(row, list):
        print(f"Row {i}: {len(row)} elements")
        if row:
            print(f"  First element: {row[0]}")
    elif isinstance(row, tuple):
        print(f"Row {i}: tuple of length {len(row)}")
        if len(row) > 0:
            print(f"  First element type: {type(row[0])}")

# Load k3_newHrows data
print("\n\nLoading k3_newHrows.pickle...")
with open(base_path / 'A3f3_level8_compute/k3_newHrows.pickle', 'rb') as f:
    k3_data = pickle.load(f)

print(f"Number of Kummer rows: {len(k3_data['krows'])}")
if k3_data['krows']:
    print(f"First Kummer row keys: {list(k3_data['krows'][0].keys())[:5]}")
    first_key = list(k3_data['krows'][0].keys())[0]
    print(f"First Kummer row value for key {first_key}: {k3_data['krows'][0][first_key]}")

# Load five_relations data
print("\n\nLoading five_relations.pickle...")
with open(base_path / 'A3f4_glanois_basis/five_relations.pickle', 'rb') as f:
    rel_data = pickle.load(f)

print(f"Free words: {len(rel_data['freewords'])}")
print(f"Relations: {len(rel_data['relations'])}")
print(f"First relation type: {type(rel_data['relations'][0])}")
print(f"First relation length: {len(rel_data['relations'][0])}")
print(f"First few values of first relation: {rel_data['relations'][0][:10]}")
print(f"Free words first few: {rel_data['freewords'][:10]}")
