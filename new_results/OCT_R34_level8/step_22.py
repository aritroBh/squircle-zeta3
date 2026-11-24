# Check pickle file content
import pickle
V = pickle.load(open('research/08_theory/A3f6_level16/vals8_150digit.pickle', 'rb'))
print(f"Type of V: {type(V)}")
print(f"Number of entries: {len(V)}")
print(f"Sample keys: {list(V.keys())[:10]}")
# Check format
sample_key = list(V.keys())[0]
print(f"Sample entry: key={sample_key}, value type={type(V[sample_key])}, value={V[sample_key]}")