import pickle
from pathlib import Path
from fractions import Fraction
import datetime

base_path = Path('research/08_theory/')

# First implement the embedding function
def embed_level8_to_level16(word):
    """
    Embed level-8 word (7, b, c) to level-16 word.
    According to j ↦ 2j mapping:
    - Letter -1 stays -1
    - Letter 7 becomes 14 (2*7)
    - Letter b ∈ {4,5,6,7} becomes 2b ∈ {8,10,12,14}
    - Letter c ∈ {0..7} becomes 2c ∈ {0,2,4,...,14}
    So (7, b, c) ↦ (14, 2b, 2c)
    """
    if isinstance(word, tuple) and len(word) == 3:
        a, b, c = word
        # a should be 7 for level-8 freewords
        return (2*a, 2*b, 2*c)
    else:
        raise ValueError(f"Unexpected word format: {word}")

# Test the embedding
test_word = (7, 4, 0)
embedded = embed_level8_to_level16(test_word)
print(f"Test embedding: {test_word} → {embedded}")

# Load all necessary data
print("\nLoading data...")
with open(base_path / 'A3f6_level16/l16_w3.pickle', 'rb') as f:
    l16_data = pickle.load(f)

with open(base_path / 'A3f4_glanois_basis/five_relations.pickle', 'rb') as f:
    rel_data = pickle.load(f)

# Check embedding for all freewords
freewords = rel_data['freewords']
print(f"\nChecking embedding for freewords:")
for i, word in enumerate(freewords[:5]):
    embedded_word = embed_level8_to_level16(word)
    print(f"  {word} → {embedded_word}")

# Now check if these embedded words exist in level-16 vocabulary
l16_vocab = set(l16_data['cols'])
print(f"\nChecking if embedded words are in level-16 vocabulary:")
all_found = True
for word in freewords[:10]:
    embedded_word = embed_level8_to_level16(word)
    if embedded_word not in l16_vocab:
        print(f"  WARNING: {word} → {embedded_word} NOT FOUND in level-16 vocabulary")
        all_found = False
    else:
        print(f"  ✓ {word} → {embedded_word} found")

if all_found:
    print("\nAll test words embed successfully into level-16 vocabulary")
else:
    print("\nSome words do not embed properly")

# Let's also check what level-16 words look like
print(f"\nSample level-16 words (first 10):")
for word in l16_data['cols'][:10]:
    print(f"  {word}")

print(f"\nSample level-16 words (last 10):")
for word in l16_data['cols'][-10:]:
    print(f"  {word}")