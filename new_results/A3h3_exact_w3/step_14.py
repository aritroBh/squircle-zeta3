import pickle
from pathlib import Path
from fractions import Fraction

base_path = Path('research/08_theory/')

# Load level-16 data to understand word format
with open(base_path / 'A3f6_level16/l16_w3.pickle', 'rb') as f:
    l16_data = pickle.load(f)

print("Level-16 words format:")
print(f"Number of words: {len(l16_data['cols'])}")
print(f"First few words: {l16_data['cols'][:5]}")

# Load level-8 relations to understand their word format
with open(base_path / 'A3f4_glanois_basis/five_relations.pickle', 'rb') as f:
    rel_data = pickle.load(f)

print(f"\nLevel-8 freewords format:")
print(f"Number of freewords: {len(rel_data['freewords'])}")
print(f"First few freewords: {rel_data['freewords'][:5]}")

# According to the problem statement, the embedding is letter j ↦ 2j
# Words are tuples (7, b, c) where b ∈ {4,5,6,7}, c ∈ {0..7}
# In level-16, the letter set should be different

# Let me check if there's pattern in level-16 words
l16_words = l16_data['cols']
print(f"\nAnalyzing level-16 words...")
print(f"First word: {l16_words[0]}")
print(f"Last word: {l16_words[-1]}")

# Check what letters appear in level-16 words
letters_set = set()
for word in l16_words[:100]:  # Sample
    if isinstance(word, tuple):
        letters_set.update(word)

print(f"\nLetters appearing in first 100 level-16 words: {sorted(letters_set)}")

# Now check the mapping
# In level-8: word format is (7, b, c) where b ∈ {4,5,6,7}, c ∈ {0..7}
# The embedding j ↦ 2j means letter j becomes 2j in level-16
# But need to check what letters mean...

# Let me try to understand the letter mapping by reading more of l16.py
l16_py_path = base_path / 'A3f6_level16/l16.py'
with open(l16_py_path, 'r') as f:
    content = f.read(2000)

# Look for mapping information
import re
# Search for patterns about letter mapping
matches = re.findall(r'letter.*?\n|mapping.*?\n|embed.*?\n', content)
if matches:
    print("\nFound in l16.py:")
    for m in matches[:10]:
        print(f"  {m.strip()}")

# Check result.md for embedding information
result_md_path = base_path / 'A3f6_level16/result.md'
if result_md_path.exists():
    with open(result_md_path, 'r') as f:
        result_content = f.read(1000)
    print(f"\nFirst 1000 chars of result.md:")
    print(result_content)