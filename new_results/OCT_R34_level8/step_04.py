# Look at the reductions file
import pickle
reductions_path = 'research/08_theory/OCT6_octagon_closed/reductions.pkl'
with open(reductions_path, 'rb') as f:
    reductions = pickle.load(f)
print("Keys:", reductions.keys())
print("\nLam expression:")
print(reductions['Lam'])
print("\nswap expression:")
print(reductions['swap'])