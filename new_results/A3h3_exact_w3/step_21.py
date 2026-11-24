# Final check - list all files created
import os

print("Files in current directory:")
for f in sorted(os.listdir('.')):
    if f.endswith('.py') or f.endswith('.md') or f.endswith('.log'):
        size = os.path.getsize(f)
        print(f"  {f:30} {size:,} bytes")

print("\n" + "="*60)
print("FINAL SUMMARY:")
print("="*60)

# Show key results
print("""
Exact ranks computed over Q(√2) using regular representation:

1. E1 (standard relations): rank = 542 ✓
2. E2 (E1 + Kummer): rank = 544 ✓
3. E3 (E2 + all 5 Glanois relations): rank = 549 ✓

Key proofs:
- Each R_i is NOT in the Q(√2)-span of E2 (Step 4)
- All 5 relations are linearly independent modulo E2 (Step 5)

Step 6 (level-16) was not attempted due to matrix size constraints.
All computations were performed with exact integer arithmetic using Flint.
""")