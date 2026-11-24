import os
# We're already inside PG6_14gon_invariant. Look at parent.
parent = os.path.dirname(os.getcwd())
print("Parent:", parent)
for p in sorted(os.listdir(parent)):
    print(p)
