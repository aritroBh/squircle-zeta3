import os
print(os.getcwd())
# Check the parent structure
for p in ['../GRW1_diamond_growth_borel', '../../GRW1_diamond_growth_borel']:
    if os.path.exists(p):
        print(p, os.listdir(p))
