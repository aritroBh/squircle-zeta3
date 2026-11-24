base='<project>/research/10_frontier/W1_corner_family/B_taylor_pi4/'
for fn in ['zsum.py','hpl.py']:
    print("="*80); print("FILE:",fn); print("="*80)
    with open(base+fn) as f:
        print(f.read())