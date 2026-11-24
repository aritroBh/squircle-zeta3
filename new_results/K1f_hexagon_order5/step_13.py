
base = '<project>/research'
for f in ['10_frontier/E1_edge_law/verify/vjet.py', '10_frontier/E1_edge_law/verify/t3_hexagon.py']:
    print(f'===== {f} =====')
    with open(f'{base}/{f}') as fh:
        print(fh.read())
    print()
