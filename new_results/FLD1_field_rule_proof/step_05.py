for fn in ['flat_exact.py','flat_symbolic_proof.py','flat_single_zeta_k5.py','corner_series.py']:
    print("="*80); print("FILE:",fn); print("="*80)
    with open(f'<project>/research/10_frontier/G1_zeta3_any_angle/{fn}') as f:
        print(f.read())