base = '<project>/research/'
import json
tab = json.load(open(base+'08_theory/COR6_true_point_symbols/cor6_table.json'))
for row in tab:
    print(f"{row['label']:42s} deg={row['deg']} dimV={row['dimV']} span={row['letter_span_dim']} "
          f"true_zero={row['true_zero']} naive_zero={row['naive_zero']} differ={row['differ']} "
          f"ncomp_true={len(row['true_d22'])} time={row['time']}")
