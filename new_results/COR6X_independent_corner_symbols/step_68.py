import json, os
tbl = json.load(open('<project>/research/08_theory/COR6_true_point_symbols/cor6_table.json'))
print(f"{'label':28s} {'deg':>3s} {'dimV':>4s} {'lspan':>5s} {'true0':>6s} {'naive0':>6s} {'differ':>6s}")
for e in tbl:
    print(f"{e['label']:28s} {e['deg']:>3d} {e['dimV']:>4d} {e['letter_span_dim']:>5d} {str(e['true_zero']):>6s} {str(e['naive_zero']):>6s} {str(e['differ']):>6s}")
print("total entries:", len(tbl))
