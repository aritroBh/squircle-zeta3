import json
tbl = json.load(open('<project>/research/08_theory/COR6_true_point_symbols/cor6_table.json'))
for e in tbl:
    if e['label'] in ('c=1/3', 'c=2/5', 'c=3/7'):
        print("=== ", e['label'], " true_d22:")
        for k, v in e['true_d22'].items():
            print("   ", k, v)
