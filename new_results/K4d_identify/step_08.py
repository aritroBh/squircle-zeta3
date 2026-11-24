import json
print("STORED COEFFICIENTS FROM K4c (key 'Sdk' = [d^k]S, k=0..4)")
print("="*80)
for f in ['modelA_h32.json','modelA_h48.json','modelA_h64.json','modelA_h96.json',
          'modelB_h32.json','modelB_h48.json','modelB_h96.json']:
    d = json.load(open(f))
    print("\n###", f, " keys:", list(d.keys()))
    for k,v in d.items():
        if k=='Sdk':
            print("  Sdk:")
            for i,val in enumerate(v):
                print("     [d^%d]S = %s"%(i,val))
        else:
            print("  %s = %s"%(k,v))
