for f in ['modelB_h48.json','modelB_h96.json','modelA_h96.json']:
    print("="*60)
    print(f, "RAW content:")
    raw = open(f).read()
    print(raw[:1200])
    print()
