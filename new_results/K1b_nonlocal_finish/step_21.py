for i in [9,10,11,13,14,15]:
    fn = d + f"step_{i:02d}.py"
    try:
        with open(fn) as f:
            print("="*30, fn)
            print(f.read()[:2200])
    except FileNotFoundError:
        print(fn, "missing")
