for i in [1,2,3,4,5]:
    fn = d + f"step_{i:02d}.py"
    with open(fn) as f:
        print("="*30, fn)
        print(f.read()[:1500])
