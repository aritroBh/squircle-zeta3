import pathlib

outputs_dir = pathlib.Path("<project>/paper/squircle_note/outputs")

# Read k=1,2,3 file
k123_file = outputs_dir / "ibp_reduce_k123.out"
if k123_file.exists():
    with open(k123_file, 'r') as f:
        content = f.read()
        print("First 500 chars of ibp_reduce_k123.out:")
        print(content[:500])
        print("\n" + "="*80 + "\n")

# Read k=4 file  
k4_file = outputs_dir / "ibp_reduce_k4.out"
if k4_file.exists():
    with open(k4_file, 'r') as f:
        content = f.read()
        print("First 500 chars of ibp_reduce_k4.out:")
        print(content[:500])