import os, sys, pathlib

# Look for paper/squircle_note/outputs
paper_base = "<project>/paper/squircle_note"
outputs_dir = pathlib.Path(paper_base) / "outputs"

if outputs_dir.exists():
    print(f"Files in {outputs_dir}:")
    for f in sorted(outputs_dir.iterdir()):
        print(f"  {f.name}")
else:
    print(f"Outputs directory doesn't exist: {outputs_dir}")
    
# Try to find from parent
print("\nSearching for reference files...")
for root, dirs, files in os.walk("<project>"):
    if "ibp_reduce_k3.out" in files or "ibp_reduce_k4.out" in files:
        print(f"Found in: {root}")
        for f in files:
            if f.startswith("ibp_reduce_k"):
                print(f"  {f}")