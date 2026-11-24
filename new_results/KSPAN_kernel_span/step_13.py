import pathlib, re

outputs_dir = pathlib.Path("<project>/paper/squircle_note/outputs")

# Parse k=3 from k123 file
k123_file = outputs_dir / "ibp_reduce_k123.out"
if k123_file.exists():
    with open(k123_file, 'r') as f:
        content = f.read()
    
    # Find k=3 section
    k3_start = content.find("===== k = 3:")
    if k3_start >= 0:
        # Find next k section or end
        k3_end = content.find("===== k =", k3_start + 1)
        k3_section = content[k3_start:k3_end] if k3_end >= 0 else content[k3_start:]
        
        print("k=3 section (first 1000 chars):")
        print(k3_section[:1000])
        print("\nParsing degree-2 pieces...")
        
        # Parse format: weight 4 l1^2 l2^1 or l1^3 l2^0
        lines = k3_section.split('\n')
        for line in lines:
            if 'weight 4' in line:
                print(line.strip())