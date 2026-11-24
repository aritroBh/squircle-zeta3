import pathlib, re

k4_file = outputs_dir / "ibp_reduce_k4.out"
if k4_file.exists():
    with open(k4_file, 'r') as f:
        content = f.read()
    
    print("k=4 section (first 1500 chars):")
    print(content[:1500])
    
    # Parse degree-2 pieces (weight 4)
    print("\n\nParsing degree-2 pieces (weight 4) for k=4:")
    lines = content.split('\n')
    for i, line in enumerate(lines):
        if 'weight 4' in line:
            # Get the weight 4 section
            key_line = line.strip()
            print(key_line)
            # Print next few lines until empty line or next weight
            j = i + 1
            while j < len(lines) and lines[j].strip() and 'weight' not in lines[j]:
                print(lines[j].strip())
                j += 1
            print()