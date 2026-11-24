with open('../WT2_weight2_all_k/step_44.py', 'r') as f:
    content = f.read()
    # Get first 200 lines
    lines = content.split('\n')
    for i, line in enumerate(lines[:200]):
        print(f"{i+1:3}: {line}")