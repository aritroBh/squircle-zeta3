with open('../WT2_weight2_all_k/fastred.py', 'r') as f:
    # Get first 100 lines
    content = f.read().split('\n')
    for i, line in enumerate(content[:100]):
        print(f"{i+1:3}: {line}")