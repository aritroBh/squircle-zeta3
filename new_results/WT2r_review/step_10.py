base = '<project>/research/08_theory/WT2b_stabilization'
with open(f'{base}/fastred.py') as f:
    src = f.read()
print(len(src))
# print the middle truncated part
print(src[4000:11000])
