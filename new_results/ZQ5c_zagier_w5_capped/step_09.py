base = '<project>/research/08_theory'
with open(f'{base}/ZQ5b_zagier_w5_larger_S/w5blib.py') as f:
    src = f.read()
print("TOTAL LEN:", len(src))
# print the truncated middle portion
print(src[9000:16500])
