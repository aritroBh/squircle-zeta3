src = open('../DEC2_decagon_exact_assembly/verify_dec2.py').read()
import re
# find section markers
for m in re.finditer(r'# ={3,}.*?={3,}|PART|Part', src):
    pass
print(src[15000:23000])