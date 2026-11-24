
base='research/08_theory/CX1_octagon_depth2_proof'
src = open(base+'/verify.py').read()
# print the part that was truncated: from "def check_symbolic" to "def check_functional"
i = src.find('def check_symbolic')
j = src.find('def check_functional')
print(src[i:j])
