# step_92: fix the removal properly
src = open("verify_pgen1.py").read()
i = src.find("        ]\n")  # the dangling closing bracket line
assert i > 0
# check context
print(repr(src[i-200:i+50]))
