cwd = "<project>/research/08_theory/G2_M4alpha"
src = open(cwd + "/step_06.py").read()
i = src.find("f.write(f'# alpha")
print(repr(src[i-50:i+150]))
