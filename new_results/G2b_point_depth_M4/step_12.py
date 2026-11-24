
base = "<project>/research/08_theory"
src = open(base + "/G2_M4alpha/step_10.py").read()
i = src.find("def load_terms")
print(src[i:i+2500])
print("======= g2_fit nz/neg =======")
src2 = open(base + "/G2_M4alpha/g2_fit.py").read()
for name in ["def nz", "def neg"]:
    j = src2.find(name)
    print(src2[j:j+400])
