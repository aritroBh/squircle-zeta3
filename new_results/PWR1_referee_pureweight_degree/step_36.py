import os
# locate D11 exact coefficients file
for cand in ['../D11_rhombus_resurgence/step1_exact.txt']:
    print(cand, os.path.exists(cand))
print(open('../D11_rhombus_resurgence/step1_exact.txt').read()[:1500])
