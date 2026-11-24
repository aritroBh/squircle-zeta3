c2 = open('../K3_fourier_all_orders/result.md').read()
# print section 1.3 remainder and section 2 start
i = c2.find('odd tops match')
print(c2[i:i+4500])
