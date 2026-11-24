import os, json
m5a='<project>/research/08_theory/M5a_next_coefficient'
m4e='<project>/research/08_theory/M4e_w5_lchi'
print(open(os.path.join(m5a,'n6_values.json')).read())
print(open(os.path.join(m4e,'Lchi8_5.json')).read()[:400])
print(open(os.path.join(m5a,'quad5_100.json')).read())
print(open(os.path.join(m5a,'W6_400.txt')).read())
