
import pickle
fs = pickle.load(open('../COR8_depth_two_any_angle/full_sols_v3.pkl','rb'))
alln = fs['alln']
from collections import Counter
prods = [x for x in alln if x[0]=='prod']
print("sample prod entries:")
for x in prods[:15]: print("   ", x)
print("...")
li4s = [x for x in alln if x[0]=='Li4']
print("Li4 args:", [x[1] for x in li4s])
