
import pickle
h = pickle.load(open('cert_D7_0.pkl','rb'))
print(type(h))
items = list(h.items())[:5]
for k, v in items:
    print(k, type(v), repr(v)[:80])
