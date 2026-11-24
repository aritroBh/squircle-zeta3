import pickle, os
os.chdir('<project>/research/08_theory/ALPHA2_angle_leading_law')
with open('beta_tables.pkl','rb') as f: bt = pickle.load(f)
with open('beta_squircle.pkl','rb') as f: bs = pickle.load(f)
print(type(bt), list(bt.keys()) if hasattr(bt,'keys') else bt)
print(type(bs), list(bs.keys()) if hasattr(bs,'keys') else bs)
