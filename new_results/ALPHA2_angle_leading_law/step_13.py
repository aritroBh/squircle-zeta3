import pickle
with open('beta_tables.pkl','rb') as f: bt = pickle.load(f)
for cname in ['c12','c35','c13','c45']:
    tab = bt[cname]
    print(f"===== {cname} ===== keys sample:", list(tab.keys())[:3])
