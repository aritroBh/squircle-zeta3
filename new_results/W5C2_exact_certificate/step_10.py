import pickle, sys
with open('../W5C1_certificate_scoping/sector_meta.pkl','rb') as f:
    meta = pickle.load(f)
print(type(meta))
if isinstance(meta, dict):
    for k in list(meta)[:20]:
        v = meta[k]
        print(k, type(v), (len(v) if hasattr(v,'__len__') else v))