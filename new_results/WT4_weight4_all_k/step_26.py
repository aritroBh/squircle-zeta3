import os, json
base='<project>/research/08_theory'
print(json.load(open(os.path.join(base,'M4a_decomposition/hp_results.json'))))
