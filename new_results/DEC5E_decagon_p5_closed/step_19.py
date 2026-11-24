import pickle, subprocess
print(open('../DEC5D_decagon_order5_depth/p5_run170.pkl.txt').read())
print(open('../DEC5D_decagon_order5_depth/verify_p5.pkl.txt').read())
lb = pickle.load(open('../DEC5D_decagon_order5_depth/local_block_exact.pkl','rb'))
print(type(lb), list(lb.keys()) if isinstance(lb, dict) else len(lb))
gold = pickle.load(open('../DEC5D_decagon_order5_depth/p5_gold.pkl','rb'))
print(gold.keys())
out = subprocess.run(['grep','-rn','73125','../DEC5D_decagon_order5_depth/step_*.py'],capture_output=True,text=True)
print(out.stdout[:3000])