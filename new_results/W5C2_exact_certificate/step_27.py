import subprocess, sys, pickle
print(subprocess.run(['bash','-c','ls ../A3c_claude_elementary/Q/ | head -30; find .. -maxdepth 2 -name "qsqrt2*" 2>/dev/null | head'], capture_output=True, text=True).stdout)
sys.path.insert(0, '../A3c_claude_elementary/Q')
vw = pickle.load(open('../W5C1_certificate_scoping/validated_w5.pkl','rb'))
print('validated_w5:', type(vw), list(vw.keys()) if isinstance(vw, dict) else len(vw))