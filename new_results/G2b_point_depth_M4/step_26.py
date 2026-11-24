
base = "<project>/research/08_theory"
s = open(base + "/G2_M4alpha/step_10.py").read()
# module-level statements (not inside def)
import ast
tree = ast.parse(s)
for node in tree.body:
    if isinstance(node, (ast.Import, ast.ImportFrom)):
        print(ast.dump(node)[:120])
    elif isinstance(node, ast.FunctionDef):
        print("def", node.name)
    else:
        print("MODULE-LEVEL:", type(node).__name__, ast.dump(node)[:200])
