
with open('step_03.py') as f:
    txt = f.read()
# fix outward orientation in fnrm
txt = txt.replace("""def fnrm(a,b,c):
    n=np.cross(b-a,c-a); return n/np.linalg.norm(n)""",
"""def fnrm(a,b,c):
    n=np.cross(b-a,c-a); n=n/np.linalg.norm(n)
    if np.dot(n,(a+b+c)/3) < 0: n = -n
    return n""")
# replace sympy Triangle area with cross product
txt = txt.replace("""import sympy.geometry as sg
# triangle (A1,P1,Q1) exact
a1=sp.Point3D(0,0,1/(2*sqrt(2))); p1=sp.Point3D(1/(3*sqrt(2)),1/(3*sqrt(2)),1/(3*sqrt(2)))
q1=sp.Point3D(1/(3*sqrt(2)),-1/(3*sqrt(2)),1/(3*sqrt(2)))
ar = sg.Triangle(a1,p1,q1).area""",
"""a1=(0,0,1/(2*sqrt(2))); p1=(1/(3*sqrt(2)),1/(3*sqrt(2)),1/(3*sqrt(2)))
q1=(1/(3*sqrt(2)),-1/(3*sqrt(2)),1/(3*sqrt(2)))
u=[p1[i]-a1[i] for i in range(3)]; v=[q1[i]-a1[i] for i in range(3)]
cr=[u[1]*v[2]-u[2]*v[1], u[2]*v[0]-u[0]*v[2], u[0]*v[1]-u[1]*v[0]]
ar = sp.sqrt(sum(c**2 for c in cr))/2""")
with open('step_03.py','w') as f:
    f.write(txt)
import subprocess
r = subprocess.run(['<home>)/.venv/bin/python','step_03.py'],
                   capture_output=True, text=True)
print(r.stdout)
print(r.stderr[-3000:] if r.returncode else "exit 0")
