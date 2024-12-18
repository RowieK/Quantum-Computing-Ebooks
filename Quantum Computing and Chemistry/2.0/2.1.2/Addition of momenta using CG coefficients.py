import sympy
from sympy import S
from sympy.physics.quantum.cg import CG, cg_simp

CG(S(1)/2, S(1)/2, S(1)/2, -S(1)/2, 1, 0).doit()
CG(S(1)/2, -S(1)/2, S(1)/2, S(1)/2, 1, 0).doit()

print(CG(S(1)/2, S(1)/2, S(1)/2, -S(1)/2, 1, 0).doit())
print(CG(S(1)/2, -S(1)/2, S(1)/2, S(1)/2, 1, 0).doit())

# The output of the code is:
# sqrt(2)/2
# sqrt(2)/2

CG(1, 0, S(1)/2, S(1)/2, S(1)/2, S(1)/2).doit()
CG(1, 1, S(1)/2, -S(1)/2, S(1)/2, S(1)/2).doit()
CG(1, -1, S(1)/2, S(1)/2, S(1)/2, S(1)/2).doit()

print(CG(1, 0, S(1)/2, S(1)/2, S(1)/2, S(1)/2).doit())
print(CG(1, 1, S(1)/2, -S(1)/2, S(1)/2, S(1)/2).doit())
