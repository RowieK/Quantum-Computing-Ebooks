import sympy
from sympy import S
from sympy.physics.quantum.cg import CG, cg_simp
import matplotlib.pyplot as plt
from matplotlib import cm

T00 = {0: (1,-1, 1,0,  1,-1, 1,1,  0,0), 
       1: (1,-1, 1,1,  1,0,  1,0,  0,0),
       2: (1,0,  1,-1, 1,-1, 1,1,  0,0),
       3: (1,0,  1,1,  1,1,  1,-1, 0,0),
       4: (1,1,  1,-1, 1,0,  1,0,  0,0),
       5: (1,1,  1,0,  1,1,  1,-1, 0,0)}
 
def comp_CG(T, k, display = None):
    CGk = CG(*T[k][0:6]) * CG(*T[k][4:10])
    if display:
        print('CG(', *T[k][0:6], ') = ', CG(*T[k][0:6]).doit())
        print('CG(', *T[k][4:10], ') = ', CG(*T[k][4:10]).doit())
        print("CG{} =".format(k), 'CG(', *T[k][0:6], ') * CG(', 
    *T[k][4:10], ') = ', CGk.doit())
    return CGk

for k in range(0, len(T00)):
    s = 'CG' + str(k) +' = comp_CG(T00, ' + str(k) + ')'
    exec(s)
s00 = ["CG0: {}, CG1: {}, CG2: {}, CG3: {}, CG4: {}, CG5: {}".
    format(CG0.doit(), CG1.doit(), CG2.doit(), CG3.doit(), CG4.doit(), CG5.doit())]
print(s00)

CG0 = comp_CG(T00, 0, display=True)

def Y_phase(theta, phi):
    Y10a = comb_Y(1, 0, theta, phi)
    Y11a = comb_Y(1, 1, theta, phi)
    Y1m1a = comb_Y(1, -1, theta, phi)
    Y10b = comb_Y(1, 0, theta, phi+1*np.pi/3)
    Y11b = comb_Y(1, 1, theta, phi+1*np.pi/3)
    Y1m1b = comb_Y(1, -1, theta, phi+1*np.pi/3)
    Y10c = comb_Y(1, 0, theta, phi+2*np.pi/3)
    Y11c = comb_Y(1, 1, theta, phi+2*np.pi/3)
    Y1m1c = comb_Y(1, -1, theta, phi+2*np.pi/3)
    return(Y10a, Y11a, Y1m1a, Y10b, Y11b, Y1m1b, Y10c, Y11c, Y1m1c)

def compute_00_Y(ax_lim, cmap, title,  fig_name):
  fig = plt.figure(figsize=plt.figaspect(1.))
  (theta, phi, xyz) = setup_grid()
  ax = fig.add_subplot(projection='3d')
  (Y10a, Y11a, Y1m1a, Y10b, Y11b, Y1m1b, Y10c, Y11c, Y1m1c) = Y_phase(theta, phi)
  Y_00 = float(CG0.doit()) * Y1m1a * Y10b * Y11c
  Y_01 = float(CG1.doit()) * Y1m1a * Y11b * Y10c
  Y_02 = float(CG2.doit()) * Y10a * Y1m1b * Y11c
  Y_03 = float(CG3.doit()) * Y10a * Y11b * Y1m1c
  Y_04 = float(CG4.doit()) * Y11a * Y1m1b * Y10c
  Y_05 = float(CG5.doit()) * Y11a * Y10b * Y1m1c
  Y = Y_00 + Y_01 + Y_02 + Y_03 + Y_04 + Y_05
  Yx, Yy, Yz = np.abs(Y) * xyz
  colour_plot(ax, Y, Yx, Yy, Yz, cmap)
  draw_axes(ax, ax_lim, title)
  plt.savefig(fig_name)
  plt.show()
  return

title = '$Nitrogen\ with\ 3p\ electrons\ (L=0,\ M=0)$'
fig_name ='Nitrogen_3p_L0_M0.png'
compute_00_Y(0.01, 'autumn', title, fig_name)
