import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.special import sph_harm
import sympy
from sympy import S
from sympy.physics.quantum.cg import CG, cg_simp

# We define a function called setup_grid() that creates a grid of polar coordinates and 
# the corresponding cartesian coordinates with the following Python functions:
#  • numpy.linspace: Returns evenly spaced numbers over a specified interval
#  • numpy.meshgrid: Returns coordinate matrices from coordinate vectors
#  T
#  he setup_grid() function has one input parameter, num, which is a positive integer 
# that is the number of distinct values of polar coordinates. 

#  It returns the following:
#  • theta, phi: Two-dimensional NumPy arrays of shape num x num
#  • xyz: Three-dimensional NumPy array of shape (3,num,num)

def setup_grid(num=100):
    theta = np.linspace(0, np.pi, num)
    phi = np.linspace(0, 2*np.pi, num)
    # Create a 2D meshgrid from two 1D arrays of theta, phi coordinates
    theta, phi = np.meshgrid(theta, phi)
    # Compute the Cartesion coordintes with radius = 1
    xyz = np.array([np.sin(theta) * np.sin(phi),
                    np.sin(theta) * np.cos(phi),
                    np.cos(theta)])
    return (theta, phi, xyz)

(theta, phi, xyz) = setup_grid()
print("Shape of meshgrid arrays, theta: {}, phi: {}, xyz: {}".
format(theta.shape, phi.shape, xyz.shape))

def colour_plot(ax, Y, Yx, Yy, Yz, cmap):
    cmap = plt.cm.ScalarMappable(cmap=plt.get_cmap(cmap))
    cmap.set_clim(-0.5, 0.5)
    ax.plot_surface(Yx, Yy, Yz,
                    facecolors=cmap.to_rgba(Y.real),
                    rstride=2, cstride=2)
    return

def draw_axes(ax, ax_lim, title):
    ax.plot([-ax_lim, ax_lim], [0,0], [0,0], c='0.5', lw=1, zorder=10)
    ax.plot([0,0], [-ax_lim, ax_lim], [0,0], c='0.5', lw=1, zorder=10)
    ax.plot([0,0], [0,0], [-ax_lim, ax_lim], c='0.5', lw=1, zorder=10)
    # Set the limits, set the title and then turn off the axes 
    ax.set_title(title)
    ax.set_xlim(-ax_lim, ax_lim)
    ax.set_ylim(-ax_lim, ax_lim)
    ax.set_zlim(-ax_lim, ax_lim)
    ax.axis('off')
    return

def comb_Y(l, m, theta, phi):
    Y = sph_harm(abs(m), l, phi, theta)
    if m < 0:
        Y = np.sqrt(2) * (-1)**m * Y.imag
    elif m > 0:
        Y = np.sqrt(2) * (-1)**m * Y.real
    return Y

def plot_orbitals(k, cmap = 'autumn'):
  for l in range(0, k+1):
    for m in range(-l, l+1):
      fig = plt.figure(figsize=plt.figaspect(1.))
      (theta, phi, xyz) = setup_grid()
      ax = fig.add_subplot(projection='3d')
      Y = comb_Y(l, m, theta, phi)
      title = r'$l={{{}}}, m={{{}}}$'.format(l, m)
      Yx, Yy, Yz = np.abs(Y) * xyz
      colour_plot(ax, Y, Yx, Yy, Yz, cmap)
      draw_axes(ax, 0.5, title)
      fig_name = 'Hydrogen_l'+str(l)+'_m'+str(m)
      plt.savefig(fig_name)
      plt.show()
  return

plot_orbitals(2)

#------------------------------------------------------------
# The following code snippet is from Quantum%20Computing%20and%20Chemistry/2.0/2.1.2/all%20things%20combined.py
#------------------------------------------------------------
CG(S(1)/2, S(1)/2, S(1)/2, -S(1)/2, 1, 0).doit()
CG(S(1)/2, -S(1)/2, S(1)/2, S(1)/2, 1, 0).doit()

CG(1, 0, S(1)/2, S(1)/2, S(1)/2, S(1)/2).doit()
CG(1, 1, S(1)/2, -S(1)/2, S(1)/2, S(1)/2).doit()
CG(1, -1, S(1)/2, S(1)/2, S(1)/2, S(1)/2).doit()


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
        print("CG{} =".format(k), 'CG(', *T[k][0:6], ') * CG(', *T[k][4:10], ') = ', CGk.doit())
    return CGk

CG0 = comp_CG(T00, 0, display=True)
CG1 = comp_CG(T00, 1, display=True)
CG2 = comp_CG(T00, 2, display=True)
CG3 = comp_CG(T00, 3, display=True)
CG4 = comp_CG(T00, 4, display=True)
CG5 = comp_CG(T00, 5, display=True)

for k in range(0, len(T00)):
    s = 'CG' + str(k) +' = comp_CG(T00, ' + str(k) + ')'
    exec(s)
s00 = ["CG0: {}, CG1: {}, CG2: {}, CG3: {}, CG4: {}, CG5: {}".
       format(CG0.doit(), CG1.doit(), CG2.doit(), CG3.doit(), CG4.doit(), CG5.doit())]
print(s00)

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

#------------------------------------------------------------
def compute_3d_Y(ax_lim, cmap, title, fig_name):
    fig = plt.figure(figsize=plt.figaspect(1.))
    (theta, phi, xyz) = setup_grid()
    ax = fig.add_subplot(projection='3d')
    (Y20a, Y21a, Y2m1a, Y22a, Y2m2a) = (
        comb_Y(2, 0, theta, phi),
        comb_Y(2, 1, theta, phi),
        comb_Y(2, -1, theta, phi),
        comb_Y(2, 2, theta, phi),
        comb_Y(2, -2, theta, phi)
    )
    # Compute a sample linear combination of spherical harmonics
    Y = 0.2 * Y20a + 0.2 * Y21a + 0.2 * Y2m1a + 0.2 * Y22a + 0.2 * Y2m2a
    Yx, Yy, Yz = np.abs(Y) * xyz
    colour_plot(ax, Y, Yx, Yy, Yz, cmap)
    draw_axes(ax, ax_lim, title)
    plt.savefig(fig_name)
    plt.show()
    return

# Update the title and filename
title = '$Copper\ with\ 3d\ electrons\ (L=2,\ M=all)$'
fig_name = 'Copper_3d_L2_M_all.png'

compute_3d_Y(0.5, 'autumn', title, fig_name)

#------------------------------------------------------------
def compute_gold_wormhole(ax_lim, cmap, title, fig_name):
    fig = plt.figure(figsize=plt.figaspect(1.))
    (theta, phi, xyz) = setup_grid()
    ax = fig.add_subplot(projection='3d')

    # Create spherical harmonics for Gold's 5d orbitals
    Y20 = comb_Y(2, 0, theta, phi)
    Y21 = comb_Y(2, 1, theta, phi + np.pi / 6)
    Y2m1 = comb_Y(2, -1, theta, phi + np.pi / 6)
    Y22 = comb_Y(2, 2, theta, phi + np.pi / 3)
    Y2m2 = comb_Y(2, -2, theta, phi + np.pi / 3)

    # Weighted sum with phase shifts for a unique Gold visualization
    Y = 0.4 * Y20 + 0.3 * Y21 - 0.2 * Y2m1 + 0.1 * Y22 - 0.1 * Y2m2
    Yx, Yy, Yz = np.abs(Y) * xyz

    colour_plot(ax, Y, Yx, Yy, Yz, cmap)
    draw_axes(ax, ax_lim, title)
    plt.savefig(fig_name)
    plt.show()
    return

# Update the title and filename
title = '$Gold\ with\ customized\ 5d\ electrons\ (wormhole-like)$'
fig_name = 'Gold_5d_wormhole.png'

# Run the visualization
compute_gold_wormhole(0.5, 'plasma', title, fig_name)

#------------------------------------------------------------
def compute_gold_classic(ax_lim, cmap, title, fig_name):
    fig = plt.figure(figsize=plt.figaspect(1.))
    (theta, phi, xyz) = setup_grid()
    ax = fig.add_subplot(projection='3d')

    # Create spherical harmonics for Gold's 5d orbitals
    Y20 = comb_Y(2, 0, theta, phi)
    Y21 = comb_Y(2, 1, theta, phi)
    Y2m1 = comb_Y(2, -1, theta, phi)
    Y22 = comb_Y(2, 2, theta, phi)
    Y2m2 = comb_Y(2, -2, theta, phi)

    # Symmetric combination of harmonics
    Y = Y20 + Y21 + Y2m1 + Y22 + Y2m2
    Yx, Yy, Yz = np.abs(Y) * xyz

    # Plot the orbital
    colour_plot(ax, Y, Yx, Yy, Yz, cmap)
    draw_axes(ax, ax_lim, title)
    plt.savefig(fig_name)
    plt.show()
    return

# Set the title and filename
title = '$Gold\ with\ 5d\ electrons\ (classic\ symmetric\ visualization)$'
fig_name = 'Gold_5d_classic.png'

# Run the visualization
compute_gold_classic(0.5, 'inferno', title, fig_name)
