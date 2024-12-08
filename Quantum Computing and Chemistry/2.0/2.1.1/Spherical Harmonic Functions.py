import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.special import sph_harm

def setup_grid(num=100):
    theta = np.linspace(0, np.pi, num)
    phi = np.linspace(0, 2*np.pi, num)
    #Create a 2D meshgrid from two 1D arrays of theta, phi coordinates
    theta, phi = np.meshgrid(theta, phi)
    #compute the Cartesian coordinates with radius = 1
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
                    facecolors=cmap.to.rgba(Y.real),
                    rstride=2, cstride=2)
    return