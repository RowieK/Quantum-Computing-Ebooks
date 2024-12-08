import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.special import sph_harm

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