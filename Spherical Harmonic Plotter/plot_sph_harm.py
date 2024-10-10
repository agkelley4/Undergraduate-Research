import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import math
from matplotlib import cm

from coeffic_to_plot import readfile


# The following import configures Matplotlib for 3D plotting.
from mpl_toolkits.mplot3d import Axes3D
from scipy.special import sph_harm
plt.rc('text', usetex=True)

# Grids of polar and azimuthal angles
theta = np.linspace(0, 2*np.pi, 100)
phi = np.linspace(0, np.pi, 100)
# Create a 2-D meshgrid of (theta, phi) angles.
theta, phi = np.meshgrid(theta, phi)
# Calculate the Cartesian coordinates of each point in the mesh.
xyz = np.array([np.cos(theta) * np.sin(phi),
                np.sin(theta) * np.sin(phi),
                np.cos(phi)])


def get_Y(el, m, a):
#"""Plot the spherical harmonic of degree el and order m on Axes ax."""

# NB In SciPy's sph_harm function the azimuthal coordinate, theta,
# comes before the polar coordinate, phi.
    Y = sph_harm(abs(m), el, theta, phi)
    
    # Linear combination of Y_l,m and Y_l,-m to create the real form.
    if m < 0:
        Y = np.sqrt(2) * (-1)**m * Y.imag
    elif m > 0:
        Y = np.sqrt(2) * (-1)**m * Y.real

  
    #Y[0:24, :] = 0
    #Y[75:100, :] = 0
    Yx, Yy, Yz = Y * xyz * a

    return Yx, Yy, Yz, Y
   
def plot_Y(ax, Yx, Yy, Yz):

    # Colour the plotted surface according to the sign of Y.
    cmap = plt.cm.ScalarMappable(cmap=plt.get_cmap('copper'))
    cmap.set_clim(-0.5, 0.5)

    ax.plot_surface(Yx.real, 0, Yz.real, cmap='winter')


    # Draw a set of x, y, z axes for reference.
    ax_lim = 0.5
   # ax.plot([-ax_lim, ax_lim], [0,0], [0,0], c='0.5', lw=1, zorder=10)
   # ax.plot([0,0], [-ax_lim, ax_lim], [0,0], c='0.5', lw=1, zorder=10)
   # ax.plot([0,0], [0,0], [-ax_lim, ax_lim], c='0.5', lw=1, zorder=10)
    # Set the Axes limits and title, turn off the Axes frame.


    ax.set_title(r'${}$'.format(titlestring))
    ax_lim = .5
    ax.set_xlim(-ax_lim, ax_lim)
    ax.set_ylim(-ax_lim, ax_lim)
    ax.set_zlim(-ax_lim, ax_lim)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.axis('on')

#call read file from my code
el_values, m_values, alpha_values = readfile()
titlestring = []

for i in range(len(el_values)):
    l, m, a = el_values[i], m_values[i], alpha_values[i]
    print(l, m, a)

    Yx_i, Yy_i, Yz_i, Y_i = get_Y(l, m, a)
    if (i == 0) :
        Yx = Yx_i
        Yy = Yy_i
        Yz = Yz_i
        titlestring = '{}Y_{{{},{}}}'.format(round(a,2), int(l), int(m))

    else:
        Yx += Yx_i
        Yy += Yy_i
        Yz += Yz_i
        titlestring += ' + {}Y_{{{},{}}}'.format(round(a,2), int(l), int(m))
        



fig = plt.figure(figsize=plt.figaspect(1.))
ax = fig.add_subplot(projection='3d')    
plot_Y(ax, Yx, Yy, Yz)

plt.savefig('{}.png'.format(titlestring))
plt.show()


