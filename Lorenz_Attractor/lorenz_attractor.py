"""
# =============================================================================
# The Lorenz Attractor
# =============================================================================

This is a simple code that simulates the famous Lorenz Model.
I simulate it to collect experience in Python

"""
# =============================================================================
# Keep the house clean
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# Define a function and call it "lorenz"
# =============================================================================
def lorenzAt(x, y, z, s=10, r=28, b=2.667):
    '''
    Given:
       x, y, z: a point of interest in three dimensional space
       s, r, b: parameters defining the lorenz attractor
    Returns:
       x_dot, y_dot, z_dot: values of the lorenz attractor's partial
           derivatives at the point x, y, z
    '''
    
# Define the ODEs    
    x_dot = s*(y - x)
    y_dot = r*x - y - x*z
    z_dot = x*y - b*z
    return x_dot, y_dot, z_dot


# Determine time increment and time horizon
dt = 0.01
num_steps = 10000

# Need one more for the initial values
xs = np.empty(num_steps + 1)
ys = np.empty(num_steps + 1)
zs = np.empty(num_steps + 1)

# Set initial values
xs[0], ys[0], zs[0] = (0., 1., 1.05)

# Step through "time", calculating the partial derivatives at the current point
# and using them to estimate the next point

# =============================================================================
# Finite Difference Method
# =============================================================================
for i in range(num_steps):
    x_dot, y_dot, z_dot = lorenzAt(xs[i], ys[i], zs[i])
    xs[i + 1] = xs[i] + (x_dot * dt)
    ys[i + 1] = ys[i] + (y_dot * dt)
    zs[i + 1] = zs[i] + (z_dot * dt)
# Note that the function output is obtained without the use of paranthesis

# Plot as function

def fifug(xx, yy, zz):
    fig = plt.figure()
    ax1 = fig.gca(projection='3d')
    ax1.plot(xx, yy, zz, lw=0.5)
    ax1.set_xlabel("X Axis")
    ax1.set_ylabel("Y Axis")
    ax1.set_zlabel("Z Axis")
    ax1.set_xlim([-20, 20])
    ax1.set_ylim([-20, 20])
    ax1.set_title("Lorenz Attractor")
    return fig
fifug(xs, ys, zs)   
