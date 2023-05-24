#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
from scipy.interpolate import griddata
from scipy.interpolate import splprep, splev
from scipy.interpolate import interp1d

# Define the input data
data = np.array([[ 0.6,  1.6,  1. ],
                        [ 0.2,  2. ,  0.5],
                        [ 0.2,  0.8,  0. ],
                        [ 0.8,  1. , -0.5],
                        [ 2. ,  1.6, -1. ],
                        [ 0.8,  0. , -1.5],
                        [ 3. ,  0.6, -2. ],
                        [ 2.8,  1. , -2.5],
                        [ 3.2,  3.2, -3. ],
                        [ 2. ,  0.2, -3.5],
                        [ 1. ,  1.6, -4. ],
                        [ 0.2,  1. , -4.5],
                        [ 2.2,  0.2, -5. ],
                        [ 1.8,  0.8, -5.5]])

x = data[:, 0]
y = data[:, 1]
z = data[:, 2]


max_z = np.max(z)
min_z = np.min(z)
tck, u = splprep([x, y, z], s=0)
u_new = np.linspace(0, 1, num=1000)
smoothed_coords = splev(u_new, tck)
smoothed_x, smoothed_y, smoothed_z = smoothed_coords

# Create a 3D plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x, y, z, c='b', label='Original Points')
ax.plot(smoothed_x, smoothed_y, smoothed_z, c='r', label='Smoothed Line')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_title('Smoothed Line Visualization')
ax.set_box_aspect([np.ptp(x), np.ptp(y), np.ptp(z)])
ax.legend()
plt.show()
