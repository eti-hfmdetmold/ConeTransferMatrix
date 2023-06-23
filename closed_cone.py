# -*- coding: utf-8 -*-
"""
@date: 30 March 2023

@author: Augustin Ernoult
"""

import numpy as np
import matplotlib.pyplot as plt

from openwind import ImpedanceComputation

"""
This file computes the input impedance curve of a weakly tapered, dissipative conical frustum
using a 1D finite-elements method.
It necessitates the library Openwind, freely available with pip or here: https://openwind.inria.fr
Computed with openwind [0.9.3]
"""

# %% constants and modelling options

# Frequencies of interest:
fmin, fmax, df = (40, 300, 0.5)
fs = np.arange(fmin, fmax+df, df)

model_options = dict(temperature=20,  # temperature in °C
                     spherical_waves=False, #Are the waves with spherical front? If False: plane waves are assumed
                     radiation_category = 'closed', # closed condition (0 flow) is imposed at the end of the pipe (or Z2=inf)
                     losses = True, # Zwikker and Kosten losses model is chosen
                     )

shape = [[0,  0.0021],[3, 0.0235]]# the geometry: length 3 m, input radius 2.1 mm, output radius 23.5 mm
#%% computations with FEM

l_ele = 5e-2 # size of elements in the mesh
order = 10 # order of the elements: it corresponds more or less to 10 nodes per elements

result_FEM = ImpedanceComputation(fs, shape, **model_options, l_ele=l_ele, order=order)

# Plot the instrument geometry
result_FEM.plot_instrument_geometry()

# Plot the impedance
fig = plt.figure()
result_FEM.plot_impedance(figure=fig, dbscale=False, normalize=False, label='FEM')
plt.show()

# Write the computed impedance in a file.
result_FEM.write_impedance('closed_cone(ow).txt')

# The output file reads
# first column: frequencies in Hz; second column: impedance, real part in Pa s/m^3; third column: impedance, imaginary part in Pa s/m^3 
