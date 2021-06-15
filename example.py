'''
This script is employed to test the features of the package ViewReturnMapping.

Igor A. Rodrigues Lopes, Jun 2021, Initial coding.
'''

import material
import matplotlib.pyplot as plt

# Define the material steel
# initialise with elastic properties
steel = material.isotropic(E=210e9, mu=0.3)
# set the von Mises yield surface
steel.set_yield_function('vonMises', sigma_y0 = 240e6)
# set the hardening law properties.
# for Linear hardening only H is needed: sigma_y = sigma_y0 + H*epbar
# for Swift hardening, define e_0 and n: sigma_y = K*(e_0 + epbar)^n, with K = sigma_y0/(e_0^n)
steel.set_hardening_law('Swift', n=0.15, e_0=0.0005)
# Define a stress tensor (must be numpy array)
import numpy as np
stress = 1.0e6 * np.array([[140, 120, 0],
	                       [120, 40 , 0],
	                       [0  ,0   ,20]])
# Get equivalent stress for this stress state
print(steel.get_equiv_stress(stress))
# Set this stress state as the initial material state
steel.set_initial_stress(1e6*np.array([[140,120,0],[120,40,0],[0,0,20]]))
# and set the initial accumulated plastic strain as zero
steel.set_initial_epbar(0.0)
# Define an incremental strain tensor
DE=np.zeros([3,3])
DE[0,0] = 0.0002
# Apply it to the material
steel.apply_increment_of_strain(DE)
# The current stress after applying the increment of strain, assuming elastic strain
stress_cur = steel.get_current_stress()
# And the corresponding equivalent stress is
print(steel.get_equiv_stress(stress_cur))
# It is not admissible. Therefore lets perform return mapping
steel.return_mapping()
#
# Now lets create a figure to plot the initial yield function
fig, ax = plt.subplots(figsize = (5,5.5))
steel.plot_yield_function(ax, 240e6, color = 'blue')
# Get the old stress, trial stress and current stress
stress_old = steel.get_initial_stress()
stress_trial = steel.get_trial_stress()
stress_new = steel.get_current_stress()
# Plot the stress states
steel.plot_stress_state(ax, stress_old, color = 'black')
steel.plot_stress_state(ax, stress_trial, color = 'red')
steel.plot_stress_state(ax, stress_new, color = 'green')
ax.axis('equal')
ax.axis('off')
steel.plot_princ_stress_axis(ax, scaling = 1.2*240e6)
ax.legend(['Initial yield surface', 'Initial stress state', 'Trial stress state', 'Current stress state'], loc = 'lower center', ncol = 2)
plt.show()
