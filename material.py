'''
This file defines the class "materials" that is used in the package
ViewReturnMapping.
Different sub-classes are also defined here.

Igor A. Rodrigues Lopes, Jun 2021, Initial coding.
'''

#from abc import ABC, abstractmethod 
import numpy as np
import yield_function
import hardening_law
from utils import *

# Define the global class material
class materials:
    """docstring for materials"""

    def set_yield_function(self, name, **kwargs):
        self.yield_func = yield_function.link_name[name](kwargs)

    def set_hardening_law(self, name, **kwargs):
        self.hard_law = hardening_law.link_name[name](kwargs)

    def set_initial_stress(self, Cauchy_old):
        self.Cauchy_old = Cauchy_old

    def set_current_stress(self, Cauchy):
        self.Cauchy = Cauchy
        DCauchy = Cauchy - self.Cauchy_old
        DCauchy_vec = symMatrix2vec(DCauchy)
        DE_vec = np.matmul(np.linalg.inv(self.Delast), DCauchy_vec)
        self.DE = vec2symMatrix(DE_vec)

    def set_initial_epbar(self, epbar_old):
        self.epbar_old = epbar_old

    def get_equiv_stress(self, Cauchy):
        '''
        Returns the equivalent stress for a given stress state.
        '''
        sigma_eq = self.yield_func.get_equivalent_stress(Cauchy)
        return sigma_eq

    def get_current_stress(self):
        '''
        Returns the current stress state.
        '''
        return self.Cauchy

    def get_initial_stress(self):
        '''
        Returns the current stress state.
        '''
        return self.Cauchy_old

    def get_trial_stress(self):
        '''
        Returns the trial stress state.
        '''
        return self.Cauchy_trial

    def get_initial_epbar(self):
        '''
        Returns the initial accumulated plastic strain.
        '''
        return self.epbar_old

    def get_current_epbar(self):
        '''
        Returns the current accumulated plastic strain.
        '''
        return self.epbar

    def apply_increment_of_strain(self, DE):
        '''
        Applies an increment of strain to the material and updates the current stress state
        according to the elastic law.
        '''
        # get the increment of stress resulting from elastic relation
        DE_vec = symMatrix2vec(DE)
        DCauchy_vec = np.matmul(self.Delast, DE_vec)
        DCauchy = vec2symMatrix(DCauchy_vec)
        self.Cauchy = self.Cauchy_old + DCauchy
        self.DE = DE

    def return_mapping(self, max_iter=20, tolerance=1e-8):

        self.Cauchy, self.epbar, self.Cauchy_trial = \
            self.yield_func.return_mapping(self.Cauchy_old,
                                           self.epbar_old,
                                           self.DE,
                                           self.Delast,
                                           self.E,
                                           self.mu,
                                           self.hard_law,
                                           max_iter,
                                           tolerance)

    def plot_yield_function(self, the_axes, sigma_y, color = None):
        return self.yield_func.plot_function(the_axes, sigma_y, color)

    def plot_stress_state(self, the_axes, Cauchy, color = None):
        if color == None:
            color = 'blue'
        # Get the vector with principal stresses
        princ = np.linalg.eigvals(Cauchy)
        princ_dev_plane = coordsPiPlane(princ)
        # Plot the point
        return the_axes.scatter(princ_dev_plane[0], princ_dev_plane[1], color = color)

    def plot_princ_stress_axis(self, the_axes, scaling):
        # Define the vectors with principal stress directions
        sigma_1 = np.array([1,0,0])
        sigma_2 = np.array([0,1,0])
        sigma_3 = np.array([0,0,1])
        #
        Tmatrix = 1/np.sqrt(6) * np.array([[np.sqrt(3), 0         , -np.sqrt(3)],
                                           [-1        , 2         , -1],
                                           [np.sqrt(2),-np.sqrt(2), np.sqrt(2)]])
        s1_dev, s2_dev, s3_dev = coordsPiPlane(sigma_1), coordsPiPlane(sigma_2), coordsPiPlane(sigma_3)
        # Array with vectors to plot and origin
        origin = np.array([[0,0,0],[0,0,0]])
        vecs = np.array([[s1_dev[0], s1_dev[1]],
                         [s2_dev[0], s2_dev[1]],
                         [s3_dev[0], s3_dev[1]]])
        # Plot axis vectors
        return the_axes.quiver(*origin, vecs[:,0], vecs[:,1], color = 'black', scale = 1/scaling, scale_units = 'xy'),\
               the_axes.text(s1_dev[0]*scaling, s1_dev[1]*scaling, '$\\sigma_1$', usetex = True, family = 'serif', size = 18),\
               the_axes.text(s2_dev[0]*scaling, s2_dev[1]*scaling, '$\\sigma_2$', usetex = True, family = 'serif', size = 18),\
               the_axes.text(s3_dev[0]*scaling, s3_dev[1]*scaling*1.1, '$\\sigma_3$', usetex = True, family = 'serif', size = 18)



class isotropic(materials):
    """docstring for isotropic"""
    def __init__(self, **kwargs):
        #super(isotropic, self).__init__()
        # Young modulus
        self.E = E = kwargs['E']
        # Poisson ratio
        self.mu = mu = kwargs['mu']
        # Define the elastic matrix
        self.Delast = E/(1+mu)/(1-2*mu) * np.array([[1-mu, mu, mu, 0, 0, 0],
                                                    [mu, 1-mu, mu, 0, 0, 0],
                                                    [mu, mu, 1-mu, 0, 0, 0],
                                                    [0, 0, 0, 1-2*mu, 0, 0],
                                                    [0, 0, 0, 0, 1-2*mu, 0],
                                                    [0, 0, 0, 0, 0, 1-2*mu]])
        # set the initial stress state to zero
        self.Cauchy_old = self.Cauchy = self.Cauchy_trial = np.array(np.zeros([3,3]))
        self.epbar_old = self.epbar = 0.0
        self.DE = np.zeros([3,3])
