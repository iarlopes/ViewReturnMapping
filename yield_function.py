'''
This file defines the class "yield_functions" that is used in the package
ViewReturnMapping.
Different sub-classes are also defined here.

Igor A. Rodrigues Lopes, Jun 2021, Initial coding.
'''

import numpy as np
from utils import *

# Define the global class yield_functions
class yield_functions:
    """docstring for yield_functions"""

    def get_yield_function_value(self, Cauchy, sigma_y):
        pass

    def return_mapping(self, Cauchy, epbar):
        pass


class vonMises(yield_functions):
    """docstring for isotropic"""
    def __init__(self, kwargs):
        #super(isotropic, self).__init__()
        # Initial yielding stress
        self.sigma_y0 = kwargs['sigma_y0']

    def get_equivalent_stress(self, Cauchy):
        '''
        Determine the von Mises equivalent stress:
        sigma_eq = sqrt(3/2) norm(s)
        '''
        # Determine deviatoric component of the stress tensor
        s, p = additive_decomposition(Cauchy)
        # equivalent von Mises stress
        sigma_eq = np.sqrt(3.0/2.0) * np.linalg.norm(s)
        return sigma_eq

    def get_yield_function_value(self, Cauchy, sigma_y):
        '''
        Define the von Mises yield function:
        f = sqrt(3/2) norm(s) - sigma_y
        '''
        # Determine the equivalent von Mises stress
        sigma_eq = self.get_equivalent_stress(Cauchy)
        # Set the function value
        f = sigma_eq - sigma_y
        return f

    def plot_function(self, the_axes, sigma_y, color = None):
        """Plot von Mises yield function in the deviatoric plane.
           It is a circumference of radius sqrt(2/3)*sigma_y."""
        # Angle from 0 to 2*pi
        theta = np.linspace(0, 2*np.pi, num=100)
        x = np.sqrt(2/3) * sigma_y * np.cos(theta)
        y = np.sqrt(2/3) * sigma_y * np.sin(theta)
        if color == None:
            color = 'blue'
        return the_axes.plot(x,y,color = color)

    def return_mapping(self, Cauchy_old, epbar_old, DEtotal, Delast, E, mu, hard_law, miter, tol):
        '''
        Perform return mapping for a von Mises yield surface and associative flow rule
        '''
        # Perform elastic predictor and define elastic trial state
        DEtotal_vec = symMatrix2vec(DEtotal)
        DCauchy_vec = np.matmul(Delast, DEtotal_vec)
        Cauchy_trial = Cauchy_old + vec2symMatrix(DCauchy_vec)
        sigma_y0 = self.sigma_y0
        sigma_y = hard_law.get_current_sigma_y(sigma_y0, epbar_old)
        sigma_eq_trial = self.get_equivalent_stress(Cauchy_trial)
        f_trial = self.get_yield_function_value(Cauchy_trial, sigma_y)
        # outputs
        print('----------------------------------------------------')
        print('Trial state:')
        print('Cauchy trial = ')
        print(Cauchy_trial)
        print('Equivalent von Mises stress = ')
        print(sigma_eq_trial)
        print('von Mises yield function value = ')
        print(f_trial)
        
        # Perform return mapping if yield criteria is satisfied
        # solve the non-linear equation: sqrt(3/2)*norm(s_trial) - 3G*Dgamma - sigma_y = 0
        if (f_trial > 0):
            Dgamma = 0.0
            epbar = epbar_old + Dgamma
            res = f_trial
            iiter = 0
            G = E/(2*(1+mu))
            # loop over Newton-Raphson iterations
            while (iiter < miter and abs(res)/sigma_y0 > tol):
                iiter = iiter + 1
                dsigma_y_dDgamma = hard_law.derivative(sigma_y0, epbar)
                dres_dDgamma = -3*G - dsigma_y_dDgamma
                #
                dDgamma = - res / dres_dDgamma
                Dgamma = Dgamma + dDgamma
                epbar = epbar_old + Dgamma
                #
                sigma_y = hard_law.get_current_sigma_y(sigma_y0, epbar)
                res = sigma_eq_trial - 3*G*Dgamma - sigma_y
                # print iteration value
                print('----------------------------------------------------')
                print(f'Iteration nr. {iiter}')
                print(f'Dgamma^({iiter}) = {Dgamma}')
                print(f'res^({iiter}) = {res}')
            if (abs(res)/sigma_y0 > tol):
                print('The local iterative scheme did not converge!')

        # Updated stress and strain state
        # Increment of plastic strain
        s_trial, p_trial = additive_decomposition(Cauchy_trial)
        DEp = Dgamma*np.sqrt(3/2)*s_trial / np.linalg.norm(s_trial)
        DEp_vec = symMatrix2vec(DEp)
        # Increment of elastic strain
        DEe_vec = DEtotal_vec - DEp_vec
        # Increment of stress
        DCauchy_vec = np.matmul(Delast, DEe_vec)
        # Final stress
        Cauchy = Cauchy_old + vec2symMatrix(DCauchy_vec)
        # Update yield stress for current epbar
        sigma_y = hard_law.get_current_sigma_y(sigma_y0, epbar)
        # Update final values of equivalent stress and yield function value
        sigma_eq = self.get_equivalent_stress(Cauchy)
        f = self.get_yield_function_value(Cauchy, sigma_y)

        # print final values
        print('----------------------------------------------------')
        print('Final state:')
        print('epbar = ')
        print(epbar)
        print('Cauchy = ')
        print(Cauchy)
        print('Equivalent von Mises stress = ')
        print(sigma_eq)
        print('von Mises yield function value = ')
        print(f)
        #
        return Cauchy, epbar, Cauchy_trial

# Define the dictionary that relates yield function names to their sub-class
link_name = {'vonMises': vonMises}