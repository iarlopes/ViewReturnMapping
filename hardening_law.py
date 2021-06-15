'''
This file defines the class "hard_laws" that is used in the package
ViewReturnMapping.
Different sub-classes are also defined here.

Igor A. Rodrigues Lopes, Jun 2021, Initial coding.
'''


# Define the global class yield_functions
class hard_laws:
    """Objects if thus class must define a function describing the evolution of
       of the yield stress with the accumulated plastic strain, i.e., the hardening
       law. The function derivative must also be defined."""

    def get_current_sigma_y(self, sigma_y0, epbar):
        pass

    def derivative(self, sigma_y0, epbar):
        pass


class linear(hard_laws):
    """
    Defines the linear hardening law.
    It is expressed by:
    sigma_y = sigma_y0 + H*epbar
    """
    def __init__(self, kwargs):
        # Set hardening modulus
        self.H = kwargs['H']

    def get_current_sigma_y(self, sigma_y0, epbar):
        return sigma_y0 + self.H * epbar

    def derivative(self, sigma_y0, epbar):
        return self.H


class Swift(hard_laws):
    """
    Defines the Swift hardening law.
    It is expressed by:
    sigma_y = K*(e_0 + epbar)^n.
    K is set according to the initial yield stress:
    K = sigma_y0 / e_0^n
    """
    def __init__(self, kwargs):
        # Set hardening properties
        self.e_0 = kwargs['e_0']
        self.n   = kwargs['n']

    def get_current_sigma_y(self, sigma_y0, epbar):
        K = sigma_y0 / self.e_0**self.n
        return K * (self.e_0 + epbar)**self.n

    def derivative(self, sigma_y0, epbar):
        K = sigma_y0 / self.e_0**self.n
        return K * self.n * (self.e_0 + epbar)**(self.n-1)
        
# Define the dictionary that relates yield function names to their sub-class
link_name = {'Linear': linear,
             'Swift' : Swift}
