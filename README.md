# ViewReturnMapping
A Python Package to Perform and Visualise Return Mapping to Yield Surface

## The package

This package is built on an _object-oriented_ framework. Therefore, each material defined here is an object of the class _materials_.
The class _materials_, where the general methods are defined, is available in `material.py`. Material sub-classes are also defined in this file. At this moment, the subclass _isotropic_ is provided, but extensions to anisotropic materials are expected.

The yield function associated with a material is an object of the class _yield_functions_, defined in `yield_function.py`.
The von Mises yield criterion is defined in the sub-class _vonMises_, also defined in `yield_function.py`. The yield criterion name is linked to the corresponding subclass through the dictionary `link_name`, defined at the end of this file.

The treatment of the hardening law is similar. The class _hard_laws_ is defined in `hardening_law.py`, along with the sub-classes _linear_ and _Swift_, and the linking dictionary `link_name`.

Some util functions are provided in `util.py`.

## Example script and Jupyter notebook

A python script is provided in `example.py` as an example of the package usage.
Additionally, a Jupyter notebook is also available in `ReturnMappingNotebook.ipynb`, with the corresponding output being provided in `ReturnMappingNotebook.pdf`.
This Jupyter notebook is intended to be used for pedagogical purposes, in the framework of _Computational Plasticity_. The students can easily adjust the materials properties and strain/stress conditions to assess their impact on the inelastic behaviour.

## Usage
A brief guide describing the usage of this package is presented in what follows.

In the first place, the module _material_ must be imported to enable the access to the class _materials_ and corresponding sub-classes.
_NumPy_ arrays must be employed to define stress and strain states, and _Matplotlib_ is used to generate plots. Hence, these modules must also be imported:
```
import material
import numpy as np
import matplotlib.pyplot as plt
```

### Material Initialisation
To initialise an object for an isotropic material, the Young modulus and Poisson ratio must be defined:
```
my_material = material.isotropic(E=<Young>, mu=<Poisson>)
```
The elastic matrix is initialised according to these values. By default, the accumulated plastic strain and initial stress state are set to zero.

### Definition of Yield Function
An yield function is attached to the material through the method `set_yield_function(name, list_of_properties)`.

For the von Mises yield function, the name is _vonMises_ and the initial yield stress, _sigma_y0_, is the only required property:
```
my_material.set_yield_function('vonMises', sigma_y0=<YieldStressValue>)
```

### Definition of Hardening Law
A hardening law is attached to the material through the method `set_hardening_law(name, list_of_properties)`.

#### Linear hardening law
For the linear hardening law, the name is _Linear_ and the hardening modulus, _H_, must be defined:
```
my_material.set_hardening_law('Linear', H=<value>)
```

#### Swift hardening law
For the Swift hardening law, the name is _Swift_ and the parameters, _n_ and _e_0_, must be defined:
```
my_material.set_hardening_law('Swift', n=<value>, e_0=<value>)
```

### Definition of Initial State
In order to define an initial state other than null stress and accumulated plstic strain (initialised by default), the methods `set_initial_stress(stress_array)` and `set_initial_epbar(value)` may be employed:
```
my_material.set_initial_stress(stress_array)
my_material.set_initial_epbar(value)
```
Note that `stress_array` must be a _NumPy_ array with dimension 3x3.

These values may be recovered through
```
stress_initial = my_material.get_initial_stress()
epbar_initial = my_material.get_initial_epbar()
```
### Applying an Increment of Strain
An increment of strain may be applied to the material by
```
my_material.apply_increment_of_strain(strain_array)
```
with `strain_array` being a _NumPy_ array of dimension 3x3, representing the incremental strain tensor.

The resulting stress state, obtained considering elastic behaviour, can be recovered through
```
stress_current = my_material.get_current_stress()
```
### Specifying a Current Stress State
Alternatively, the current stress state may be specified through
```
my_material.set_current_stress(stress_array)
```
with `stress_array` being a _NumPy_ array of dimension 3x3, representing the specified stress tensor.

The incremental strain tensor leading to such stress state, considering elastic behaviour, can be recovered through
```
inc_strain = my_material.DE
```

### Performing Return Mapping
If the current stress state does not comply with the yield criterion, return mapping can be performed through
```
my_material.return_mapping(max_iter, tolerance)
```
The optional arguments `max_iter` and `tolerance` specify the maximum number of iterations and the convergence tolerance for the iterative procedure. If not defined, the default values are 20 and 1E-8, respectively.

The elastic trial stress and current stress tensor after return mapping can be recovered through
```
stress_trial = my_material.get_trial_stress()
stress_current = my_material.get_current_stress()
```

To check the equivalent stress associated with a given stress state, the following method may be employed
```
sigma_equivalent = my_material.get_equiv_stress(stress_array)
```

### Plotting Stress States and Yield Functions
To plot stress states and yield curves in the deviatoric plane, _Matplotlib_ figure and axes must be created in the first place.
```
fig, ax = plt.subplots(figsize = (<size_x>,<size_y>))
```
The yield curve is plotted through
```
my_material.plot_yield_function(ax, <yield_stress>, color = <string_color>)
```
A specific stress state can be represented in the figure through
```
my_material.plot_stress_state(ax, stress_array, color = <string_color>)
```
The principal stress axes vectors may be included in the representation by
```
my_material.plot_princ_stress_axis(ax, scaling = <value>)
```
I suggest that the scaling value is related to the yield stress value, for instance `1.2*yield_stress`, for better visualisation.

## Contacts
Send me an e-mail to ilopes@fe.up.pt if you have any question or comment.
Contributions and collaborations are welcomed.

I want to express my gratitude to [António Carneiro](https://github.com/amcc1996) for his valuable suggestions.
