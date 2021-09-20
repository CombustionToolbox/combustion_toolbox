# Combustion Toolbox
![repo size](https://img.shields.io/github/repo-size/AlbertoCuadra/combustion_toolbox) ![last modified](https://img.shields.io/github/last-commit/AlbertoCuadra/combustion_toolbox)

A MATLAB-GUI based open-source tool for solving combustion problems.

Website: https://combustiontoolbox.netlify.app/


## Introduction
As a first step towards the development of a wider-scope thermochemical tool, in this work we present a thermochemical code with application to gaseous combustion problems recently implemented by the authors in MATLAB and Python. The MATLAB version solves six chemical equilibrium problems (`TP, HP, SP, TV, EV and SV transformations`; where T denotes temperature, P pressure, H enthalpy, S entropy, E internal energy and V volume), `incident and reflected planar shock waves`, as well as `ideal detonations according to Chapman-Jouguet theory and overdriven detonations`, assuming always ideal gases in all cases.

---
⚠️ **NOTE**

- The overdriven detonation routine are based on Caltech’s Shock and Detonation Toolbox [Browne et. al (2008)](https://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.181.1150&rank=1&q=Numerical%20Solution%20Methods%20for%20Shock%20and%20Detonation%20Jump%20Condit&osm=&ossid=). This is going to be replaced in a specific in-house set of routines based on NASA's reports.
- At the moment, the Python version does not have all the capabilities that the MATLAB version has. I will continue with the development of this version adding all the remaining capabilities. I will also add a GUI using Qt6 and Pyside6.

---

The code computes the equilibrium composition by minimization of the Gibbs–Helmholtz free energy by using Lagrange multipliers, and employs NASA’s 9-coefficient polynomial fits to evaluate the thermodynamic properties. Results computed with **Combustion Toolbox** have been validated against, and are in good agreement with, NASA’s Chemical Equilibrium with Applications (CEA) program, CANTERA and Caltech’s Shock and Detonation Toolbox. Along with the plain code, the new tool has been `equipped with a Graphical User Interface` developed in MATLAB 2021 under AppDesigner.

This project is also part of my PhD.

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Authors

* **[Alberto Cuadra-Lara](https://albertocuadra.netlify.app/)<sup>1</sup>** - *Main Developer and App designer*
* **Marcos Vera<sup>1</sup>** - *Developer*  
<sup>1</sup>  Grupo de Mecánica de Fluidos, Universidad Carlos III, Av. Universidad 30, 28911, Leganés, Spain

See also the list of [contributors](https://github.com/AlbertoCuadra/combustion_toolbox/blob/master/CONTRIBUTORS.md) who participated in this project.

## Acknowledgments

A.C. is deeply grateful to Samuel Delbarre for his help in the validation and debuging process of the segregated method (deprecated).
