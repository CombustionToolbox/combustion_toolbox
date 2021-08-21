# Combustion Toolbox
![repo size](https://img.shields.io/github/repo-size/AlbertoCuadra/combustion_toolbox) ![last modified](https://img.shields.io/github/last-commit/AlbertoCuadra/combustion_toolbox)

A MATLAB-GUI based open-source tool for solving combustion problems.

Website: https://combustiontoolbox.netlify.app/

Table of contents
=================

<!--ts-->
   * [Introduction](#Introduction)
   * [Getting Started](#Getting-Started)
   * [Prerequisites](#Prerequisites)
   * [Installing](#Installing)
   * [Running the tests](#Running-the-tests)
   * [Built With](#Built-With)
   * [Contributing](#Contributing)
   * [Versioning](#Versioning)
   * [Authors](#Authors)
   * [License](#License)
   * [Acknowledgments](#Acknowledgments)
   
<!--te-->

## Introduction
As a first step towards the development of a wider-scope thermochemical tool, in this work we present a thermochemical code with application to gaseous combustion problems recently implemented by the authors in MATLAB. The code solves six chemical equilibrium problems (TP, HP, SP, TV, EV and SV transformations; where T denotes temperature, P pressure, H enthalpy, S entropy, E internal energy and V volume), incident and reflected planar shock waves, as well as ideal detonations according to Chapman-Jouguet theory and overdriven detonations, assuming always ideal gases in all cases.

The code computes the equilibrium composition using equilibrium constants rather than by minimization of the Gibbs–Helmholtz free energy, and employs NASA’s 9-coefficient polynomial fits to evaluate the thermodynamic properties. Along with the plain code, the new tool has been equipped with a Graphical User Interface (hereafter **Combustion Toolbox**) developed in MATLAB 2018 under AppDesigner.

Results computed with **Combustion Toolbox** have been validated against, and are in good agreement with, NASA’s Chemical Equilibrium with Applications (CEA) program, CANTERA, and Caltech’s Shock and Detonation Toolbox. Moreover, the time required for the computations is comparable to that of other existing codes. **Combustion Toolbox** has teaching and research aspirations and will be distributed as open source package as soon as it has been fully tested.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

MATLAB

```
Give examples
```

### Installing



## Running the tests


## Built With

* [MATLAB-AppDesigner](https://www.mathworks.com/products/matlab/app-designer.html) - Used to developt the Graphical User Interface

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning



## Authors

* **[Alberto Cuadra-Lara](https://albertocuadra.netlify.app/)<sup>1</sup>** - *Main Developer and App designer*
* **Marcos Vera<sup>1</sup>** - *Developer*  
<sup>1</sup>  Grupo de Mecánica de Fluidos, Universidad Carlos III, Av. Universidad 30, 28911, Leganés, Spain

See also the list of [contributors](https://github.com/AlbertoCuadra/combustion_toolbox/blob/master/CONTRIBUTORS.md) who participated in this project.

## License


## Acknowledgments

A.C. is deeply grateful to Samuel Delbarre for his help in the validation and debuging process.
