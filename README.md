# Combustion Toolbox
A MATLAB MATLAB/GUI based thermochemical code
## Introduction
As a first step towards the development of a wider-scope thermochemical tool, in this work we present a thermochemical code with application to gaseous combustion problems recently implemented by the authors in MATLAB\textsuperscript{\textregistered}. The code solves six chemical equilibrium problems (TP, HP, SP, TV, EV and SV transformations; where T denotes temperature, P pressure, H enthalpy, S entropy, E internal energy and V volume), incident and reflected planar shock waves, as well as ideal detonations according to Chapman-Jouguet theory, assuming always ideal gases in all cases.\newline

The code computes the equilibrium composition using equilibrium constants rather than by minimization of the Gibbs–Helmholtz free energy, and employs NASA’s 9-coefficient polynomial fits to evaluate the thermodynamic properties. Along with the plain code, the new tool has been equipped with a Graphical User Interface (hereafter **Combustion-Toolbox**) developed in MATLAB\textsuperscript{\textregistered} 2018 under AppDesigner.\newline

Results computed with **Combustion-Toolbox** have been validated against, and are in good agreement with, NASA’s Chemical Equilibrium with Applications (CEA) program \cite{cite1}, CANTERA, and Caltech’s Shock and Detonation Toolbox \cite{sdtoolbox}. Moreover, the time required for the computations is comparable to that of other existing codes. **Combustion-Toolbox** has teaching and research aspirations and will be distributed as open source package as soon as it has been fully tested.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

MATLAB

```
Give examples
```

### Installing



## Running the tests



### Break down into end to end tests


## Built With

* [MATLAB-AppDesigner](https://es.mathworks.com/products/matlab/app-designer.html) - Used to developt the Graphical User Interface

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning



## Authors

* **Alberto Cuadra-Lara** - *Initial work* - [AlbertoCuadra](https://github.com/AlbertoCuadra)
* **Marcos Vera** - *Initial work*  

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License


## Acknowledgments

* Hat tip to anyone whose code was used
* Inspiration
* etc
