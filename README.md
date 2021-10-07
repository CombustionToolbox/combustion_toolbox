# Combustion Toolbox
[![DOI](https://zenodo.org/badge/214183123.svg)](https://doi.org/10.5281/zenodo.5554912) ![repo size](https://img.shields.io/github/repo-size/AlbertoCuadra/combustion_toolbox) ![last modified](https://img.shields.io/github/last-commit/AlbertoCuadra/combustion_toolbox)

A MATLAB-GUI based open-source tool for solving combustion problems.

<!-- Website: https://combustiontoolbox.netlify.app/ -->
:top: There is also a (less complete) [Python version](https://github.com/AlbertoCuadra/Combustion-PyToolbox)

## Introduction
As a first step towards the development of a wider-scope thermochemical tool, in this work we present a thermochemical code with application to gaseous combustion problems recently implemented by the authors in MATLAB and Python. The MATLAB version solves six chemical equilibrium problems (`TP, HP, SP, TV, EV and SV transformations`; where T denotes temperature, P pressure, H enthalpy, S entropy, E internal energy and V volume), `incident and reflected planar shock waves`, as well as `ideal detonations according to Chapman-Jouguet theory and overdriven detonations`, assuming always ideal gases in all cases. The code also computes equilibrium properties of ideal plasmas, i.e., no coulombic interactions are considered. Moreover, part of our theoretical work "Thermochemical effects on hypersonic shock waves interacting with weak turbulence" is included. This allow us to compute from a theoretical perspective the jump conditions of an incident shock wave of a diatomic species (e.g., N2, O2, H2, F2) considering only dissociation of the species.

<!--
---
⚠️ **NOTE**

- At the moment, the Python version does not have all the capabilities that the MATLAB version has. I will continue with the development of this version adding all the remaining capabilities. I will also add a GUI using Qt6 and Pyside6.

---
-->

The code computes the equilibrium composition by minimization of the Gibbs–Helmholtz free energy by using Lagrange multipliers, and employs NASA’s 9-coefficient polynomial fits to evaluate the thermodynamic properties. Results computed with **Combustion Toolbox** have been validated against, and are in good agreement with, NASA’s Chemical Equilibrium with Applications (CEA) program, CANTERA and Caltech’s Shock and Detonation Toolbox. Along with the plain code, the new tool has been `equipped with a Graphical User Interface` developed in MATLAB 2021 under AppDesigner.

This project is also part of the PhD of [Alberto Cuadra-Lara](https://www.acuadralara.com/).

---
⚠️ **NOTE**

- At the moment, this is a preview of the tool and does **not** include the new GUI. The first final version v1.0.0 is expected to be released in late 2021. Check out the  upcoming features of [Combustion Toolbox v1.0.0](https://github.com/AlbertoCuadra/combustion_toolbox/projects/2).

---


## Examples

<p align="center">
    <img src="https://github.com/AlbertoCuadra/combustion_toolbox/blob/master/Validations/Hugoniot_benchmarking.svg" width="400">
</p>
    
<h1 align="center">Figure 1: Hugoniot curves for different molecular gases at pre-shock temperature T1 = 300 K and pressure p1 = 1 atm \[numerical results obtained with Combustion Toolbox (lines) and contrasted with NASA’s Chemical Equilibrium with Applications (CEA) code excluding ionization (symbols)\].</h1>
    
<p align="center">
    <img src="https://github.com/AlbertoCuadra/combustion_toolbox/blob/master/GUI/Snapshots/snapshot_1.svg" width="600">
</p>
<h1 align="center">Figure 2: Snapshot of the GUI.</h1>

Figure 1 is postprocesed with [Veusz](https://github.com/veusz/veusz) and [Inkscape](https://inkscape.org/). 

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

Please send feedback or inquiries to [acuadra@ing.uc3m.es](mailto:acuadra@ing.uc3m.es)

Thank you for testing Combustion Toolbox!

## Developers

* **[Alberto Cuadra-Lara](https://acuadralara.com/)** - *Main Developer and App designer*
* **César Huete** - *Developer*
* **Marcos Vera** - *Developer*

Grupo de Mecánica de Fluidos, Universidad Carlos III, Av. Universidad 30, 28911, Leganés, Spain

See also the list of [contributors](https://github.com/AlbertoCuadra/combustion_toolbox/blob/master/CONTRIBUTORS.md) who participated in this project.

## Citing Combustion Toolbox

```bibtex
@misc{combustiontoolbox,
    author = "Alberto Cuadra, César Huete and Marcos Vera",
    title = "Combustion Toolbox: A MATLAB-GUI based open-source tool for solving combustion problems",
    year = 2021,
    note = "Version 0.3.541",
    doi = {https://doi.org/10.5281/zenodo.5554912}
}
```
