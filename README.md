<p align="left">
    <img alt="UC3M" style="border-width:0" src="https://github.com/AlbertoCuadra/combustion_toolbox/blob/master/GUI/Icons/logo.svg" width="1500"/></a>
</p>

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6369301.svg)](https://doi.org/10.5281/zenodo.6369301)
[![File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://es.mathworks.com/matlabcentral/fileexchange/101088-combustion-toolbox)
[![CI](https://github.com/AlbertoCuadra/combustion_toolbox/actions/workflows/CI.yml/badge.svg)](https://github.com/AlbertoCuadra/combustion_toolbox/actions/workflows/CI.yml)
[![CD](https://github.com/AlbertoCuadra/combustion_toolbox/actions/workflows/CD.yml/badge.svg)](https://github.com/AlbertoCuadra/combustion_toolbox/actions/workflows/CD.yml)
![last modified](https://img.shields.io/github/last-commit/AlbertoCuadra/combustion_toolbox)
[![license](https://img.shields.io/github/license/AlbertoCuadra/combustion_toolbox)](https://www.gnu.org/licenses/gpl-3.0.html)

A MATLAB-GUI based open-source tool for solving gaseous combustion problems.

<!-- Website: https://combustiontoolbox.netlify.app/ -->
:top: There is also a (less complete) [Python version](https://github.com/AlbertoCuadra/Combustion-PyToolbox)


## Features

Combustion Toolbox is a a MATLAB-GUI based tool for solving gaseous combustion problems.
Features
  - The code stems from the minimization of the free energy of the system by using Lagrange multipliers combined with a Newton-Raphson method, upon condition that initial gas properties are defined by two functions of states (e.g., temperature and pressure)
  - When temperature is not externally imposed, the code retrieves a routine also based on Newton-Raphson method to find the equilibrium temperature
  - Solve processes that involve strong changes in the dynamic pressure, such as detonations and shock waves in the steady state
  - Find the equilibrium conditions of the different phenomena undergoing behind the shock: molecular vibrational excitation up to dissociation, and electronic excitation up to ionization, thereby providing the `properties of the gas in plasma state` within the temperature range given by the NASA’s 9-coefficient polynomial fits.
  - The corresponding thermodynamic properties of the species are modelled with `NASA’s 9-coefficient polynomial fits`, which ranges `up to 20000 K`, and the ideal gas equation of state
  - Results are in `excellent agreement with NASA’s Chemical Equilibrium with Applications (CEA) program`, CANTERA and Caltech’s Shock and Detonation Toolbox
  - All the routines and computations are encapsulated in a more comprehensive and user-friendly GUI
  - The code is in it’s transition to Python
  - Display predefined plots (e.g., molar fraction vs equilence ratio)
  - Export results in a spreadsheet (requires Excel)
  - Export results as a .mat format
* `Chemical equilibrium problems`
  - TP: Equilibrium composition at defined temperature and pressure
  - HP: Adiabatic temperature and composition at constant pressure
  - SP: Isentropic compression/expansion to a specified pressure
  - TV: Equilibrium composition at defined temperature and constant volume
  - EV: Adiabatic temperature and composition at constant volume
  - SV: Isentropic compression/expansion to a specified volume
* `Shock calculations:`
  - Pre-shock and post shock states
  - Equilibrium or frozen composition
  - Incident or reflected shocks
  - Chapman-Jouguet detonations and overdriven detonations
  - Reflected detonations
  - Oblique shocks/detonations
  - Shock polar solutions incident and reflected states, considering dissociation, ionization, and recombinantion in multi-species mixtures
  - Hugoniot curves
  - Ideal jump conditions for a given adiabatic index and pre-shock Mach number
* `Rocket propellant performance assuming:`
  - Infinite-Area-Chamber model (IAC)
  - Finite-Area-Chamber model (FAC) - under development -
* All the routines and computations are encapsulated in a more comprehensive and `user-friendly GUI`
* The code `is in it’s transition to Python`
* Export results in a spreadsheet
* Export results as a .mat format
* `Display predefined plots` (e.g., molar fraction vs equilence ratio)


<!--
---
⚠️ **NOTE**
- At the moment, the Python version does not have all the capabilities that the MATLAB version has. I will continue with the development of this version adding all the remaining capabilities. I will also add a GUI using Qt6 and Pyside6.
---
-->

This project is also part of the PhD of [Alberto Cuadra-Lara](https://www.acuadralara.com/).

---
⚠️ **NOTE**

- The first final version v1.0.0 is expected to be released in April 2022. Check out the  upcoming features of [Combustion Toolbox v1.0.0](https://github.com/AlbertoCuadra/combustion_toolbox/projects/2).
- [Wiki](https://github.com/AlbertoCuadra/combustion_toolbox/wiki) under construction :building_construction:

---
## Start here!
The [tutorial](https://github.com/AlbertoCuadra/combustion_toolbox/wiki/Tutorial) will help you get started using Combustion Toolbox on your pc.

## Gallery
We have several examples of what Combustion Toolbox can do. Here we show a preview of the GUI and some results obtained from Combustion Toolbox.

<p align="left">
    <img src="https://github.com/AlbertoCuadra/combustion_toolbox/blob/master/GUI/Snapshots/snapshot_1.svg" width="500">
</p>

**Figure 1:** *Snapshot of the GUI*.

<p align="left">
    <img src="https://github.com/AlbertoCuadra/combustion_toolbox/blob/master/Validations/Figures/Hugoniot_benchmarking.svg" width="400">
</p>
    
**Figure 2:** *Hugoniot curves for different molecular gases at pre-shock temperature T1 = 300 K and pressure p1 = 1 atm \[numerical results obtained with Combustion Toolbox (lines) and contrasted with NASA’s Chemical Equilibrium with Applications (CEA) code excluding ionization (symbols)\]*.
    
<p align="left">
    <img src="https://github.com/AlbertoCuadra/combustion_toolbox/blob/master/Validations/Figures/run_validation_DET_CEA_3_molar.svg" width="1200">
</p>

**Figure 3:** *Example CJ detonation for lean to rich CH4-air mixtures at standard conditions: (a) variation of molar fraction, (b) variation of temperature. The computational time was of 9.25 seconds using a Intel(R) Core(TM) i7-8700 CPU @ 3.20GHz for a set of 24 species considered and a total of 351 case studies.*

<p align="left">
    <img src="https://github.com/AlbertoCuadra/combustion_toolbox/blob/master/Validations/Figures/polar_shock_full_and_frozen_both_air_complete.svg" width="1000">
</p>

**Figure 4:** *Pressure-deflection shock polar (left) and wave angle-deflection shock polar (right) for an air mixture (78.084% N2, 20.9476% O2, 0.9365% Ar, 0.0319% CO2) at pre-shock temperature T1 = 300 K and pressure p1 = 1 atm, and a range of preshock Mach numbers M1 = [2, 14]; line: considering dissociation, ionization, and recombination in multi-species mixtures; dashed: considering a thermochemically frozen air mixture.*

## Contributing

Please read [CONTRIBUTING.md](https://github.com/AlbertoCuadra/combustion_toolbox/blob/master/CONTRIBUTING.md) for details on our code of conduct, and the process for submitting pull requests to us.

Please send feedback or inquiries to [acuadra@ing.uc3m.es](mailto:acuadra@ing.uc3m.es)

Thank you for testing Combustion Toolbox!

## Acknowledgements
* Combustion Toolbox's color palette is obtained from the following repository: Stephen (2021). ColorBrewer: Attractive and Distinctive Colormaps (https://github.com/DrosteEffect/BrewerMap), GitHub. Retrieved December 3, 2021.
* For validations, Combustion Toolbox uses CPU Info from the following repository: Ben Tordoff (2022). CPU Info (https://github.com/BJTor/CPUInfo/releases/tag/v1.3), GitHub. Retrieved March 22, 2022.
## Developers

* **[Alberto Cuadra-Lara](https://acuadralara.com/)** - *Core Developer and App designer*
* **César Huete** - *Developer*
* **Marcos Vera** - *Developer*

Grupo de Mecánica de Fluidos, Universidad Carlos III, Av. Universidad 30, 28911, Leganés, Spain

See also the list of [contributors](https://github.com/AlbertoCuadra/combustion_toolbox/blob/master/CONTRIBUTORS.md) who participated in this project.

## Citing Combustion Toolbox

```bibtex
@misc{combustiontoolbox,
    author = "Cuadra, A and Huete, C and Vera, M",
    title = "Combustion Toolbox: A MATLAB-GUI based open-source tool for solving combustion problems",
    year = 2022,
    note = "Version 0.8.9",
    doi = {https://doi.org/10.5281/zenodo.6369301}
}
```
