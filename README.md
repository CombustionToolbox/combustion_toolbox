<img alt="Combustion Toolbox" align="left" href="https://combustion-toolbox-website.readthedocs.io" style="border-width:0" src="https://github.com/CombustionToolbox/combustion_toolbox/blob/master/gui/assets/logo_CT_noversion_matlab.png" width="115"/>

## Combustion Toolbox: A MATLAB-GUI based open-source tool for solving gaseous combustion problems

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5554911.svg)](https://doi.org/10.5281/zenodo.5554911)
[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=CombustionToolbox/combustion_toolbox&file=CONTENTS.m)
[![File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://es.mathworks.com/matlabcentral/fileexchange/101088-combustion-toolbox)
[![CD](https://github.com/CombustionToolbox/combustion_toolbox/actions/workflows/CD.yml/badge.svg)](https://github.com/CombustionToolbox/combustion_toolbox/actions/workflows/CD.yml)
[![Documentation](https://readthedocs.org/projects/combustion-toolbox-website/badge/?version=latest)](https://combustion-toolbox-website.readthedocs.io/en/latest/?badge=latest)
[![license](https://img.shields.io/github/license/CombustionToolbox/combustion_toolbox)](https://www.gnu.org/licenses/gpl-3.0.html)

<br>

<p align=center>
    <img src="https://github.com/CombustionToolbox/combustion_toolbox_website/blob/main/docs/source/_static/gif/example_det_overdriven_gui.gif" width="345">
    &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
    <img src="https://github.com/CombustionToolbox/combustion_toolbox_website/blob/main/docs/source/_static/gif/example_det_overdriven.gif" width="360">
</p>

## Main features

The Combustion Toolbox is a MATLAB-based thermochemical framework designed to solve problems involving chemical equilibrium for both gas- and condensed-phase species. The toolbox is composed of several modules, each of which is designed to solve a specific class of problems:

- **CT-EQUIL** computes the equilibrium composition of multi-component gas mixtures undergoing thermochemical transformations. The final equilibrium state is determined by a predefined set of chemical species (gaseous—including ions—or condensed phases) and two thermodynamic state functions, such as enthalpy and pressure, e.g., for isobaric combustion processes.
- **CT-SD** solves steady-state shock and detonation wave problems for both normal and oblique incidence.  
- **CT-ROCKET** estimates the theoretical performance of rocket engines under highly idealized conditions.
- **CT-TURBULENCE** performs detailed analyses of turbulent flows, including turbulent statistics computations, Helmholtz decomposition, and spectral analyses.

The framework also features a **user-friendly graphical user interface (GUI)**.

> For installation instructions and usage guidelines, visit the [Combustion Toolbox website](https://combustion-toolbox-website.readthedocs.io).


## Citing Combustion Toolbox

If you use the Combustion Toolbox in a publication, please cite it using the following references:

* *Cuadra, A., Huete, C., & Vera, M. (2024). Combustion Toolbox: An open-source thermochemical code for gas- and condensed-phase problems involving chemical equilibrium. [arXiv:2409.15086](https://doi.org/10.48550/arXiv.2409.15086).*
* *Cuadra, A., Huete, C., & Vera, M. (2024). Combustion Toolbox: A MATLAB-GUI based open-source tool for solving gaseous combustion problems. Version 1.2.0. Zenodo. [doi:10.5281/zenodo.5554911](https://doi.org/10.5281/zenodo.5554911).*

It can be handy the BibTeX format:

```bibtex
@article{cuadra2024a_preprint,
    title         = {{Combustion Toolbox: An open-source thermochemical code for gas- and condensed-phase problems involving chemical equilibrium}},
    author        = {Cuadra, A. and Huete, C. and Vera, M.},
    journal       = {{arXiv preprint arXiv:2409.15086}},
    year          = {2024},
    eprint        = {2409.15086},
    archivePrefix = {arXiv},
    primaryClass  = {physics.chem-ph},
    doi           = {10.48550/arXiv.2409.15086}
}

@misc{combustiontoolbox,
    author  = "Cuadra, A. and Huete, C. and Vera, M.",
    title   = "{Combustion Toolbox: A MATLAB-GUI based open-source tool for solving gaseous combustion problems}",
    year    = 2024,
    note    = "Version 1.2.0",
    doi     = {https://doi.org/10.5281/zenodo.5554911}
}
```
