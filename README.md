<img alt="Combustion Toolbox" align="left" href="https://combustion-toolbox-website.readthedocs.io" style="border-width:0" src="https://github.com/CombustionToolbox/combustion_toolbox/blob/master/gui/assets/logo_CT_noversion_matlab.png" width="115"/>

## Combustion Toolbox: A MATLAB-based framework for solving combustion and high-speed flow problems

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

- **CT-EQUIL** computes the thermochemical equilibrium composition of multi-component gas mixtures for prescribed set of chemical species (gaseous—including ions—or condensed phases) and thermodynamic constraints (e.g. pressure-enthalpy).
- **CT-SD** solves steady-state shock and detonation wave problems for both normal and oblique incidence.  
- **CT-ROCKET** estimates the theoretical performance of rocket engines under highly idealized conditions.
- **CT-LIA** predicts shock-turbulence interaction statistics using linear theory, accounting for thermochemical effects.
- **CT-TURBULENCE** performs spectral and statistical analysis of turbulent flows, including energy spectra, Helmholtz decomposition, and turbulence diagnostics.

The framework also includes an intuitive **graphical user interface (GUI)**, with a **royalty-free standalone version** available for Windows, macOS, and Linux.

> For installation instructions and usage guidelines, visit the [Combustion Toolbox website](https://combustion-toolbox-website.readthedocs.io).


## Citing Combustion Toolbox

If you use the Combustion Toolbox in a publication, please cite it using the following references:

* *Cuadra, A., Huete, C., & Vera, M. (2026). Combustion Toolbox: An open-source thermochemical code for gas-and condensed-phase problems involving chemical equilibrium. Computer Physics Communications 320, 110004. [doi:10.1016/j.cpc.2025.110004.](https://doi.org/10.1016/j.cpc.2025.110004).*
* *Cuadra, A., Huete, C., & Vera, M. (2026). Combustion Toolbox: A MATLAB-GUI based open-source tool for solving gaseous combustion problems. Version 1.2.8. Zenodo. [doi:10.5281/zenodo.5554911](https://doi.org/10.5281/zenodo.5554911).*

It can be handy the BibTeX format:

```bibtex
@article{cuadra2026a,
    title   = {{Combustion Toolbox: An open-source thermochemical code for gas- and condensed-phase problems involving chemical equilibrium}},
    author  = {A. Cuadra and C. Huete and M. Vera},
    journal = {Computer Physics Communications},
    volume  = {320},
    pages   = {110004},
    year    = {2026},
    issn    = {0010-4655},
    doi     = {https://doi.org/10.1016/j.cpc.2025.110004},
}

@misc{combustiontoolbox,
    title   = "{{Combustion Toolbox: A MATLAB-GUI based open-source tool for solving gaseous combustion problems}}",
    author  = {A. Cuadra and C. Huete and M. Vera},
    year    = {2026},
    note    = {Version 1.2.8},
    doi     = {https://doi.org/10.5281/zenodo.5554911}
}
```
