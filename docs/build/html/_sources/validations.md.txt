# Validations

A set of the results obtained using Combustion Toolbox, [CEA-NASA](https://cearun.grc.nasa.gov/), [CANTERA](https://cantera.org/) & [SD-Toolbox](https://shepherd.caltech.edu/EDL/PublicResources/sdt/)<sup>1</sup>, and [TEA](https://github.com/dzesmin/TEA).

<sup>1</sup> The Shock & Detonation Toolbox uses the Cantera software package as kernel for the thermochemical calculations.

***
**_For the sake of clarity, we only show a reduced set of species in the validation of the mole fractions._**
To run all the validations contrasted with CEA at once, at the prompt type:
```matlab
>> run_validations_CEA
```
## Validation TP 1

* Contrasted with: NASA's Chemical Equilibrium with Applications software
* Problem type: Equilibrium composition at defined T and p
* Temperature [K]   = 2500
* Pressure    [bar] = 1
* Initial mixture [moles]:

  <img src="https://latex.codecogs.com/svg.image?\text{C}_6\text{H}_6&space;&plus;&space;\frac{7.5}{\phi}(\text{O}_2&space;&plus;&space;3.76\text{N}_2),&space;\text{with&space;a&space;equivalence&space;ratio&space;}&space;\phi&space;\in&space;[0.5,&space;4]&space;" title="\text{C}_6\text{H}_6 + \frac{7.5}{\phi}(\text{O}_2 + 3.76\text{N}_2),&space;\text{with&space;a&space;equivalence&space;ratio&space;}&space;\phi&space;\in&space;[0.5,&space;4]&space;" />

* List of species considered = ListSpecies('Soot formation Extended')
* URL Folder Results CEA: https://github.com/AlbertoCuadra/combustion_toolbox/tree/master/Validations/CEA/Data/TP

To repeat the results, run:
```matlab
>> run_validation_TP_CEA_1.m
```

<p align="center">
    <img src="https://github.com/AlbertoCuadra/combustion_toolbox/blob/master/Validations/Figures/run_validation_TP_CEA_1_molar.svg" width="1000">
    <img src="https://github.com/AlbertoCuadra/combustion_toolbox/blob/master/Validations/Figures/run_validation_TP_CEA_1_properties.svg" width="600">
</p>

## Validation TP 2

* Contrasted with: NASA's Chemical Equilibrium with Applications software
* Problem type: Equilibrium composition at defined T and p
* Temperature [K]   = 2500
* Pressure    [bar] = 1
* Initial mixture [moles]:

  <img src="https://latex.codecogs.com/svg.image?\text{C}_6\text{H}_6&space;&plus;&space;\frac{7.5}{\phi}\text{O}_2,&space;\text{with&space;a&space;equivalence&space;ratio&space;}&space;\phi&space;\in&space;[0.5,&space;4]&space;" title="\text{C}_6\text{H}_6 + \frac{7.5}{\phi}\text{O}_2,&space;\text{with&space;a&space;equivalence&space;ratio&space;}&space;\phi&space;\in&space;[0.5,&space;4]&space;" />
* List of species considered = ListSpecies('Soot formation Extended')
* URL Folder Results CEA: https://github.com/AlbertoCuadra/combustion_toolbox/tree/master/Validations/CEA/Data/TP

To repeat the results, run:
```matlab
>> run_validation_TP_CEA_2.m
```

<p align="center">
    <img src="https://github.com/AlbertoCuadra/combustion_toolbox/blob/master/Validations/Figures/run_validation_TP_CEA_2_molar.svg" width="1000">
    <img src="https://github.com/AlbertoCuadra/combustion_toolbox/blob/master/Validations/Figures/run_validation_TP_CEA_2_properties.svg" width="600">
</p>

## Validation TP 3

* Contrasted with: NASA's Chemical Equilibrium with Applications software
* Problem type: Equilibrium composition at defined T and p
* Temperature [K]   = 2500
* Pressure    [bar] = 1
* Initial mixture [moles]:

  <img src="https://latex.codecogs.com/svg.image?\text{C}\text{H}_3\text{OH}&space;&plus;&space;\frac{1.5}{\phi}(\text{O}_2&space;&plus;&space;3.76\text{N}_2),&space;\text{with&space;a&space;equivalence&space;ratio&space;}&space;\phi&space;\in&space;[0.5,&space;4]&space;" title="\text{C}\text{H}_3\text{OH} + \frac{1.5}{\phi}(\text{O}_2 + 3.76\text{N}_2),&space;\text{with&space;a&space;equivalence&space;ratio&space;}&space;\phi&space;\in&space;[0.5,&space;4]&space;" />
* List of species considered = ListSpecies('Soot formation Extended')
* URL Folder Results CEA: https://github.com/AlbertoCuadra/combustion_toolbox/tree/master/Validations/CEA/Data/TP

To repeat the results, run:
```matlab
>> run_validation_TP_CEA_3.m
```

<p align="center">
    <img src="https://github.com/AlbertoCuadra/combustion_toolbox/blob/master/Validations/Figures/run_validation_TP_CEA_3_molar.svg" width="1000">
    <img src="https://github.com/AlbertoCuadra/combustion_toolbox/blob/master/Validations/Figures/run_validation_TP_CEA_3_properties.svg" width="600">
</p>

## Validation TP 4

* Contrasted with: NASA's Chemical Equilibrium with Applications software
* Problem type: Equilibrium composition at defined T and p
* Temperature [K]   = 2500
* Pressure    [bar] = 1
* Initial mixture [moles]:

  <img src="https://latex.codecogs.com/svg.image?\text{C}\text{H}_3\text{OH}&space;&plus;&space;\frac{1.5}{\phi}\text{O}_2,&space;\text{with&space;a&space;equivalence&space;ratio&space;}&space;\phi&space;\in&space;[0.5,&space;4]&space;" title="\text{C}\text{H}_3\text{OH} + \frac{1.5}{\phi}\text{O}_2,&space;\text{with&space;a&space;equivalence&space;ratio&space;}&space;\phi&space;\in&space;[0.5,&space;4]&space;" />
* List of species considered = ListSpecies('Soot formation Extended')
* URL Folder Results CEA: https://github.com/AlbertoCuadra/combustion_toolbox/tree/master/Validations/CEA/Data/TP

To repeat the results, run:
```matlab
>> run_validation_TP_CEA_4.m
```

<p align="center">
    <img src="https://github.com/AlbertoCuadra/combustion_toolbox/blob/master/Validations/Figures/run_validation_TP_CEA_4_molar.svg" width="1000">
    <img src="https://github.com/AlbertoCuadra/combustion_toolbox/blob/master/Validations/Figures/run_validation_TP_CEA_4_properties.svg" width="600">
</p>

## Validation HP 1

* Contrasted with: NASA's Chemical Equilibrium with Applications software
* Problem type: Adiabatic T and composition at constant p
* Temperature [K]   = 300
* Pressure    [bar] = 1
* Initial mixture [moles]:

  <img src="https://latex.codecogs.com/svg.image?\text{C}_2\text{H}_2\text{acetylene}&space;&plus;&space;\frac{2.5}{\phi}(\text{O}_2&space;&plus;&space;3.76\text{N}_2),&space;\text{with&space;a&space;equivalence&space;ratio&space;}&space;\phi&space;\in&space;[0.5,&space;4]&space;" title="\text{C}_2\text{H}_2\text{acetylene} + \frac{2.5}{\phi}(\text{O}_2 + 3.76\text{N}_2),&space;\text{with&space;a&space;equivalence&space;ratio&space;}&space;\phi&space;\in&space;[0.5,&space;4]&space;" />
* List of species considered = ListSpecies('Soot formation Extended')
* URL Folder Results CEA: https://github.com/AlbertoCuadra/combustion_toolbox/tree/master/Validations/CEA/Data/HP

To repeat the results, run:
```matlab
>> run_validation_HP_CEA_1.m
```

<p align="center">
    <img src="https://github.com/AlbertoCuadra/combustion_toolbox/blob/master/Validations/Figures/run_validation_HP_CEA_1_molar.svg" width="1000">
    <img src="https://github.com/AlbertoCuadra/combustion_toolbox/blob/master/Validations/Figures/run_validation_HP_CEA_1_properties.svg" width="600">
</p>

## Validation HP 2

* Contrasted with: NASA's Chemical Equilibrium with Applications software
* Problem type: Adiabatic T and composition at constant p
* Temperature [K]   = 300
* Pressure    [bar] = 1
* Initial mixture [moles]:

  <img src="https://latex.codecogs.com/svg.image?\text{C}_2\text{H}_2\text{acetylene}&space;&plus;&space;\frac{2.5}{\phi}\text{O}_2,&space;\text{with&space;a&space;equivalence&space;ratio&space;}&space;\phi&space;\in&space;[0.5,&space;4]&space;" title="\text{C}_2\text{H}_2\text{acetylene} + \frac{2.5}{\phi}\text{O}_2,&space;\text{with&space;a&space;equivalence&space;ratio&space;}&space;\phi&space;\in&space;[0.5,&space;4]&space;" />
* List of species considered = ListSpecies('Soot formation Extended')
* URL Folder Results CEA: https://github.com/AlbertoCuadra/combustion_toolbox/tree/master/Validations/CEA/Data/HP

To repeat the results, run:
```matlab
>> run_validation_HP_CEA_2.m
```

<p align="center">
    <img src="https://github.com/AlbertoCuadra/combustion_toolbox/blob/master/Validations/Figures/run_validation_HP_CEA_2_molar.svg" width="1000">
    <img src="https://github.com/AlbertoCuadra/combustion_toolbox/blob/master/Validations/Figures/run_validation_HP_CEA_2_properties.svg" width="600">
</p>

## Validation HP 3

* Contrasted with: NASA's Chemical Equilibrium with Applications software
* Problem type: Adiabatic T and composition at constant p
* Temperature [K]   = 300
* Pressure    [bar] = 1
* Initial mixture [moles]:

  <img src="https://latex.codecogs.com/svg.image?\text{C}\text{H}_4&space;&plus;&space;\frac{2}{\phi}(\text{O}_2&space;&plus;&space;3.76\text{N}_2),&space;\text{with&space;a&space;equivalence&space;ratio&space;}&space;\phi&space;\in&space;[0.5,&space;4]&space;" title="\text{C}\text{H}_4 + \frac{2}{\phi}(\text{O}_2 + 3.76\text{N}_2),&space;\text{with&space;a&space;equivalence&space;ratio&space;}&space;\phi&space;\in&space;[0.5,&space;4]&space;" />
* List of species considered = ListSpecies('Soot formation Extended')
* URL Folder Results CEA: https://github.com/AlbertoCuadra/combustion_toolbox/tree/master/Validations/CEA/Data/HP

To repeat the results, run:
```matlab
>> run_validation_HP_CEA_3.m
```

<p align="center">
    <img src="https://github.com/AlbertoCuadra/combustion_toolbox/blob/master/Validations/Figures/run_validation_HP_CEA_3_molar.svg" width="1000">
    <img src="https://github.com/AlbertoCuadra/combustion_toolbox/blob/master/Validations/Figures/run_validation_HP_CEA_3_properties.svg" width="600">
</p>

## Validation HP 4

* Contrasted with: NASA's Chemical Equilibrium with Applications software
* Problem type: Adiabatic T and composition at constant p
* Temperature [K]   = 300
* Pressure    [bar] = 1
* Initial mixture [moles]:

  <img src="https://latex.codecogs.com/svg.image?\text{C}\text{H}_4&space;&plus;&space;\frac{2}{\phi}\text{O}_2,&space;\text{with&space;a&space;equivalence&space;ratio&space;}&space;\phi&space;\in&space;[0.5,&space;4]&space;" title="\text{C}\text{H}_4 + \frac{2}{\phi}\text{O}_2,&space;\text{with&space;a&space;equivalence&space;ratio&space;}&space;\phi&space;\in&space;[0.5,&space;4]&space;" />
* List of species considered = ListSpecies('Soot formation Extended')
* URL Folder Results CEA: https://github.com/AlbertoCuadra/combustion_toolbox/tree/master/Validations/CEA/Data/HP

To repeat the results, run:
```matlab
>> run_validation_HP_CEA_4.m
```

<p align="center">
    <img src="https://github.com/AlbertoCuadra/combustion_toolbox/blob/master/Validations/Figures/run_validation_HP_CEA_4_molar.svg" width="1000">
    <img src="https://github.com/AlbertoCuadra/combustion_toolbox/blob/master/Validations/Figures/run_validation_HP_CEA_4_properties.svg" width="600">
</p>

## Validation DET 1

* Contrasted with: NASA's Chemical Equilibrium with Applications software
* Problem type: Chapman-Jouget Detonation
* Temperature [K]   = 300
* Pressure    [bar] = 1
* Initial mixture [moles]:

  <img src="https://latex.codecogs.com/svg.image?\text{C}_2\text{H}_2\text{acetylene}&space;&plus;&space;\frac{2.5}{\phi}(\text{O}_2&space;&plus;&space;3.76\text{N}_2),&space;\text{with&space;a&space;equivalence&space;ratio&space;}&space;\phi&space;\in&space;[0.5,&space;4]&space;" title="\text{C}_2\text{H}_2\text{acetylene} + \frac{2.5}{\phi}(\text{O}_2 + 3.76\text{N}_2),&space;\text{with&space;a&space;equivalence&space;ratio&space;}&space;\phi&space;\in&space;[0.5,&space;4]&space;" />
* List of species considered = ListSpecies('Soot formation Extended')
* URL Folder Results CEA: https://github.com/AlbertoCuadra/combustion_toolbox/tree/master/Validations/CEA/Data/HP

To repeat the results, run:
```matlab
>> run_validation_DET_CEA_1.m
```

<p align="center">
    <img src="https://github.com/AlbertoCuadra/combustion_toolbox/blob/master/Validations/Figures/run_validation_DET_CEA_1_molar.svg" width="1000">
</p>

## Validation DET 2

* Contrasted with: NASA's Chemical Equilibrium with Applications software
* Problem type: Chapman-Jouget Detonation
* Temperature [K]   = 300
* Pressure    [bar] = 1
* Initial mixture [moles]:

  <img src="https://latex.codecogs.com/svg.image?\text{C}_2\text{H}_2\text{acetylene}&space;&plus;&space;\frac{2.5}{\phi}\text{O}_2,&space;\text{with&space;a&space;equivalence&space;ratio&space;}&space;\phi&space;\in&space;[0.5,&space;4]&space;" title="\text{C}_2\text{H}_2\text{acetylene} + \frac{2.5}{\phi}\text{O}_2,&space;\text{with&space;a&space;equivalence&space;ratio&space;}&space;\phi&space;\in&space;[0.5,&space;4]&space;" />
* List of species considered = ListSpecies('Soot formation Extended')
* URL Folder Results CEA: https://github.com/AlbertoCuadra/combustion_toolbox/tree/master/Validations/CEA/Data/DET

To repeat the results, run:
```matlab
>> run_validation_DET_CEA_2.m
```

<p align="center">
    <img src="https://github.com/AlbertoCuadra/combustion_toolbox/blob/master/Validations/Figures/run_validation_DET_CEA_2_molar.svg" width="1000">
</p>

## Validation DET 3

* Contrasted with: NASA's Chemical Equilibrium with Applications software
* Problem type: Chapman-Jouget Detonation
* Temperature [K]   = 300
* Pressure    [bar] = 1
* Initial mixture [moles]:

  <img src="https://latex.codecogs.com/svg.image?\text{C}\text{H}_4&space;&plus;&space;\frac{2}{\phi}(\text{O}_2&space;&plus;&space;3.76\text{N}_2),&space;\text{with&space;a&space;equivalence&space;ratio&space;}&space;\phi&space;\in&space;[0.5,&space;4]&space;" title="\text{C}\text{H}_4 + \frac{2}{\phi}(\text{O}_2 + 3.76\text{N}_2),&space;\text{with&space;a&space;equivalence&space;ratio&space;}&space;\phi&space;\in&space;[0.5,&space;4]&space;" />
* List of species considered = ListSpecies('Soot formation Extended')
* URL Folder Results CEA: https://github.com/AlbertoCuadra/combustion_toolbox/tree/master/Validations/CEA/Data/DET

To repeat the results, run:
```matlab
>> run_validation_DET_CEA_3.m
```

<p align="center">
    <img src="https://github.com/AlbertoCuadra/combustion_toolbox/blob/master/Validations/Figures/run_validation_DET_CEA_3_molar.svg" width="1000">
</p>

## Validation DET 4

* Contrasted with: NASA's Chemical Equilibrium with Applications software
* Problem type: Chapman-Jouget Detonation
* Temperature [K]   = 300
* Pressure    [bar] = 1
* Initial mixture [moles]:

  <img src="https://latex.codecogs.com/svg.image?\text{C}\text{H}_4&space;&plus;&space;\frac{2}{\phi}\text{O}_2,&space;\text{with&space;a&space;equivalence&space;ratio&space;}&space;\phi&space;\in&space;[0.5,&space;4]&space;" title="\text{C}\text{H}_4 + \frac{2}{\phi}\text{O}_2,&space;\text{with&space;a&space;equivalence&space;ratio&space;}&space;\phi&space;\in&space;[0.5,&space;4]&space;" />
* List of species considered = ListSpecies('Soot formation Extended')
* URL Folder Results CEA: https://github.com/AlbertoCuadra/combustion_toolbox/tree/master/Validations/CEA/Data/DET

To repeat the results, run:
```matlab
>> run_validation_DET_CEA_4.m
```

<p align="center">
    <img src="https://github.com/AlbertoCuadra/combustion_toolbox/blob/master/Validations/Figures/run_validation_DET_CEA_4_molar.svg" width="1000">
</p>

## Validation SHOCK IONIZATION 1

* Contrasted with: NASA's Chemical Equilibrium with Applications software
* Problem type: Shock incident
* Temperature [K]   = 300
* Pressure    [bar] = 1
* Initial mixture [moles]:

  <img src="https://latex.codecogs.com/svg.image?3.7276&space;\text{N}_2&space;&plus;&space;\text{O}_2&space;&plus;&space;0.0447&space;\text{Ar}&space;&plus;&space;0.0015&space;\text{CO}_2"/>
* List of species considered = {'eminus', 'Ar', 'Arplus', 'C', 'Cplus', 'Cminus', 'CN', 'CNplus', 'CNminus', 'CNN', 'CO', 'COplus', 'CO2', 'CO2plus', 'C2', 'C2plus', 'C2minus', 'CCN', 'CNC', 'OCCN', 'C2N2', 'C2O', 'C3', 'N', 'Nplus', 'Nminus', 'NCO', 'NO', 'NOplus', 'NO2', 'NO2minus', 'NO3', 'NO3minus', 'N2', 'N2plus', 'N2minus', 'N2O', 'NCN', 'N2Oplus', 'N2O3', 'N3', 'O', 'Oplus', 'Ominus', 'O2', 'O2plus', 'O2minus', 'O3'}
* URL Folder Results CEA: https://github.com/AlbertoCuadra/combustion_toolbox/tree/master/Validations/CEA/Data/Shocks

To repeat the results, run:
```matlab
>> run_validation_SHOCK_IONIZATION_CEA_1.m
```

<p align="center">
    <img src="https://github.com/AlbertoCuadra/combustion_toolbox/blob/master/Validations/Figures/run_validation_SHOCK_IONIZATION_CEA_1_properties_2.svg" width="600">
    <img src="https://github.com/AlbertoCuadra/combustion_toolbox/blob/master/Validations/Figures/run_validation_SHOCK_IONIZATION_CEA_1_properties_3.svg" width="600">
</p>

## Validation SHOCK POLAR 1

* Contrasted with: Caltech's SD Toolbox and CANTERA
* Problem type: Shock Polar
* Temperature [K]   = 300
* Pressure    [bar] = 1.01325
* Initial mixture [moles]:

  <img src="https://latex.codecogs.com/svg.image?3.7619&space;\text{N}_2&space;&plus;&space;\text{O}_2"/>
* List of species considered = Frozen
* URL Folder Results SDToolbox: https://github.com/AlbertoCuadra/combustion_toolbox/tree/master/Validations/SDToolbox/Data

To repeat the results, run:
```matlab
>> run_validation_SHOCK_POLAR_SDToolbox_1.m
```

<p align="center">
    <img src="https://github.com/AlbertoCuadra/combustion_toolbox/blob/master/Validations/Figures/run_validation_SHOCK_POLAR_SDToolbox_1_properties_1.svg" width="400">
    <img src="https://github.com/AlbertoCuadra/combustion_toolbox/blob/master/Validations/Figures/run_validation_SHOCK_POLAR_SDToolbox_1_properties_2.svg" width="400">
</p>

## Validation TEA 1
* Contrasted with: Thermochemical Equilibrium Abundances of chemical species software
* Problem type: Equilibrium composition at defined T and p
* Temperature [K]   = linspace(500, 5000)
* Pressure    [bar] = logspace(-5, 2)
* Initial mixture [moles]:
  + H  = 1.0000000000e+00
  + He = 8.5113803820e-02
  + C  = 2.6915348039e-04
  + N  = 6.7608297539e-05
  + O  = 4.8977881937e-04
* List of species considered = {'C', 'CH4', 'CO2', 'CO', 'H2', 'H', 'H2O', 'He', 'N2', 'N', 'NH3', 'O'}
* URL Results TEA: https://github.com/dzesmin/TEA/blob/master/doc/examples/quick_example/results/quick_example.tea

To repeat the results, run:
```matlab
>> run_validation_TP_TEA.m
```

<p align="left">
    <img src="https://github.com/AlbertoCuadra/combustion_toolbox/blob/master/Validations/Figures/run_validation_TP_TEA_1_molar.svg" width="1000">
</p> 