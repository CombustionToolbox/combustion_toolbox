# Thermodynamics properties
Functions to obtain thermodynamic properties from a given mixture:

***
## adiabaticIndex(mix)

Get the adiabatic index [-].
* **Input:** (mix)
  * instance of struct mixture. 
* **Output:** [val]
  * adiabatic index of the mixture [-].
***
## cp_mass(mix)

Get the mass-basis specific heat at constant pressure [kJ/kg-K].
* **Input:** (mix)
  * instance of struct mixture. 
* **Output:** [val]
  * mass basis specific heat of the mixture at constant pressure [kJ/kg-K].
***
## cp_mole(mix)

Get the mole-basis specific heat at constant pressure [kJ/kmol-K].
* **Input:** (mix)
  * instance of struct mixture. 
* **Output:** [val]
  * mole basis specific heat of the mixture at constant pressure [kJ/kmol-K].
***
## cv_mass(mix)

Get the mass-basis specific heat at constant volume [kJ/kg-K].
* **Input:** (mix)
  * instance of struct mixture. 
* **Output:** [val]
  * mass basis specific heat of the mixture at constant volume [kJ/kg-K].
***
## cv_mole(mix)

Get the mole-basis specific heat at constant volume [kJ/kmol-K].
* **Input:** (mix)
  * instance of struct mixture. 
* **Output:** [val]
  * mole basis specific heat of the mixture at constant volume [kJ/kmol-K].
***
## density(mix)

Get the density [kg/m3].
* **Input:** (mix)
  * instance of struct mixture. 
* **Output:** [val]
  * density of the mixture [kg/m3].
***
## enthalpy_mass(mix)

Get the mass specific enthalpy [kJ/kg].
* **Input:** (mix)
  * instance of struct mixture. 
* **Output:** [val]
  * mass specific enthalpy of the mixture [kJ/kg].
***
## enthalpy_mole(mix)

Get the mole specific enthalpy [kJ/kmol].
* **Input:** (mix)
  * instance of struct mixture. 
* **Output:** [val]
  * mole specific enthalpy of the mixture [kJ/kmol].
***
## entropy_mass(mix)

Get the mass specific entropy [kJ/kg-K].
* **Input:** (mix)
  * instance of struct mixture. 
* **Output:** [val]
  * mass specific entropy of the mixture [kJ/kmol].
***
## entropy_mole(mix)

Get the mole specific entropy [kJ/kmol-K].
* **Input:** (mix)
  * instance of struct mixture. 
* **Output:** [val]
  * mole specific entropy of the mixture [kJ/kmol-K].
***
## equivalenceRatio(mix)

Get the equivalence ratio of the initial mixture [-].
* **Input:** (mix)
  * instance of struct mixture. 
* **Output:** [val]
  * equivalence ratio of the initial mixture [-].
***
## intEnergy_mass(mix)

Get the mass specific internal energy [kJ/kg].
* **Input:** (mix)
  * instance of struct mixture. 
* **Output:** [val]
  * mass specific internal energy of the mixture [kJ/kg].
***
## intEnergy_mole(mix)

Get the mole specific internal energy [kJ/kmol].
* **Input:** (mix)
  * instance of struct mixture. 
* **Output:** [val]
  * mole specific internal energy of the mixture [kJ/kmol].
***
## massFractions(mix)

Get the mass fractions of all the species in the mixture [-].
* **Input:** (mix)
  * instance of struct mixture. 
* **Output:** [val]
  * vector with the mass fractions of all the species in the mixture [kJ/kmol].
***
## meanMolecularWeight(mix)

Get the mean molecular weight [kg/kmol].
* **Input:** (mix)
  * instance of struct mixture. 
* **Output:** [val]
  * mean molecular weight of the mixture [kg/kmol].
***
## pressure(mix)

Get the pressure [bar].
* **Input:** (mix)
  * instance of struct mixture. 
* **Output:** [val]
  * pressure of the mixture [bar].
***
## temperature(mix)

Get the temperature [K].
* **Input:** (mix)
  * instance of struct mixture. 
* **Output:** [val]
  * temperature of the mixture [K].
***
## velocity_relative(mix)

Get the velocity of the gases relative to the shock front [m/s]. Only for shocks and detonations problems.
* **Input:** (mix)
  * instance of struct mixture.  
* **Output:** [val]
  * velocity of the gases relative to the shock front [m/s].
***