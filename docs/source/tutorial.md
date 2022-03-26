
# Tutorial

## Getting Started 
Start MATLAB and browse for folder where you have downloaded Combustion Toolbox. To include files in PATH run this command in the command window: 
```matlab
>> addpath(genpath(pwd))
```
First, using Combustion Toolbox, you have to initialize the tool (load databases, set default variables, ...). To do that at the prompt type:
```matlab
>> self = App
```
If files contained in Combustion Toolbox are correctly declared, you should see something like this:
```matlab
self = 

  struct with fields:

            E: [1×1 struct]
            S: [1×1 struct]
            C: [1×1 struct]
         Misc: [1×1 struct]
           PD: [1×1 struct]
           PS: [1×1 struct]
           TN: [1×1 struct]
    DB_master: [1×1 struct]
           DB: [1×1 struct]
```
## Setting the state
Indicate temperature and pressure of the initial mixture
 ```matlab
>> self = set_prop(self, 'TR', 300, 'pR', 1 * 1.01325);
```
Indicate species and nº moles of each species in the initial mixture
#### Individual case
For example, for a stochiometric CH4-ideal_air mixture:
```matlab
>> self.PD.S_Fuel     = {'CH4'};
>> self.PD.N_Fuel     = 1;
>> self.PD.S_Oxidizer = {'O2'};
>> self.PD.N_Oxidizer = 2;
>> self.PD.S_Inert    = {'N2'};
>> self.PD.N_Inert    = 7.52;
```
This is the same as specifying a unit value for the equivalence ratio:
```matlab
>> self.PD.S_Fuel     = {'CH4'};
>> self.PD.S_Oxidizer = {'O2'};
>> self.PD.S_Inert    = {'N2'};
>> self = set_prop(self, 'phi', 1);
>> self.PD.proportion_inerts_O2 = 79/21;
```
The last two lines of code establish the equivalence relation and the proportion of the inert species over the oxidazer, respectively. The number of moles is calculated considering that the number of moles of fuel is one.
#### Several cases
CT also allows the computation of a range of values of different properties. For example, in case we want to compute a range of values of the equivalence ratio, e.g., phi = 0.5:0.01:5, do this:
```matlab
>> self.PD.S_Fuel     = {'CH4'};
>> self.PD.S_Oxidizer = {'O2'};
>> self.PD.S_Inert    = {'N2'};
>> self = set_prop(self, 'phi', 0.5:0.01:5);
>> self.PD.proportion_inerts_O2 = 79/21;
```
## Chemical Equilibrium

Depending on the problem you want to solve, you may need to configure additional inputs. For example, to compute the equilibrium composition at a defined temperature and pressure (TP) we have to set these values as
```matlab
>> self = set_prop(self, 'pP', self.PD.pR.value, 'TP', 3000);
```
and to solve the aforementioned problem, run
```matlab
>> self = SolveProblem(self, 'TP');
```
The results are contained in self.PS. By default, this routine print the results through the command window (default: self.Misc.FLAG_RESULTS=true) which gives for the stoichiometric case (phi=1):
```
***********************************************************
-----------------------------------------------------------
Problem type: TP  | phi = 1.000
-----------------------------------------------------------
              |   REACTANTS     |       PRODUCTS
T [K]         |	    300.000	|	3000.000
p [bar]       |	      1.013	|	   1.013
r [kg/m3]     |	      1.123	|	   0.103
h [kJ/kg]     |	   -254.958	|	2573.851
e [kJ/kg]     |	   -345.224	|	1589.064
s [kJ/(kg-K)] |	     24.865	|	  37.636
cp [kJ/(kg-K)]|	      1.080	|	   5.562
gamma [-]     |	      1.386	|	   1.168
-----------------------------------------------------------
REACTANTS        Xi [-]
N2           	 7.1483e-01
O2           	 1.9011e-01
CH4          	 9.5057e-02
MINORS[+22]      0.0000e+00

TOTAL  		 1.0000e+00
-----------------------------------------------------------
PRODUCTS	 Xi [-]
N2           	 6.4761e-01
H2O          	 1.1166e-01
CO           	 5.8504e-02
OH           	 3.5944e-02
H2           	 3.0832e-02
CO2          	 2.8626e-02
H            	 2.7588e-02
O2           	 2.5917e-02
O            	 1.8097e-02
NO           	 1.5206e-02
N            	 1.1122e-05
HO2          	 9.5922e-06
NO2          	 3.1370e-06
N2O          	 7.8817e-07
HCO          	 1.2470e-07
NH2          	 1.1863e-07
NH3          	 3.0276e-08
HCN          	 4.4118e-09
CN           	 4.0118e-10
CH           	 5.7134e-13
CH3          	 1.5607e-13
CH4          	 1.1368e-14
MINORS[+3]       1.3720e-15

TOTAL  		 1.0000e+00
-----------------------------------------------------------
***********************************************************
```
## Postprocessed results (predefined plots for several cases)
There are some predefined charts based on the selected problem, in case you have calculated multiple cases. Just calling the routine
```matlab
>> postResults(self);
```
will reproduce **Figure 1** which represents the variation of the molar fraction with the equivalence ratio for lean to rich CH4-ideal_air mixtures at 3000 [K] and 1.01325 [bar]. 

<p align="left">
    <img src="https://github.com/AlbertoCuadra/combustion_toolbox/blob/master/Validations/Tutorial_1.svg" width="1000">
</p>

**Figure 1:** *Example TP: variation of molar fraction for lean to rich CH4-ideal_air mixtures at 3000 [K] and 1.01325 [bar], a set of 24 species considered and a total of 451 case studies.*

## Thermodynamic Properties

## Congratulations!
Congratulations you have finished the Combustion Toolbox Matlab tutorial! You should now be ready to begin using Combustion Toolbox on your own.