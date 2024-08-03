from reaktoro import *

# define this list of species 
gasList = [
    "O2", "N2", "H2O", "CH4", "CO", "CO2", "H2", "OH"
]

condensedList = ["CaCO3(cd)", "CaO(cd)", "Fe(cd)", "Fe3O4(cd)", "C(gr)"]

displaySpecoes = [
    "O2", "N2", "H2O", "CH4", "CO", "CO2", "H2", "OH", "CaCO3(cd)", "CaO(cd)", "Fe(cd)", "Fe3O4(cd)", "C(gr)"
]

db = NasaDatabase("nasa-cea")

gases = GaseousPhase(gasList)
condensed = CondensedPhases(condensedList)

condensed = CondensedPhases(speciate("Ca C H O Fe"))
        
system = ChemicalSystem(db, gases, condensed)

print("Condensed Phases")
for species in system.species():
    if species.aggregateState() == AggregateState.CondensedPhase:
        print(":: " + species.name())
        
state = ChemicalState(system)
state.temperature(1050, "kelvin")
state.pressure(1, "atm")
state.set("O2", 20.46, "mol")
state.set("N2", 187.19, "mol")
state.set("H2O", 1.775, "mol")
state.set("CH4", 2.2554, "mol")
state.set("Fe3O4(cd)", 13.1, "mol")
state.set("Fe(cd)", 3.527, "mol")
state.set("C(gr)", 85.59, "mol")
state.set("CaCO3(cd)", 0.1499, "mol")
state.set("CaO(cd)", 0.6063, "mol")


specs = EquilibriumSpecs(system)
specs.temperature()
specs.pressure()

options = EquilibriumOptions()
options.epsilon = 1e-30

solver = EquilibriumSolver(specs)
solver.setOptions(options)

solver.solve(state)  # equilibrate the `state` object!

print("=== FINAL STATE ===")
print(state)