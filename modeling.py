import scipy.integrate
import sympy as sp
import json
import numpy as np
import scipy
from tqdm import tqdm
import matplotlib.pyplot as plt
from time import time

"""
MAKE SURE THAT FOR THE JSON:
    1. THE SYMBOLS USED IN THE EQUATION ARE ALL IN THE PARAMETERS
    2. THE EQUATION IS FORMATED CORRECTLY
    3. VARIABLES LIKE HILL CONSTANT SHOULDN'T BE NULL OR ELSE THEY WILL BE SET TO 100 IN THE CODE
"""

# default value if anything is NULL
DEFAULT_VALUE = 100

with open("entire.json") as f:
    data = json.load(f)


tracking_variables = ["ATP", "ADP", "NADH", "NADPH"]

for chunk in data:
    substrate, product = chunk["substrate"], chunk["product"]
    tracking_variables += [substrate, product]

tracking_variables = list(set(tracking_variables))


print(tracking_variables)
intial_values = [DEFAULT_VALUE for _ in range(len(tracking_variables))]


enzyme_info = {}
for chunk in data:
    enzyme = chunk["enzyme"]
    equation = chunk["equation"]
    parameters = chunk["parameters"]
    substrate = chunk["substrate"]
    product = chunk["product"]

    parameters_names = list(parameters.keys())
    local_dict = {name: sp.symbols(name) for name in parameters_names}
    expr = sp.sympify(equation, locals=local_dict)

    for symbol in expr.free_symbols: 
        if str(symbol) not in parameters_names:
            raise Exception(f"Equation for enzyme {enzyme} uses parameters that aren't specified")

    # for example maps D-Glucose -> A 
    tracked_variable_to_symbol = {}

    # these are all the symbols that are used in the equation
    for symbol in expr.free_symbols:
        symbol = str(symbol)

        value_symbol = parameters[symbol]["start_val"]
        type_symbol = parameters[symbol]["type"]
        species_symbol = parameters[symbol]["species"]

        if value_symbol == "NULL" or type_symbol == "concentration":
            value_symbol = DEFAULT_VALUE
        else:
            value_symbol = float(sp.sympify(value_symbol))


        # not a constant since we are tracking that over time
        if type_symbol == "concentration" and species_symbol in tracking_variables:
            tracked_variable_to_symbol[species_symbol] = symbol
            intial_values[tracking_variables.index(species_symbol)] = value_symbol
        # is a constant (i.e. Vmax, Km) so we substitue into the expression
        else:            
            expr = expr.subs({symbol: value_symbol})
    

    enzyme_info[enzyme] = (expr, tracked_variable_to_symbol, substrate, product)



def modeling(t, concentrations):
    changes = [0 for _ in range(len(concentrations))]

    for expr, tracked_variable_to_symbol, substrate, product in enzyme_info.values():
        substrate_idx, product_idx = tracking_variables.index(substrate), tracking_variables.index(product)

        symbol_to_concentration = {}
        
        for idx, concentration in enumerate(concentrations):
            tracked_variable = tracking_variables[idx]
            if tracked_variable in tracked_variable_to_symbol.keys():
                symbol_to_concentration[tracked_variable_to_symbol[tracked_variable]] = concentration

        rxn_rate = expr.evalf(subs=symbol_to_concentration)
        
        # conservation of matter (rxn_rate can't be greater than the existing concentration of substrate)
        rxn_rate = min(rxn_rate, concentrations[substrate_idx])

        changes[substrate_idx] -= rxn_rate
        changes[product_idx] += rxn_rate
    
    """
    # perturbations
    if (int(t) % 20) == 0:
        glucose_idx = tracking_variables.index("D-Glucose")
        changes[glucose_idx] = 100 - concentrations[glucose_idx]
    """

    return np.array(changes)


print("Parsing Done!")


NUM_ITER = 50
t = np.linspace(0, NUM_ITER, NUM_ITER * 10)


begin = time()
sol = scipy.integrate.solve_ivp(modeling, [0, NUM_ITER], intial_values, method="LSODA", dense_output=True)
z = sol.sol(t)
end = time()

print(f"Solving Done: seconds per iter: {(end-begin)/NUM_ITER}")

for idx, metabolite in enumerate(tracking_variables):
    datapoints = z[idx]
    
    # Plot the iterative data with markers and a solid line
    plt.plot(t, datapoints, color="purple", linestyle='-', linewidth=1.5)
    
    # Add title and axis labels with LaTeX for µM symbol
    plt.title(f"Comparison for {metabolite}", fontsize=16, fontweight='bold')
    plt.xlabel("Iteration", fontsize=14)
    plt.ylabel(r'Concentration ($\mu$M)', fontsize=14)

    # Enable grid for better readability of the plot
    plt.grid(True, linestyle='--', alpha=0.7)

    # Improve the layout to avoid overlap
    plt.tight_layout()

    # Save the plot with higher resolution (300 dpi)
    plt.savefig(f"graphs/{metabolite}", dpi=600)
    
    # Clear figure for the next plot
    plt.clf()


palmitoyl_coa = z[tracking_variables.index("Palmitoyl-CoA")]
np.save("palmitoyl-coa.npy", palmitoyl_coa)
np.save("time.npy", t)