from tire_toolkit import workflows
import yaml


with open("./runs.yml", "r") as file:
    runs = yaml.safe_load(file)

# Initialize tires
tires = workflows.import_tires(
    tire_names = runs["Tires"]["Tire Names"],
    tir_paths = runs["Tires"][".tir File Paths"]
    )

# Solve for updated tire scaling coefficients
if runs["Scale Virtual Tires"]["Run"]:
    scaling_coeffs = tires.coeff_matching(
        reference_tire = runs["Scale Virtual Tires"]["Reference"],
        target_tire = runs["Scale Virtual Tires"]["Target"],
        FZ_min = runs["Scale Virtual Tires"]["Minimum FZ"],
        FZ_max = runs["Scale Virtual Tires"]["Maximum FZ"],
        weighting = runs["Scale Virtual Tires"]["Weighting"],
        mesh = runs["Scale Virtual Tires"]["Refinement"])
    
    for tire_name in runs["Scale Virtual Tires"]["Tires to Scale"]:
        tires.set_scaling_coeffs(tire_name, scaling_coeffs)

# Compare two tires numerically
if runs["Tire Assessment"]["Run"]:
    print("Tire Assessment isn't Implemented")

# Compare tires graphically
if runs["First Principles Plots"]["Run"]:
    workflows.plot_tires(
        tire_object = tires, 
        names = runs["First Principles Plots"]["Tire Names"],
        FZ_min = runs["First Principles Plots"]["Minimum FZ"], 
        FZ_max = runs["First Principles Plots"]["Maximum FZ"], 
        file_name = runs["First Principles Plots"]["Output File Name"])