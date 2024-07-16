import os
# import git
import yaml
import numpy as np
import matplotlib.pyplot as plt
from tire_toolkit import workflows
from tire_toolkit.assets.tire_model.MF52 import MF52


# cwd = os.getcwd()

# git_repo = git.Repo(cwd, search_parent_directories = True)
# git_root = git_repo.git.rev_parse("--show-toplevel")

# os.chdir(git_root)

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

# new_tire = MF52("Hoosier_16_7.5x7_R25B_12psi", "./tire_toolkit/assets/tir_files/Round_8_Hoosier_R25B_16x7p5_10_on_7in_12psi_PAC02_UM2.tir")

# SA_sweep = np.linspace(-20 * np.pi / 180, 20 * np.pi / 180, 100)
# FZ_sweep = np.linspace(100, 4000, 100)

# SA_vals = []

# print(new_tire.tire_eval(FZ = 1091, alpha = 15 * np.pi / 180, kappa = 0, gamma = 0)[1])

# for load in FZ_sweep:
#     SA_vals.append(max(new_tire.tire_eval(FZ = load, alpha = SA_sweep, kappa = 0, gamma = 0)[1]))

# plt.figure()

# plt.plot(FZ_sweep, SA_vals)
# plt.title("Peak FY vs FZ")
# plt.xlabel("Normal Load (N)")
# plt.ylabel("Peak Lateral Force (N)")

# plt.savefig("Peak_FY_vs_FZ.png")

# plt.show()