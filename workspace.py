import numpy as np
import matplotlib.pyplot as plt
from tire_toolkit import workflows
from tire_toolkit.assets.tire_model.MF52 import MF52

tires = workflows.import_tires(
    tire_names = 
        [
        "Hoosier_16_7.5x7_R25B_8psi",
        "Hoosier_16_7.5x7_R25B_10psi",
        "Hoosier_16_7.5x7_R25B_12psi",
        "Hoosier_16_7.5x7_R25B_14psi"
        ],
    tir_paths = 
        [
        "./tire_toolkit/assets/tir_files/Round_8_Hoosier_R25B_16x7p5_10_on_7in_8psi_PAC02_UM2.tir",
        "./tire_toolkit/assets/tir_files/Round_8_Hoosier_R25B_16x7p5_10_on_7in_10psi_PAC02_UM2.tir",
        "./tire_toolkit/assets/tir_files/Round_8_Hoosier_R25B_16x7p5_10_on_7in_12psi_PAC02_UM2.tir",
        "./tire_toolkit/assets/tir_files/Round_8_Hoosier_R25B_16x7p5_10_on_7in_14psi_PAC02_UM2.tir"
        ]
)

workflows.plot_tires(
    tire_object = tires, 
    names = [
        "Hoosier_16_7.5x7_R25B_8psi",
        "Hoosier_16_7.5x7_R25B_10psi",
        "Hoosier_16_7.5x7_R25B_12psi",
        "Hoosier_16_7.5x7_R25B_14psi"], 
    FZ_min = 100, 
    FZ_max = 1091, 
    file_name = "First_Principles_Metrics.pdf")

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