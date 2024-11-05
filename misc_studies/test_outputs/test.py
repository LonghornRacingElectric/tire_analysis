from tire_toolkit.tire_model.MF52 import MF52
import numpy as np

Tire = MF52(tire_name="test", file_path="./tir_files/Round_8_Hoosier_R25B_16x7p5_10_on_7in_12psi_PAC02_UM2.tir")
loads = Tire.tire_eval(FZ = -700, alpha = 10 * np.pi / 180, kappa = 0, gamma = 0)
print([float(x) for x in loads])