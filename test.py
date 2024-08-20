from tire_toolkit.assets.tire_model.MF52 import MF52
import numpy as np

Tire = MF52(tire_name="test", file_path="./tir_files/Modified_Round_8_Hoosier_R25B_16x7p5_10_on_7in_12psi_PAC02_UM2.tir")
print([float(x) for x in Tire.tire_eval(FZ=654, alpha=5 * np.pi / 180, kappa=0, gamma=-20 * np.pi / 180)])