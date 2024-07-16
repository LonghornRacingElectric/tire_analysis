from LHR_tire_toolkit.process_tire import ProcessTire
import numpy as np

Tire = ProcessTire(tire_name="test", tir_path="./tir_files/Round_8_Hoosier_R25B_16x7p5_10_on_7in_12psi_PAC02_UM2.tir")
print(Tire.get_loads(FZ = 700, alpha = 10 * np.pi / 180, kappa = 0, gamma = 0))