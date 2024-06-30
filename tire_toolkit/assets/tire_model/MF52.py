import numpy as np
import matplotlib.pyplot as plt

from .file_processing._process_tir import _Processor

# from .MF52_calculations._longitudinal_force import
from .MF52_calculations._lateral_force import get_F_y
from .MF52_calculations._aligning_moment import get_M_z


import pandas as pd
from mpl_toolkits.mplot3d import Axes3D

class MF52:
    def __init__(self, tire_name: str, file_path: str) -> None:
        self.tire_name = tire_name
        self._tire_params = _Processor(name = tire_name, file_path = file_path)

        self._init_coeffs()
        self._init_consts()
        self._pure_slip_metrics = self._param_eval(self._vertical_coeffs["FNOMIN"])

    def _init_coeffs(self) -> None:
        self._lat_coeffs = self._tire_params.get_parameters("LATERAL_COEFFICIENTS")
        self._long_coeffs = self._tire_params.get_parameters("LONGITUDINAL_COEFFICIENTS")
        self._aligning_coeffs = self._tire_params.get_parameters("ALIGNING_COEFFICIENTS")

        self.scaling_coeffs = self._tire_params.get_parameters("SCALING_COEFFICIENTS")

    def _init_consts(self) -> None:
        self._dimensions = self._tire_params.get_parameters("DIMENSION")
        self._operating_conds = self._tire_params.get_parameters("TYRE_CONDITIONS")
        self._vertical_coeffs = self._tire_params.get_parameters("VERTICAL")
        self.structural = self._tire_params.get_parameters("STRUCTURAL")
    
    def _param_eval(self, FZ: float) -> None:
        self.mu_x = self.get_mu(FZ)[0]
        self.mu_y = self.get_mu(FZ)[1]
        self.C_y = self.get_cornering_stiffness(FZ, 0, 0.25)
        self.C_mz = self.get_aligning_stiffness(FZ, 0, 0.25)

    def tire_eval(self, FZ: float, alpha: float, kappa: float, gamma: float) -> list[float]:
        
        F_y_result = get_F_y(
            lat_coeffs = self._lat_coeffs,
            scaling_coeffs = self.scaling_coeffs,
            vertical_coeffs = self._vertical_coeffs,
            dimensions = self._dimensions,
            operating_conditions = self._operating_conds, 
            FZ = FZ,
            alpha = alpha,
            kappa = kappa,
            gamma = gamma
        )[1]

        MZ_result = get_M_z(
            aligning_coeffs = self._aligning_coeffs,
            scaling_coeffs = self.scaling_coeffs,
            lat_coeffs = self._lat_coeffs,
            long_coeffs = self._long_coeffs,
            vertical_coeffs = self._vertical_coeffs,
            dimensions = self._dimensions,
            operating_conditions = self._operating_conds,
            FZ = FZ,
            alpha = alpha,
            kappa = kappa,
            gamma = gamma
        )[1]

        return [0, F_y_result, FZ, 0, 0, MZ_result]

    def get_cornering_stiffness(self, FZ: float = 0, x_1: float = 0, x_2: float = 0.25, alpha: float = 0, kappa: float = 0, gamma: float = 0) -> float:
        if FZ == 0:
            FZ = self._vertical_coeffs["FNOMIN"]
        
        F_y_1 = self.tire_eval(FZ = FZ, kappa = kappa, alpha = x_1 * np.pi / 180, gamma = gamma * np.pi / 180)[1]
        F_y_2 = self.tire_eval(FZ = FZ, kappa = kappa, alpha = x_2 * np.pi / 180, gamma = gamma * np.pi / 180)[1]

        C_y = abs((F_y_2 - F_y_1) / (x_2 - x_1))

        return C_y
    
    def get_FNOMIN(self):
        return self._vertical_coeffs["FNOMIN"]

    def get_aligning_stiffness(self, FZ: float = 0, x_1: float = 0, x_2: float = 0.25, alpha: float = 0, kappa: float = 0, gamma: float = 0) -> float:
        if FZ == 0:
            FZ = self._vertical_coeffs["FNOMIN"]
        
        M_z_1 = self.tire_eval(FZ = FZ, kappa = kappa, alpha = x_1 * np.pi / 180, gamma = gamma * np.pi / 180)[5]
        M_z_2 = self.tire_eval(FZ = FZ, kappa = kappa, alpha = x_2 * np.pi / 180, gamma = gamma * np.pi / 180)[5]

        C_mz = abs((M_z_2 - M_z_1) / (x_2 - x_1))

        return C_mz

    def get_mu(self, FZ: float = 0, alpha: float = 0, kappa: float = 0, gamma: float = 0) -> list[float]:
        if FZ == 0:
            FZ = self._vertical_coeffs["FNOMIN"]
        
        # mu_x_result = 

        mu_y_result = F_y_result = get_F_y(
            lat_coeffs = self._lat_coeffs,
            scaling_coeffs = self.scaling_coeffs,
            vertical_coeffs = self._vertical_coeffs,
            dimensions = self._dimensions,
            operating_conditions = self._operating_conds, 
            FZ = FZ,
            alpha = alpha,
            kappa = kappa,
            gamma = gamma
        )[0]

        return [0, mu_y_result]

    def get_lateral_stiffness(self):
        return self.structural["LATERAL_STIFFNESS"]

    def plot_text(self):
        output_str = \
f"""
Properties at FNOMIN = {round(self._vertical_coeffs["FNOMIN"])} N

C_y (N/deg): {round(self.C_y)}
C_mz (Nm/deg): {round(self.C_mz, 3)}
mu_x (-): {round(self.mu_x, 3)}
mu_y (-): {round(self.mu_y, 3)}
K_y (N/m): {round(self.structural["LATERAL_STIFFNESS"])}
K_z (N/m): {round(self._vertical_coeffs["VERTICAL_STIFFNESS"])}
"""

        return output_str