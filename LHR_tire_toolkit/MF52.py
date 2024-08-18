import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D

from tire_model.file_processing._process_tir import _Processor

from tire_model.MF52_calculations._longitudinal_force import get_F_x
from tire_model.MF52_calculations._lateral_force import get_F_y
from tire_model.MF52_calculations._overturning_moment import get_M_x
from tire_model.MF52_calculations._rolling_resistance import get_M_y
from tire_model.MF52_calculations._aligning_moment import get_M_z


class MF52:
    """
    ## MF52

    Pacejka Magic Formula 5.2 tire model

    Parameters
    ----------
    tire_name : str
        Name of tire for internal use
    file_path : str
        File path to .tir file
    """
    def __init__(self, tire_name: str, file_path: str) -> None:
        self.tire_name = tire_name
        self._tire_params = _Processor(name = tire_name, file_path = file_path)

        self._init_coeffs()
        self._init_consts()

    def _init_coeffs(self) -> None:
        self._lat_coeffs = self._tire_params.get_parameters("LATERAL_COEFFICIENTS")
        self._long_coeffs = self._tire_params.get_parameters("LONGITUDINAL_COEFFICIENTS")
        self._overturning_coeffs = self._tire_params.get_parameters("OVERTURNING_COEFFICIENTS")
        self._rolling_coeffs = self._tire_params.get_parameters("ROLLING_COEFFICIENTS")
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
        """
        ## Tire Evaluation

        Evaluates tire forces and moments at a given FZ, alpha, kappa, and gamma

        Parameters
        ----------
        FZ : float
            Normal load in Newtons
        alpha : float
            Slip angle in radians
        kappa : float
            Slip ratio (unitless)
        gamma : float
            Inclination angle in radians
            
        Returns
        -------
        list[float]
            Tire forces and moments in the form: [F_x, F_y, F_z, M_x, M_y, M_z]
        """
        mu_x, F_x_result = get_F_x(
            long_coeffs = self._long_coeffs,
            scaling_coeffs = self.scaling_coeffs,
            vertical_coeffs = self._vertical_coeffs,
            dimensions = self._dimensions,
            operating_conditions = self._operating_conds, 
            FZ = FZ,
            alpha = alpha,
            kappa = kappa,
            gamma = gamma
        )

        mu_y, F_y_result = get_F_y(
            lat_coeffs = self._lat_coeffs,
            scaling_coeffs = self.scaling_coeffs,
            vertical_coeffs = self._vertical_coeffs,
            dimensions = self._dimensions,
            operating_conditions = self._operating_conds, 
            FZ = FZ,
            alpha = alpha,
            kappa = kappa,
            gamma = gamma
        )

        Mx_result = get_M_x(
            overturning_coeffs = self._overturning_coeffs,
            scaling_coeffs = self.scaling_coeffs,
            vertical_coeffs = self._vertical_coeffs,
            dimensions = self._dimensions,
            operating_conditions = self._operating_conds, 
            FZ = FZ,
            alpha = alpha,
            kappa = kappa,
            gamma = gamma,
            Fy = F_y_result
        )

        My_result = get_M_y(
            rolling_coeffs = self._rolling_coeffs,
            long_coeffs = self._long_coeffs,
            scaling_coeffs = self.scaling_coeffs,
            vertical_coeffs = self._vertical_coeffs,
            dimensions = self._dimensions,
            operating_conditions = self._operating_conds, 
            FZ = FZ,
            alpha = alpha,
            kappa = kappa,
            gamma = gamma
        )

        self.pneu_trail, MZ_result = get_M_z(
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
        )

        return [F_x_result, F_y_result, FZ, Mx_result, My_result, MZ_result]

    def get_cornering_stiffness(self, FZ: float = 0, x_1: float = 0, x_2: float = 0.25, alpha: float = 0, kappa: float = 0, gamma: float = 0) -> float:
        if FZ == 0:
            FZ = self._vertical_coeffs["FNOMIN"]
        
        F_y_1 = self.tire_eval(FZ = FZ, kappa = kappa, alpha = x_1 * np.pi / 180, gamma = gamma * np.pi / 180)[1]
        F_y_2 = self.tire_eval(FZ = FZ, kappa = kappa, alpha = x_2 * np.pi / 180, gamma = gamma * np.pi / 180)[1]

        C_y = abs((F_y_2 - F_y_1) / (x_2 - x_1))

        return C_y

    def get_camber_stiffness(self, FZ: float = 0, x_1: float = 0, x_2: float = 0.1, alpha: float = 0, kappa: float = 0, gamma: float = 0) -> float:
        if FZ == 0:
            FZ = self._vertical_coeffs["FNOMIN"]
        
        F_y_1 = self.tire_eval(FZ = FZ, kappa = kappa, alpha = 0 * np.pi / 180, gamma = x_1 * np.pi / 180)[1]
        F_y_2 = self.tire_eval(FZ = FZ, kappa = kappa, alpha = 0 * np.pi / 180, gamma = x_2 * np.pi / 180)[1]

        C_gamma = abs((F_y_2 - F_y_1) / (x_2 - x_1))

        return C_gamma

    def get_slip_stiffness(self, FZ: float = 0):
        if FZ == 0:
            FZ = self._vertical_coeffs["FNOMIN"]

        F_x_1 = self.tire_eval(FZ = FZ, kappa = -0.001, alpha = 0 * np.pi / 180, gamma = 0 * np.pi / 180)[0]
        F_x_2 = self.tire_eval(FZ = FZ, kappa = 0.001, alpha = 0 * np.pi / 180, gamma = 0 * np.pi / 180)[0]

        C_x = (-1 * (F_x_2 - F_x_1) / (0.001 - (-0.001)))

        return C_x

    def get_peak_F_y_alpha(self, FZ, approx = True):
        if approx:
            alpha_sweep = np.linspace(0, 20, 2000)

            FY = (self.tire_eval(FZ = FZ, alpha = alpha_sweep * np.pi / 180, kappa = 0, gamma = 0)[1])

            return list(alpha_sweep)[list(FY).index(min(FY))]

        else:
            alpha = minimize(lambda x: self.tire_eval(FZ = FZ, alpha = x[0] * np.pi / 180, kappa = 0, gamma = 0)[1], x0 = [7], bounds = [(0, 90)], method = "SLSQP").x

            return alpha[0]

    def get_peak_M_z_alpha(self, FZ, approx = True):
        if approx:
            alpha_sweep = np.linspace(0, 20, 2000)

            MZ = (self.tire_eval(FZ = FZ, alpha = alpha_sweep * np.pi / 180, kappa = 0, gamma = 0)[5])

            return list(alpha_sweep)[list(MZ).index(min(MZ))]

        else:
            alpha = minimize(lambda x: -1 * self.tire_eval(FZ = FZ, alpha = x[0] * np.pi / 180, kappa = 0, gamma = 0)[5], x0 = [2], bounds = [(0, 90)], method = "SLSQP").x

            return alpha[0]

    def get_peak_F_x_kappa(self, FZ, approx = True):
        if approx:
            kappa_sweep = np.linspace(0, 1, 2000)

            FX = (-1 * self.tire_eval(FZ = FZ, alpha = 0 * np.pi / 180, kappa = kappa_sweep, gamma = 0)[0])
            
            return list(kappa_sweep)[list(FX).index(min(FX))]

        else:
            kappa = minimize(lambda x: -1 * self.tire_eval(FZ = FZ, alpha = 0 * np.pi / 180, kappa = x[0], gamma = 0)[0], x0 = [2], bounds = [(0, 1)], method = "SLSQP").x

            return kappa[0]

    def get_FNOMIN(self):
        return self._vertical_coeffs["FNOMIN"]

    def get_pneu_trail(self, FZ: float = 0):
        if FZ == 0:
            FZ = self._vertical_coeffs["FNOMIN"]

        M_z_1 = self.tire_eval(FZ = FZ, kappa = 0, alpha = 0 * np.pi / 180, gamma = 0)[5]

        return self.pneu_trail

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
        
        mu_x_result, F_x_result = get_F_x(
            long_coeffs = self._long_coeffs,
            scaling_coeffs = self.scaling_coeffs,
            vertical_coeffs = self._vertical_coeffs,
            dimensions = self._dimensions,
            operating_conditions = self._operating_conds, 
            FZ = FZ,
            alpha = alpha,
            kappa = kappa,
            gamma = gamma
        )

        mu_y_result, F_y_result = get_F_y(
            lat_coeffs = self._lat_coeffs,
            scaling_coeffs = self.scaling_coeffs,
            vertical_coeffs = self._vertical_coeffs,
            dimensions = self._dimensions,
            operating_conditions = self._operating_conds, 
            FZ = FZ,
            alpha = alpha,
            kappa = kappa,
            gamma = gamma
        )

        return [mu_x_result, mu_y_result]

    def get_lateral_stiffness(self):
        return self.structural["LATERAL_STIFFNESS"]

    def get_vertical_stiffness(self):
        return self._vertical_coeffs["VERTICAL_STIFFNESS"]

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