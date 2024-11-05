from tire_toolkit.tire_model.file_processing.process_tir import Processor
from scipy.optimize import minimize
from typing import Sequence, Tuple
import numpy as np

from tire_toolkit.tire_model.MF52_calculations.longitudinal_force import get_Fx
from tire_toolkit.tire_model.MF52_calculations.lateral_force import get_Fy
from tire_toolkit.tire_model.MF52_calculations.overturning_moment import get_Mx
from tire_toolkit.tire_model.MF52_calculations.rolling_resistance import get_My
from tire_toolkit.tire_model.MF52_calculations.aligning_moment import get_Mz


class MF52:
    """
    ## MF52 Tire Model

    Pacejka Magic Formula 5.2 Implementation

    Parameters
    ----------
    tire_name : str
        Internal name of tire for use in analysis
    file_path : str
        Path to desired .tir file
    """
    def __init__(self, tire_name: str, file_path: str) -> None:
        self.tire_name = tire_name
        self._tire_params = Processor(name = tire_name, file_path = file_path)

        self._init_coeffs()
        self._init_consts()

    def _init_coeffs(self) -> None:
        """
        ## Initialize Coefficients

        Stores coefficients locally from .tir file

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        self._lat_coeffs = self._tire_params.get_parameters("LATERAL_COEFFICIENTS")
        self._long_coeffs = self._tire_params.get_parameters("LONGITUDINAL_COEFFICIENTS")
        self._overturning_coeffs = self._tire_params.get_parameters("OVERTURNING_COEFFICIENTS")
        self._rolling_coeffs = self._tire_params.get_parameters("ROLLING_COEFFICIENTS")
        self._aligning_coeffs = self._tire_params.get_parameters("ALIGNING_COEFFICIENTS")

        self.scaling_coeffs = self._tire_params.get_parameters("SCALING_COEFFICIENTS")

    def _init_consts(self) -> None:
        """
        ## Initialize Constants

        Initializes general constants from .tir file

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        self._dimensions = self._tire_params.get_parameters("DIMENSION")
        self._operating_conds = self._tire_params.get_parameters("TYRE_CONDITIONS")
        self._vertical_coeffs = self._tire_params.get_parameters("VERTICAL")
        self.structural = self._tire_params.get_parameters("STRUCTURAL")
    
    def _param_eval(self, Fz: float) -> None:
        """
        ## Parameter Evaluation

        Evaluates tire parameters at a given normal load

        Parameters
        ----------
        Fz : float
            Normal load in Newtons
        """
        self.mu_x = self.get_mu(Fz)[0]
        self.mu_y = self.get_mu(Fz)[1]
        self.C_y = self.get_cornering_stiffness(Fz, 0, 0.25)
        self.C_mz = self.get_aligning_stiffness(Fz, 0, 0.25)

    def tire_eval(self, Fz: float, alpha: float, kappa: float, gamma: float) -> Sequence[float]:
        """
        ## Tire Evaluation

        Evaluates forces and moments for a given tire state

        Parameters
        ----------
        Fz : float
            Normal load in Newtons
        alpha : float
            Slip angle in radians
        kappa : float
            Dimensionless slip ratio
        gamma : float
            Inclination angle in radians

        Returns
        -------
        Sequence[float]
            Forces and moments in the form: [Fx, Fy, Fz, Mx, My, Mz]
        """
        mu_x, Fx_result = get_Fx(
            long_coeffs = self._long_coeffs,
            scaling_coeffs = self.scaling_coeffs,
            vertical_coeffs = self._vertical_coeffs,
            dimensions = self._dimensions,
            operating_conditions = self._operating_conds, 
            Fz = Fz,
            alpha = alpha,
            kappa = kappa,
            gamma = gamma
        )

        mu_y, Fy_result = get_Fy(
            lat_coeffs = self._lat_coeffs,
            scaling_coeffs = self.scaling_coeffs,
            vertical_coeffs = self._vertical_coeffs,
            dimensions = self._dimensions,
            operating_conditions = self._operating_conds, 
            Fz = Fz,
            alpha = alpha,
            kappa = kappa,
            gamma = gamma
        )

        Mx_result = get_Mx(
            overturning_coeffs = self._overturning_coeffs,
            scaling_coeffs = self.scaling_coeffs,
            vertical_coeffs = self._vertical_coeffs,
            dimensions = self._dimensions,
            operating_conditions = self._operating_conds, 
            Fz = Fz,
            alpha = alpha,
            kappa = kappa,
            gamma = gamma,
            Fy = Fy_result
        )

        My_result = get_My(
            rolling_coeffs = self._rolling_coeffs,
            long_coeffs = self._long_coeffs,
            scaling_coeffs = self.scaling_coeffs,
            vertical_coeffs = self._vertical_coeffs,
            dimensions = self._dimensions,
            operating_conditions = self._operating_conds, 
            Fz = Fz,
            alpha = alpha,
            kappa = kappa,
            gamma = gamma
        )

        self.pneu_trail, Mz_result = get_Mz(
            aligning_coeffs = self._aligning_coeffs,
            scaling_coeffs = self.scaling_coeffs,
            lat_coeffs = self._lat_coeffs,
            long_coeffs = self._long_coeffs,
            vertical_coeffs = self._vertical_coeffs,
            dimensions = self._dimensions,
            operating_conditions = self._operating_conds,
            Fz = Fz,
            alpha = alpha,
            kappa = kappa,
            gamma = gamma
        )

        return [Fx_result, Fy_result, Fz, Mx_result, My_result, Mz_result]

    def get_cornering_stiffness(self, Fz: float = 0, dx: float = 0.25, alpha: float = 0, kappa: float = 0, gamma: float = 0) -> float:
        """
        ## Get Cornering Stiffness

        Calculates cornering stiffness for given tire state

        Parameters
        ----------
        Fz : float, optional
            Normal load in Newtons, by default 0
        dx : float, optional
            Slip angle delta for derivative calculation, by default 0.25
        alpha : float, optional
            Slip angle in radians, by default 0
        kappa : float, optional
            Dimensionless slip ratio, by default 0
        gamma : float, optional
            Inclination angle in radians, by default 0

        Returns
        -------
        float
            Cornering stiffness in N/deg
        """
        if Fz == 0:
            Fz = self._vertical_coeffs["FNOMIN"]
        
        Fy_1 = self.tire_eval(Fz = Fz, kappa = kappa, alpha = alpha - dx / 2, gamma = gamma)[1]
        Fy_2 = self.tire_eval(Fz = Fz, kappa = kappa, alpha = alpha + dx / 2, gamma = gamma)[1]

        C_y = abs((Fy_2 - Fy_1) / dx)

        return C_y

    def get_camber_stiffness(self, Fz: float = 0, dx: float = 0.1, alpha: float = 0, kappa: float = 0, gamma: float = 0) -> float:
        """
        ## Get Camber Stiffness

        Calculates camber stiffness for given tire state

        Parameters
        ----------
        Fz : float, optional
            Normal load in Newtons, by default 0
        dx : float, optional
            Inclination angle delta for derivative calculation, by default 0.1
        alpha : float, optional
            Slip angle in radians, by default 0
        kappa : float, optional
            Dimensionless slip ratio, by default 0
        gamma : float, optional
            Inclination angle in radians, by default 0

        Returns
        -------
        float
            Cornering stiffness in N/deg
        """
        if Fz == 0:
            Fz = self._vertical_coeffs["FNOMIN"]
        
        F_y_1 = self.tire_eval(Fz = Fz, kappa = kappa, alpha = alpha, gamma = gamma - dx / 2)[1]
        F_y_2 = self.tire_eval(Fz = Fz, kappa = kappa, alpha = alpha, gamma = gamma + dx / 2)[1]

        C_gamma = abs((F_y_2 - F_y_1) / dx)

        return C_gamma

    def get_slip_stiffness(self, Fz: float = 0, dx: float = 0.002, alpha: float = 0, kappa: float = 0, gamma: float = 0) -> float:
        """
        ## Get Slip Stiffness

        Calculates slip stiffness for given tire state

        Parameters
        ----------
        Fz : float, optional
            Normal load in Newtons, by default 0
        dx : float, optional
            Dimensionless slip ratio delta for derivative calculation, by default 0.002
        alpha : float, optional
            Slip angle in radians, by default 0
        kappa : float, optional
            Dimensionless slip ratio, by default 0
        gamma : float, optional
            Inclination angle in radians, by default 0

        Returns
        -------
        float
            Cornering stiffness in N/deg
        """
        if Fz == 0:
            Fz = self._vertical_coeffs["FNOMIN"]

        F_x_1 = self.tire_eval(Fz = Fz, kappa = kappa - dx / 2, alpha = alpha, gamma = gamma)[0]
        F_x_2 = self.tire_eval(Fz = Fz, kappa = kappa + dx / 2, alpha = alpha, gamma = gamma)[0]

        C_x = -1 * (F_x_2 - F_x_1) / dx

        return C_x

    def get_F_y_at_gamma(self, Fz: float, approx: bool = True) -> Tuple[float, float]:
        """
        ## Get Peak Fy at Gamma

        Calculates the maximum Fy possible at a given normal load and inclination angle

        Parameters
        ----------
        Fz : float
            Normal load in Newtons
        approx : bool, optional
            Whether to solve for exact value (alternatively select nearest), by default True

        Returns
        -------
        Tuple[float, float]
            Tuple of the form: (gamma, Fy)
        """
        gamma_sweep = np.linspace(-10, 10, 25) * np.pi / 180
        alpha_sweep = np.linspace(0, 20, 25) * np.pi / 180
        max_Fy = 0
        max_Fy_gamma = 0

        
        if approx:
            for gamma in gamma_sweep:
                FY = (self.tire_eval(Fz = Fz, alpha = alpha_sweep, kappa = 0, gamma = gamma)[1])
                if max(abs(FY)) > max_Fy:
                    max_Fy = max(abs(FY))
                    max_Fy_gamma = gamma
        
        else:
            max_Fy_gamma = minimize(lambda x: -1 * max(abs(self.tire_eval(Fz = Fz, alpha = alpha_sweep, kappa = 0, gamma = x[0])[1])), x0 = [0], bounds = [(-10 * np.pi / 180, 10 * np.pi / 180)], method = "SLSQP").x
            peak_force = max(abs(self.tire_eval(Fz = Fz, alpha = alpha_sweep, kappa = 0, gamma = max_Fy_gamma)[1]))

        return max_Fy_gamma, peak_force

    def get_peak_F_y_alpha(self, Fz: float, approx: bool = True) -> float:
        """
        ## Get Alpha at Peak Fy

        Calculates the alpha corresponding to the peak Fy at a given Fz

        Parameters
        ----------
        Fz : float
            Normal load in Newtons
        approx : bool, optional
            Whether to solve for exact value (alternatively select nearest), by default True

        Returns
        -------
        float
            Slip angle in radians
        """
        if approx:
            alpha_sweep = np.linspace(0, 20, 2000)

            FY = (self.tire_eval(Fz = Fz, alpha = alpha_sweep, kappa = 0, gamma = 0)[1])

            return list(alpha_sweep)[list(FY).index(min(FY))]

        else:
            alpha = minimize(lambda x: self.tire_eval(Fz = Fz, alpha = x[0], kappa = 0, gamma = 0)[1], x0 = [7], bounds = [(0, 90)], method = "SLSQP").x

            return alpha[0]

    def get_peak_M_z_alpha(self, Fz: float, approx: bool = True) -> float:
        """
        ## Get Alpha at Peak Mz

        Calculates the alpha corresponding to the peak Mz at a given Fz

        Parameters
        ----------
        Fz : float
            Normal load in Newtons
        approx : bool, optional
            Whether to solve for exact value (alternatively select nearest), by default True

        Returns
        -------
        float
            Slip angle in radians
        """
        if approx:
            alpha_sweep = np.linspace(0, 20, 2000)

            MZ = (self.tire_eval(Fz = Fz, alpha = alpha_sweep, kappa = 0, gamma = 0)[5])

            return list(alpha_sweep)[list(MZ).index(min(MZ))]

        else:
            alpha = minimize(lambda x: -1 * self.tire_eval(Fz = Fz, alpha = x[0], kappa = 0, gamma = 0)[5], x0 = [2 * np.pi / 180], bounds = [(0, 90 * np.pi / 180)], method = "SLSQP").x

            return alpha[0]

    def get_peak_F_x_kappa(self, Fz: float, approx: bool = True) -> float:
        """
        ## Get Kappa at Peak Fx

        Calculates the kappa corresponding to the peak Fx at a given Fz

        Parameters
        ----------
        Fz : float
            Normal load in Newtons
        approx : bool, optional
            Whether to solve for exact value (alternatively select nearest), by default True

        Returns
        -------
        float
            Dimensionless slip ratio
        """
        if approx:
            kappa_sweep = np.linspace(0, 1, 2000)

            FX = (-1 * self.tire_eval(Fz = Fz, alpha = 0, kappa = kappa_sweep, gamma = 0)[0])
            
            return list(kappa_sweep)[list(FX).index(min(FX))]

        else:
            kappa = minimize(lambda x: -1 * self.tire_eval(Fz = Fz, alpha = 0, kappa = x[0], gamma = 0)[0], x0 = [2], bounds = [(0, 1)], method = "SLSQP").x

            return kappa[0]

    def get_pneu_trail(self, Fz: float = 0) -> float:
        """
        ## Get Pneumatic Trail

        Calculates pneumatic trail at a given normal load

        Parameters
        ----------
        Fz : float, optional
            Normal load in Newtons, by default 0

        Returns
        -------
        float
            Pneumatic trail in meters
        """
        if Fz == 0:
            Fz = self._vertical_coeffs["FNOMIN"]

        M_z_1 = self.tire_eval(Fz = Fz, kappa = 0, alpha = 0, gamma = 0)[5]

        return self.pneu_trail

    def get_aligning_stiffness(self, Fz: float = 0, dx: float = 0.25, alpha: float = 0, kappa: float = 0, gamma: float = 0) -> float:
        """
        ## Get Aligning Stiffness

        Calculates aligning stiffness for given tire state

        Parameters
        ----------
        Fz : float, optional
            Normal load in Newtons, by default 0
        dx : float, optional
            Slip angle delta for derivative calculation, by default 0.25
        alpha : float, optional
            Slip angle in radians, by default 0
        kappa : float, optional
            Dimensionless slip ratio, by default 0
        gamma : float, optional
            Inclination angle in radians, by default 0

        Returns
        -------
        float
            Aligning stiffness in N/deg
        """
        if Fz == 0:
            Fz = self._vertical_coeffs["FNOMIN"]
        
        M_z_1 = self.tire_eval(Fz = Fz, kappa = kappa, alpha = alpha - dx / 2, gamma = gamma)[5]
        M_z_2 = self.tire_eval(Fz = Fz, kappa = kappa, alpha = alpha + dx / 2, gamma = gamma)[5]

        C_mz = abs((M_z_2 - M_z_1) / (dx))

        return C_mz

    def get_mu(self, Fz: float = 0, alpha: float = 0, kappa: float = 0, gamma: float = 0) -> Sequence[float]:
        """
        ## Get Mu

        Calculates lateral and longitudinal friction coefficients for given tire state

        Parameters
        ----------
        Fz : float, optional
            Normal load in Newtons, by default 0
        alpha : float, optional
            Slip angle in radians, by default 0
        kappa : float, optional
            Dimensionless slip ratio, by default 0
        gamma : float, optional
            Inclination angle in radians, by default 0

        Returns
        -------
        Sequence[float]
            List of friction coefficients in the form: [mu_x, mu_y]
        """
        if Fz == 0:
            Fz = self._vertical_coeffs["FNOMIN"]
        
        mu_x_result, F_x_result = get_Fx(
            long_coeffs = self._long_coeffs,
            scaling_coeffs = self.scaling_coeffs,
            vertical_coeffs = self._vertical_coeffs,
            dimensions = self._dimensions,
            operating_conditions = self._operating_conds, 
            Fz = Fz,
            alpha = alpha,
            kappa = kappa,
            gamma = gamma
        )

        mu_y_result, F_y_result = get_Fy(
            lat_coeffs = self._lat_coeffs,
            scaling_coeffs = self.scaling_coeffs,
            vertical_coeffs = self._vertical_coeffs,
            dimensions = self._dimensions,
            operating_conditions = self._operating_conds, 
            Fz = Fz,
            alpha = alpha,
            kappa = kappa,
            gamma = gamma
        )

        return [mu_x_result, mu_y_result]
    
    @property
    def Fz(self):
        return self._vertical_coeffs["FNOMIN"]

    @property
    def lateral_stiffness(self):
        return self.structural["LATERAL_STIFFNESS"]

    @property
    def vertical_stiffness(self):
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