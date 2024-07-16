import matplotlib.figure
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import time

from scipy.optimize import basinhopping
from scipy.optimize import minimize

from tire_toolkit.assets.tire_model.MF52 import MF52

class Analysis:
    def __init__(self, num: int, name: list[str], path: list[str]) -> None:
        self.num = num
        self.names = name
        self.paths = path

        self._tire_lst = []
        self._import_tires()

    def _import_tires(self) -> None:
        for i in range(len(self.names)):
            self._tire_lst.append(MF52(self.names[i], self.paths[i]))
    
    def get_loads(self, name: str, FZ: float, alpha: float, kappa: float, gamma: float) -> list[float]:
        selected_tire = self._tire_lst[self.names.index(name)]

        return selected_tire.tire_eval(FZ = FZ, alpha = alpha, kappa = kappa, gamma = gamma)

    def comparison_plot(self, tire_names: list[str], FZ_min: float, FZ_max: float) -> matplotlib.figure.Figure:

        selected_tires = []

        for name in tire_names:
            selected_tires.append(self._tire_lst[self.names.index(name)])

        FZ_sweep = np.linspace(FZ_min, FZ_max, 1000)

        fig, axes = plt.subplots(2, 3)
        fig2, axes2 = plt.subplots(2, 3)

        fig.subplots_adjust(left = 0.1, right = 0.9, bottom = 0.15, top = 0.85, wspace = 0.6, hspace = 0.5)
        fig2.subplots_adjust(left = 0.1, right = 0.9, bottom = 0.15, top = 0.85, wspace = 0.6, hspace = 0.5)

        fig.set_figwidth(12)
        fig.set_figheight(8)

        fig2.set_figwidth(12)
        fig2.set_figheight(8)

        ### Page 1

        fig.suptitle(f'First Principle Tire Metrics: Comparison')

        # Long Friction Coefficient
        axes[0, 0].set_title("Long Friction Coefficient vs FZ", pad = 10)
        axes[0, 0].set_xlabel("FZ (N)")
        axes[0, 0].set_ylabel("mu_x (-)")
        axes[0, 0].grid()

        # Lat Friction Coefficient
        axes[0, 1].set_title("Lat Friction Coefficient vs FZ", pad = 10)
        axes[0, 1].set_xlabel("FZ (N)")
        axes[0, 1].set_ylabel("mu_y (-)")
        axes[0, 1].grid()

        # Cornering Stiffness
        axes[0, 2].set_title("Cornering Stiffness vs FZ", pad = 10)
        axes[0, 2].set_xlabel("FZ (N)")
        axes[0, 2].set_ylabel("C_y (N/deg)")
        axes[0, 2].grid()

        # Aligning Moment Stiffness
        axes[1, 0].set_title("Aligning Moment Stiffness vs FZ", pad = 10)
        axes[1, 0].set_xlabel("FZ (N)")
        axes[1, 0].set_ylabel("C_mz (Nm/deg)")
        axes[1, 0].grid()

        # Relaxation Length
        axes[1, 1].set_title("Relaxation Length vs FZ", pad = 10)
        axes[1, 1].set_xlabel("FZ (N)")
        axes[1, 1].set_ylabel("Relaxation Length (m)")
        axes[1, 1].grid()

        # Pneumatic Trail
        axes[1, 2].set_title("Pneumatic Trail vs FZ", pad = 10)
        axes[1, 2].set_xlabel("FZ (N)")
        axes[1, 2].set_ylabel("Pneumatic Trail (mm)")
        axes[1, 2].grid()

        ### Page 2

        # Peak Fy Slip Angle
        axes2[0, 0].set_title("Peak Fy Alpha vs FZ", pad = 10)
        axes2[0, 0].set_xlabel("FZ (N)")
        axes2[0, 0].set_ylabel("Alpha (deg)")
        axes2[0, 0].grid()

        # Peak Mz Slip Angle
        axes2[0, 1].set_title("Peak Mz Alpha vs FZ", pad = 10)
        axes2[0, 1].set_xlabel("FZ (N)")
        axes2[0, 1].set_ylabel("Alpha (deg)")
        axes2[0, 1].grid()

        # Peak Fx Slip Ratio
        axes2[0, 2].set_title("Peak Fx Kappa vs FZ", pad = 10)
        axes2[0, 2].set_xlabel("FZ (N)")
        axes2[0, 2].set_ylabel("Kappa (-)")
        axes2[0, 2].grid()

        # Slip Stiffness
        axes2[1, 0].set_title("Slip Stiffness vs FZ", pad = 10)
        axes2[1, 0].set_xlabel("FZ (N)")
        axes2[1, 0].set_ylabel("Slip Stiffness (Fx/Kappa)")
        axes2[1, 0].grid()

        # Camber Stiffness
        axes2[1, 1].set_title("Camber Stiffness vs FZ", pad = 10)
        axes2[1, 1].set_xlabel("FZ (N)")
        axes2[1, 1].set_ylabel("Camber Stiffness (N/deg)")
        axes2[1, 1].grid()

        fig2.delaxes(axes2[1, 2])
        
        for tire in selected_tires:

            mu_x = []
            mu_y = []
            C_alpha = []
            C_kappa = []
            C_mz = []
            pneu_trail = []
            C_gamma = []
            rel_len = []
            peak_Fy_alpha = []
            peak_Mz_alpha = []
            peak_Fx_kappa = []

            for FZ in FZ_sweep:
                mu_x.append(tire.get_mu(FZ)[0])
                mu_y.append(tire.get_mu(FZ)[1])

                C_alpha.append(tire.get_cornering_stiffness(FZ, 0, 0.25))
                C_kappa.append(tire.get_slip_stiffness(FZ))
                C_mz.append(tire.get_aligning_stiffness(FZ, 0, 0.25))
                pneu_trail.append(tire.get_pneu_trail(FZ) * 1000)
                C_gamma.append(tire.get_camber_stiffness(FZ, 0, 0.25))
                rel_len.append(tire.get_cornering_stiffness(FZ, 0, 0.25) / tire.get_lateral_stiffness() * 180 / np.pi)
                peak_Fy_alpha.append(tire.get_peak_F_y_alpha(FZ, approx = False))
                peak_Mz_alpha.append(tire.get_peak_M_z_alpha(FZ, approx = False))
                peak_Fx_kappa.append(tire.get_peak_F_x_kappa(FZ, approx = False))

            ### Page 1

            # Long Friction Coefficient
            axes[0, 0].plot(FZ_sweep, mu_x)

            # Lat Friction Coefficient
            axes[0, 1].plot(FZ_sweep, mu_y)

            # Cornering Stiffness
            axes[0, 2].plot(FZ_sweep, C_alpha)

            # Aligning Moment Stiffness
            axes[1, 0].plot(FZ_sweep, C_mz)

            # Relaxation Length
            axes[1, 1].plot(FZ_sweep, rel_len)

            # Pneumatic Trail
            axes[1, 2].plot(FZ_sweep, pneu_trail)

            ### Page 2

            # Peak Fy Slip Angle
            axes2[0, 0].plot(FZ_sweep, peak_Fy_alpha)

            # Peak Mz Slip Angle
            axes2[0, 1].plot(FZ_sweep, peak_Mz_alpha)

            # Peak Fx Slip Ratio
            axes2[0, 2].plot(FZ_sweep, peak_Fx_kappa)

            # Slip Stiffness
            axes2[1, 0].plot(FZ_sweep, C_kappa)

            # Camber Stiffness
            axes2[1, 1].plot(FZ_sweep, C_gamma)

        fig.legend(tire_names, fontsize = "8", bbox_to_anchor = (1, 1))
        fig2.legend(tire_names, fontsize = "8", bbox_to_anchor = (1, 1))

        return [fig, fig2]

    def coeff_matching(self, reference_tire: str, target_tire: str, FZ_min: float, FZ_max: float, weighting: list[float], mesh: int) -> list:

        output_str = "Solving for Scaling Coefficients. View convergence progress below:"

        print("".join(["-" for x in list(output_str)]))

        print(f"\n{output_str}")

        self.reference_tire = self._tire_lst[self.names.index(reference_tire)]
        self.target_tire = self._tire_lst[self.names.index(target_tire)]

        original_coeffs = list(self.target_tire.scaling_coeffs.values()) + [self.reference_tire.structural["LATERAL_STIFFNESS"]]

        self.FZ_sweep = np.linspace(FZ_min, FZ_max, mesh)

        self.mesh = mesh
        self.previous_residual_norm = None
        
        self.resid_norm = np.inf
        self.current_coeffs = [1 for x in range(31)]
        self.weighting = weighting

        start = time.time()

        shift_factors = ["LHX", "LVX", "LHY", "LVY", "LXAL", "LVMX"]
        shape_factors = ["LCX", "LCY"]
        curvature_factors = ["LEX", "LEY"]
        bounds = []

        for key in self.target_tire.scaling_coeffs:
            if key in shift_factors:
                bounds.append((-5, 5))
            elif key in shape_factors:
                bounds.append((0.5, 2))
            elif key in curvature_factors:
                bounds.append((0.5, 2))
            else:
                bounds.append((0.2, 5))

        # gen_soln = basinhopping(
        #     func = self._full_eval,
        #     x0 = [1 for x in range(31)],
        #     minimizer_kwargs={"method": "SLSQP"}
        # ).x
        gen_soln = minimize(
            fun = self._full_eval,
            x0 = [1 for x in range(31)],
            bounds = bounds,
            method = "SLSQP",
            # method = "Nelder-Mead"
            ).x
        
        self.set_scaling_coeffs(self.target_tire.tire_name, gen_soln, output_tir = False)

        rel_len_soln = minimize(self._rel_len_eval, [1], bounds = [(0, 3)]).x

        self.set_scaling_coeffs(self.target_tire.tire_name, original_coeffs, output_tir = False)

        vert_stiff_soln = self.reference_tire._vertical_coeffs["VERTICAL_STIFFNESS"] / self.target_tire._vertical_coeffs["VERTICAL_STIFFNESS"]

        f = open("./outputs/current_coeff_soln.txt", "a")

        f.write(f"K_y_ratio = {rel_len_soln[0]}\nK_z_ratio = {vert_stiff_soln}")

        f.close()

        end = time.time()

        print("\n\nSolver Converged |", end = " ")

        print(f"Elapsted Time: {round(end - start)} seconds, Residual Norm: {self.resid_norm}\n")

        scaling_coeffs = [*gen_soln, *rel_len_soln, vert_stiff_soln]

        return scaling_coeffs

    def set_scaling_coeffs(self, tire_name: str, coeffs: list[float], output_tir: bool = True) -> None:
        self.target_tire = self._tire_lst[self.names.index(tire_name)]
        current_path = self.paths[self.names.index(tire_name)]

        index = 0
        for key in list(self.target_tire.scaling_coeffs.keys()):
            self.target_tire.scaling_coeffs[key] = coeffs[index]

            index += 1

        if len(coeffs) > 25:
            if coeffs[-1] < 100:
                self.target_tire.structural["LATERAL_STIFFNESS"] *= coeffs[-2]
                self.target_tire._vertical_coeffs["VERTICAL_STIFFNESS"] *= coeffs[-1]
            else:
                self.target_tire.structural["LATERAL_STIFFNESS"] = coeffs[-1]
            
        self._tire_lst[self.names.index(tire_name)] = self.target_tire

        if output_tir:
            self.save_tir(
                tire_names = 
                    [
                        tire_name
                    ],
                path_to_original_tirs = 
                    [
                        current_path
                    ],
                output_path = "./outputs"
            )
    
    def _full_eval(self, scaling_coeffs) -> float:

        index = 0
        for key in self.target_tire.scaling_coeffs:
            self.target_tire.scaling_coeffs[key] = scaling_coeffs[index]

            index += 1
        
        residuals = []

        for FZ in self.FZ_sweep:
            if self.weighting[0] != 0:
                ref_mu_x = self.reference_tire.get_mu(FZ)[0]
                iter_mu_x = self.target_tire.get_mu(FZ)[0]
                residuals.append((ref_mu_x - iter_mu_x) / ref_mu_x * 100 * self.weighting[0])

            if self.weighting[1] != 0:
                ref_mu_y = self.reference_tire.get_mu(FZ)[1]
                iter_mu_y = self.target_tire.get_mu(FZ)[1]
                residuals.append((ref_mu_y - iter_mu_y) / ref_mu_y * 100 * self.weighting[1])

            if self.weighting[2] != 0:
                ref_C_y = self.reference_tire.get_cornering_stiffness(FZ, 0, 0.25)
                iter_C_y = self.target_tire.get_cornering_stiffness(FZ, 0, 0.25)
                residuals.append((ref_C_y - iter_C_y) / ref_C_y * 100 * self.weighting[2])

            if self.weighting[3] != 0:
                ref_C_mz = self.reference_tire.get_aligning_stiffness(FZ, 0.001, 0.25)
                iter_C_mz = self.target_tire.get_aligning_stiffness(FZ, 0.001, 0.25)
                residuals.append((ref_C_mz - iter_C_mz) / ref_C_mz * 100 * self.weighting[3])

            if self.weighting[4] != 0:
                ref_peak_Fy_alpha = self.reference_tire.get_peak_F_y_alpha(FZ)
                iter_peak_Fy_alpha = self.target_tire.get_peak_F_y_alpha(FZ)
                residuals.append((ref_peak_Fy_alpha - iter_peak_Fy_alpha) / ref_peak_Fy_alpha * 100 * self.weighting[4])
            
            if self.weighting[5] != 0:
                ref_slip_stiff = self.reference_tire.get_slip_stiffness(FZ)
                iter_slip_stiff = self.target_tire.get_slip_stiffness(FZ)
                residuals.append((ref_slip_stiff - iter_slip_stiff) / ref_slip_stiff * 100 * self.weighting[5])

            if self.weighting[6] != 0:
                ref_peak_Mz_alpha = self.reference_tire.get_peak_M_z_alpha(FZ)
                iter_peak_Mz_alpha = self.target_tire.get_peak_M_z_alpha(FZ)
                residuals.append((ref_peak_Mz_alpha - iter_peak_Mz_alpha) / ref_peak_Mz_alpha * 100 * self.weighting[6])

            if self.weighting[7] != 0:
                ref_peak_Fx_kappa = self.reference_tire.get_peak_F_x_kappa(FZ)
                iter_peak_Fx_kappa = self.target_tire.get_peak_F_x_kappa(FZ)
                residuals.append((ref_peak_Fx_kappa - iter_peak_Fx_kappa) / ref_peak_Fx_kappa * 100 * self.weighting[7])

            if self.weighting[8] != 0:
                ref_C_gamma = self.reference_tire.get_camber_stiffness(FZ)
                iter_C_gamma = self.target_tire.get_camber_stiffness(FZ)
                residuals.append((ref_C_gamma - iter_C_gamma) / ref_C_gamma * 100 * self.weighting[8])
            
            if self.weighting[9] != 0:
                ref_pneu_trail = self.reference_tire.get_pneu_trail(FZ)
                iter_pneu_trail = self.target_tire.get_pneu_trail(FZ)
                residuals.append((ref_pneu_trail - iter_pneu_trail) / ref_pneu_trail * 100 * self.weighting[9])
        
        print(f"\rResidual Norm: {round(np.linalg.norm(residuals), 2)}          ", end = "")
        
        if np.linalg.norm(residuals) < self.resid_norm:
            self.current_coeffs = scaling_coeffs
            self.resid_norm = np.linalg.norm(residuals)

            f = open("./outputs/current_coeff_soln.txt", "w")

            index = 0
            for key in self.target_tire.scaling_coeffs:
                f.write(f"{key} = {scaling_coeffs[index]}\n")
                index += 1
            
            f.close()
        
        return np.linalg.norm(residuals)
    
    def _rel_len_eval(self, scale) -> float:

        residuals = []

        for FZ in self.FZ_sweep:
            lat_stiff = self.target_tire.structural["LATERAL_STIFFNESS"]

            lat_stiff *= scale[0]

            ref_rel_len = self.reference_tire.get_cornering_stiffness(FZ, 0, 0.25) / self.reference_tire.get_lateral_stiffness() * 180 / np.pi
            iter_rel_len = self.target_tire.get_cornering_stiffness(FZ, 0, 0.25) / lat_stiff * 180 / np.pi

            residuals.append((ref_rel_len - iter_rel_len) / ref_rel_len * 100)

        return np.linalg.norm(residuals)

    def save_tir(self, tire_names: list[str], path_to_original_tirs: list[str], output_path: str = "./outputs"):

        selected_tires = []

        for name in tire_names:
            selected_tires.append(self._tire_lst[self.names.index(name)])

        for tire in selected_tires:
            selected_tire = self._tire_lst[self.names.index(tire.tire_name)]

            f = open(path_to_original_tirs[tire_names.index(tire.tire_name)], "r")

            scaling_coeffs = selected_tire.scaling_coeffs
            structural_vals = selected_tire.structural
            vertical_vals = selected_tire._vertical_coeffs

            for line in f:
                if "=" in line:
                    line_stripped = line.replace(" ", "")

                    if "$" in line:
                        line_stripped = line_stripped[ : line_stripped.index("$")]
                    line_split = line_stripped.split("=")

                    if line_split[0] in scaling_coeffs.keys():
                        scaling_keys = list(scaling_coeffs.keys())
                        coeff_label = line_split[0]

                        coeff = list(str(round(float(scaling_coeffs[coeff_label]), 3)))
                        modified_line = list(line)
                        modified_line[modified_line.index("=") + 2 : ] = coeff

                        output_line = "".join(modified_line) + "\n"
                    
                    elif line_split[0] in structural_vals.keys():
                        structural_keys = list(structural_vals.keys())
                        coeff_label = line_split[0]

                        coeff = list(str(round(float(structural_vals[coeff_label]), 3)))
                        modified_line = list(line)
                        modified_line[modified_line.index("=") + 2 : ] = coeff

                        output_line = "".join(modified_line) + "\n"
                    
                    elif line_split[0] in vertical_vals.keys():
                        vertical_keys = list(vertical_vals.keys())
                        coeff_label = line_split[0]

                        coeff = list(str(round(float(vertical_vals[coeff_label]), 3)))
                        modified_line = list(line)
                        modified_line[modified_line.index("=") + 2 : ] = coeff

                        output_line = "".join(modified_line) + "\n"
                    
                    else:
                        output_line = line

                else:
                    output_line = line

                f = open(f"{output_path}/{tire.tire_name}.tir", "a")
                f.write(output_line)
                f.close()

            print(f".tir File Successfully Written to {output_path}/{tire.tire_name}.tir\n")