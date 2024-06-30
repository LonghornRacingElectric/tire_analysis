import matplotlib.figure
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd

from .MF52 import MF52

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

    # TODO: set default FZ_min and FZ_max as values from .tir
    def single_plot(self, name: str, FZ_min: float,FZ_max: float) -> matplotlib.figure.Figure:
        selected_tire = self._tire_lst[self.names.index(name)]

        FZ_sweep = np.linspace(FZ_min, FZ_max, 1000)

        mu_x = []
        mu_y = []
        C_y = []
        C_mz = []
        rel_len = []

        for FZ in FZ_sweep:
            mu_x.append(selected_tire.get_mu(FZ)[0])
            mu_y.append(selected_tire.get_mu(FZ)[1])

            C_y.append(selected_tire.get_cornering_stiffness(FZ, 0, 0.25))
            C_mz.append(selected_tire.get_aligning_stiffness(FZ, 0, 0.25))

            rel_len.append(selected_tire.get_cornering_stiffness(FZ, 0, 0.25) / selected_tire.get_lateral_stiffness() * 180 / (np.pi))

        fig, axes = plt.subplots(2, 3)

        plt.subplots_adjust(left = 0.1, right = 0.9, bottom = 0.15, top = 0.85, wspace = 0.6, hspace = 0.5)

        fig.set_figwidth(12)
        fig.set_figheight(8)

        fig.suptitle(f'First Principle Tire Metrics: {name}')

        # Long Friction Coefficient
        axes[0, 0].plot(FZ_sweep, mu_x)
        axes[0, 0].set_title("Long Friction Coefficient vs FZ", pad = 20)
        axes[0, 0].set_xlabel("FZ (N)")
        axes[0, 0].set_ylabel("mu_x (-)")
        axes[0, 0].grid()

        # Lat Friction Coefficient
        axes[0, 1].plot(FZ_sweep, mu_y)
        axes[0, 1].set_title("Lat Friction Coefficient vs FZ", pad = 20)
        axes[0, 1].set_xlabel("FZ (N)")
        axes[0, 1].set_ylabel("mu_y (-)")
        axes[0, 1].grid()

        # Cornering Stiffness
        axes[0, 2].plot(FZ_sweep, C_y)
        axes[0, 2].set_title("Cornering Stiffness vs FZ", pad = 20)
        axes[0, 2].set_xlabel("FZ (N)")
        axes[0, 2].set_ylabel("C_y (N/deg)")
        axes[0, 2].grid()

        # Aligning Moment Stiffness
        axes[1, 0].plot(FZ_sweep, C_mz)
        axes[1, 0].set_title("Aligning Moment Stiffness vs FZ", pad = 20)
        axes[1, 0].set_xlabel("FZ (N)")
        axes[1, 0].set_ylabel("C_mz (Nm/deg)")
        axes[1, 0].grid()

        # Relaxation Length
        axes[1, 1].plot(FZ_sweep, rel_len)
        axes[1, 1].set_title("Relaxation Length vs FZ", pad = 20)
        axes[1, 1].set_xlabel("FZ (N)")
        axes[1, 1].set_ylabel("Relaxation Length (m)")
        axes[1, 1].grid()

        # Remove sixth, unused plot
        fig.delaxes(axes[1, 2])
        fig.text(0.7, 0.15, selected_tire.plot_text(), fontsize = 12, bbox = dict(boxstyle = 'round, pad=1', facecolor = 'wheat', alpha = 1))

        return fig

    def comparison_plot(self, tire_names: list[str], FZ_min: float, FZ_max: float) -> matplotlib.figure.Figure:

        selected_tires = []

        for name in tire_names:
            selected_tires.append(self._tire_lst[self.names.index(name)])

        FZ_sweep = np.linspace(FZ_min, FZ_max, 1000)

        fig, axes = plt.subplots(2, 3)

        plt.subplots_adjust(left = 0.1, right = 0.9, bottom = 0.15, top = 0.85, wspace = 0.6, hspace = 0.5)

        fig.set_figwidth(12)
        fig.set_figheight(8)

        fig.suptitle(f'First Principle Tire Metrics: Comparison')

        # Long Friction Coefficient
        axes[0, 0].set_title("Long Friction Coefficient vs FZ", pad = 20)
        axes[0, 0].set_xlabel("FZ (N)")
        axes[0, 0].set_ylabel("mu_x (-)")
        axes[0, 0].grid()

        # Lat Friction Coefficient
        axes[0, 1].set_title("Lat Friction Coefficient vs FZ", pad = 20)
        axes[0, 1].set_xlabel("FZ (N)")
        axes[0, 1].set_ylabel("mu_y (-)")
        axes[0, 1].grid()

        # Cornering Stiffness
        axes[0, 2].set_title("Cornering Stiffness vs FZ", pad = 20)
        axes[0, 2].set_xlabel("FZ (N)")
        axes[0, 2].set_ylabel("C_y (N/deg)")
        axes[0, 2].grid()

        # Aligning Moment Stiffness
        axes[1, 0].set_title("Aligning Moment Stiffness vs FZ", pad = 20)
        axes[1, 0].set_xlabel("FZ (N)")
        axes[1, 0].set_ylabel("C_mz (Nm/deg)")
        axes[1, 0].grid()

        # Relaxation Length
        axes[1, 1].set_title("Relaxation Length vs FZ", pad = 20)
        axes[1, 1].set_xlabel("FZ (N)")
        axes[1, 1].set_ylabel("Relaxation Length (m)")
        axes[1, 1].grid()

        # Remove sixth, unused plot
        fig.delaxes(axes[1, 2])
        
        for tire in selected_tires:

            mu_x = []
            mu_y = []
            C_y = []
            C_mz = []
            rel_len = []

            for FZ in FZ_sweep:
                mu_x.append(tire.get_mu(FZ)[0])
                mu_y.append(tire.get_mu(FZ)[1])

                C_y.append(tire.get_cornering_stiffness(FZ, 0, 0.25))
                C_mz.append(tire.get_aligning_stiffness(FZ, 0, 0.25))

                rel_len.append(tire.get_cornering_stiffness(FZ, 0, 0.25) / tire.get_lateral_stiffness() * 180 / (np.pi))

            # Long Friction Coefficient
            axes[0, 0].plot(FZ_sweep, mu_x)

            # Lat Friction Coefficient
            axes[0, 1].plot(FZ_sweep, mu_y)

            # Cornering Stiffness
            axes[0, 2].plot(FZ_sweep, C_y)

            # Aligning Moment Stiffness
            axes[1, 0].plot(FZ_sweep, C_mz)

            # Relaxation Length
            axes[1, 1].plot(FZ_sweep, rel_len)

        fig.legend(tire_names, fontsize = "14", bbox_to_anchor = (0.875, 0.35))

        return fig