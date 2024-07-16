from tire_toolkit import workflows
from tire_toolkit.assets.tire_model.analysis import Analysis

class ProcessTire:
    """
    Description
    -----------
    MF5.2 Tire Processing Object

    Parameters
    ----------
    tire_name : str
        Name of tire
    tir_path : str
        Path to .tir file
    
    Returns
    -------
    None

    """
    
    def __init__(self, tire_name: str, tir_path: str) -> None:
        self.tire_name: str = tire_name
        self.tir_path: str = tir_path

        self.tires: Analysis = workflows.import_tires(
            tire_names = [tire_name],
            tir_paths = [tir_path])
    
    def get_loads(self, FZ: float=0, alpha: float=0, kappa: float=0, gamma: float=0):
        """
        Description
        -----------
        Calculates Tire Loads from MF5.2 Model

        Parameters
        ----------
        FZ : float
            Normal Load (N)
        alpha : float
            Slip Angle (rad)
        kappa : float
            Slip Ratio (-)
        gamma : float
            Inclination Angle (rad)

        Returns
        -------
        list[float]
            Tire Loads [Fx, Fy, Fz, Mx, My, Mz]

        """
        
        return [float(x) for x in self.tires.get_loads(name=self.tire_name, FZ=FZ, alpha=alpha * -1, kappa=kappa, gamma=gamma)]

