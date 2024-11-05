from tire_toolkit.tire_model.analysis import Analysis
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.figure import Figure
from typing import Sequence


def import_tires(tire_names: Sequence[str], tir_paths: Sequence[str]) -> Analysis:
    """
    ## Import Tires

    Imports specified tires

    Parameters
    ----------
    tire_names : Sequence[str]
        Internal name for use of tires in analysis
    tir_paths : Sequence[str]
        File path to corresponding .tir files in tire_names

    Returns
    -------
    Analysis
        Analysis object for qunaitfying tire metrics
    """
    tires = Analysis(num = len(tire_names), name = tire_names, path = tir_paths)

    return tires

def plot_tires(tire_object: Analysis, names: Sequence[str], FZ_min: float, FZ_max: float, file_name: str) -> None:
    """
    ## Plot Tires

    Generates PDF report of tire metrics 

    Parameters
    ----------
    tire_object : Analysis
        Tire analysis object for comparisons
    names : Sequence[str]
        Internal names of desired tires to plot
    FZ_min : float
        Minimum Fz value for sweeps
    FZ_max : float
        Maximum Fz value for sweeps
    file_name : str
        Output name of PDF file

    Returns
    -------
    None
    """
    figs = tire_object.comparison_plot(tire_names = names, FZ_min = FZ_min, FZ_max = FZ_max)

    save_pdf(figures = figs, file_name = file_name, save_location = "./outputs")

def save_pdf(figures: Sequence[Figure], file_name: str, save_location: str) -> None:
    """
    ## Save PDF

    Saves desired figuers to PDF (each figure passed becomes a separate page)

    Parameters
    ----------
    figures : Sequence[Figure]
        Sequence of figures to save to PDF
    file_name : str
        Output name of PDF file
    save_location : str
        Location to save output file
    
    Returns
    -------
    None
    """
    p = PdfPages(f"{save_location}/{file_name}")

    for fig in figures:
        fig.savefig(p, format = "pdf")
    
    p.close()