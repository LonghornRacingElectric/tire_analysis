import matplotlib
import matplotlib.figure
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from .assets.tire_model.analysis import Analysis

def import_tires(tire_names: list[str], tir_paths: list[str]) -> Analysis:
    tires = Analysis(num = len(tire_names), name = tire_names, path = tir_paths)

    return tires

def save_pdf(figures: list[matplotlib.figure.Figure], file_name: str, save_location: str) -> None:
    p = PdfPages(f"{save_location}/{file_name}")

    for fig in figures:
        fig.savefig(p, format = "pdf")
    
    p.close()