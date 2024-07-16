import matplotlib
import matplotlib.figure
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from tire_model.analysis import Analysis

def import_tires(tire_names: list[str], tir_paths: list[str]) -> Analysis:
    tires = Analysis(num = len(tire_names), name = tire_names, path = tir_paths)

    return tires

def plot_tires(tire_object: Analysis, names: list[str], FZ_min: float, FZ_max: float, file_name: str) -> None:
    # if len(names) == 1:
    #     fig = tire_object.single_plot(name = names[0], FZ_min = FZ_min, FZ_max = FZ_max)
    
    # else:
    #     fig = tire_object.comparison_plot(tire_names = names, FZ_min = FZ_min, FZ_max = FZ_max)

    figs = tire_object.comparison_plot(tire_names = names, FZ_min = FZ_min, FZ_max = FZ_max)

    save_pdf(figures = figs, file_name = file_name, save_location = "./outputs")

def save_pdf(figures: list[matplotlib.figure.Figure], file_name: str, save_location: str) -> None:
    p = PdfPages(f"{save_location}/{file_name}")

    for fig in figures:
        fig.savefig(p, format = "pdf")
    
    p.close()